%% prototype_single_water_FK.m
% Minimal 1D water-column prototype:
% - Diffusion (axial)
% - Forced convection (axial advection, upwind)
% - Optional free convection module (your FK logic adapted to single volume)
% - Injection at internal stub (z_in = 0.7H) as local energy source term
%
% Author: (you)
% Notes:
% - Robin-type heat exchange to concrete/soil is approximated as a distributed sink term.

clear; clc

%% ----------------------------- Geometry ------------------------------
H   = 20.0;          % [m] water height
dz  = 0.05;          % [m] grid size
N   = round(H/dz) + 1;
z   = linspace(0,H,N)';  % [m], z=0 bottom, z=H top

A   = 200.0;         % [m^2] cross-sectional area

% Injection (stub) position
z_in       = 0.7 * H;
[~, idx_in] = min(abs(z - z_in));

%% ----------------------------- Physics -------------------------------
rho = 997;           % [kg/m^3]
cp  = 4180;          % [J/(kg*K)]
k_w = 0.6;           % [W/(m*K)]
alpha = k_w/(rho*cp);% [m^2/s] thermal diffusivity

% Forced convection velocity (sign sets direction)
% Convention: v > 0 means upward (towards increasing z), v < 0 downward.
v = +2e-5;           % [m/s] (change sign to test direction)
Vdot = abs(v)*A;     % [m^3/s]
mdot = rho*Vdot;     % [kg/s]

% Free convection time scale (Schäfer fixed interval)
dt_mix = 900;        % [s]

% Ambient piston/soil/concrete temperature
T_env = 45.0;        % [°C]

% Side-wall heat loss: approximate Robin via UA per length.
% Equivalent h through concrete: h_eq = k_conc / t_conc (conduction-limited)
k_conc = 1.7;        % [W/(m*K)] typical concrete
t_conc = 0.30;       % [m] assumed wall thickness
h_eq   = k_conc / max(t_conc, eps);   % [W/(m^2*K)]

% Need a perimeter to convert h*A_side to UA; with only area given, assume circular tank:
% A = pi*R^2 => R = sqrt(A/pi), perimeter P = 2*pi*R, side area per length = P
R = sqrt(A/pi);
P = 2*pi*R;          % [m] perimeter
UA_per_length = h_eq * P;  % [W/(m*K)] distributed UA
k_loss = UA_per_length / (rho*cp*A);  % [1/s] coefficient for sink term: -k_loss*(T-T_env)

%% ----------------------------- Numerics ------------------------------
use_free_convection = true;   % toggle FK on/off

dt    = 30;                   % [s] internal timestep (stable for diffusion/advection here)
t_end = 6*3600;               % [s] simulate a few hours
Nt    = round(t_end/dt);

% We apply FK term each step, but it internally uses dt_mix (as in your module).
% For a strict interpretation you could apply FK only every 900 s; see flag below.
apply_FK_only_on_dtMix = true;

%% ----------------------------- Initial profiles -----------------------
% Profile 1: smooth increasing from 45°C to 75°C
T0_ramp = 45 + (75-45)*(z/H);

% Profile 2: constant 60°C
T0_flat = 60*ones(N,1);

profiles = {T0_ramp, T0_flat};
profile_names = {'Ramp 45->75', 'Flat 60'};

% Injection temperatures for each profile: warmer 70 and colder 50
Tin_list = [70, 50];  % [°C]

%% ----------------------------- Run cases ------------------------------
for p = 1:numel(profiles)
    for tinCase = 1:numel(Tin_list)

        Tin = Tin_list(tinCase);
        T = profiles{p};

        % Store snapshots every 900s for plotting
        snap_every = 900;                  % [s]
        snap_stride = round(snap_every/dt);
        snaps_T = [];
        snaps_t = [];

        for n = 1:Nt
            t = (n-1)*dt;

            % --- core PDE terms: diffusion + forced convection + side loss + injection
            dTdt = zeros(N,1);

            % Diffusion: alpha * d2T/dz2 (central)
            d2 = second_derivative_central(T, dz);

            % Advection: -v * dT/dz (upwind)
            d1 = first_derivative_upwind(T, dz, v, T_env);

            % Side losses: -k_loss*(T - T_env)
            loss = -k_loss*(T - T_env);

            % Injection as local source (well-mixed control volume at idx_in)
            inj = zeros(N,1);
            inj(idx_in) = (mdot/(rho*A*dz)) * (Tin - T(idx_in)); % [K/s]

            dTdt = alpha*d2 - v*d1 + loss + inj;

            % --- optional free convection contribution
            if use_free_convection
                doFK = true;
                if apply_FK_only_on_dtMix
                    % Apply FK only when t hits multiples of dt_mix
                    doFK = (mod(t, dt_mix) < 0.5*dt); % tolerant trigger
                end

                if doFK
                    dTdt_fk = free_convection_single_volume( ...
                        T, Tin, idx_in, dz, dt_mix, mdot, cp, A, rho);
                    dTdt = dTdt + dTdt_fk;
                end
            end

            % --- time integration
            T = T + dt*dTdt;

            % snapshots
            if mod(n-1, snap_stride) == 0
                snaps_T = [snaps_T, T]; %#ok<AGROW>
                snaps_t = [snaps_t; t]; %#ok<AGROW>
            end
        end

        % ----------------------------- Plot -----------------------------
        fig = figure('Name', sprintf('%s | Tin=%g | FK=%d', profile_names{p}, Tin, use_free_convection));
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

        nexttile;
        plot(profiles{p}, z, '--', 'LineWidth', 1.5); hold on;
        plot(snaps_T(:,end), z, '-', 'LineWidth', 2.0);
        grid on;
        xlabel('T [°C]'); ylabel('z [m]');
        title(sprintf('%s | Tin=%g°C | v=%g m/s', profile_names{p}, Tin, v), 'Interpreter','none');
        legend('initial','final','Location','best');

        nexttile;
        plot(snaps_T(:,end) - profiles{p}, z, '-', 'LineWidth', 2.0);
        grid on;
        xlabel('\DeltaT final-initial [K]'); ylabel('z [m]');
        title(sprintf('Delta (FK=%d, applyFK@%ds=%d)', use_free_convection, dt_mix, apply_FK_only_on_dtMix));
        xline(0,'--');

    end
end

%% ============================ Functions ===============================

function d2 = second_derivative_central(T, dz)
% Compute second derivative with simple boundary handling (Neumann 0 at ends).
% Comment: For prototype, we keep axial ends insulated for diffusion term.
    N = numel(T);
    d2 = zeros(N,1);

    % interior
    d2(2:N-1) = (T(3:N) - 2*T(2:N-1) + T(1:N-2)) / dz^2;

    % Neumann(0): dT/dz=0 => mirror
    d2(1)   = (T(2)   - T(1)) / dz^2;
    d2(N)   = (T(N-1) - T(N)) / dz^2;
end

function d1 = first_derivative_upwind(T, dz, v, T_bc)
% Upwind first derivative for advection term.
% Uses a Dirichlet ghost value T_bc at the inflow boundary.
    N = numel(T);
    d1 = zeros(N,1);

    if v > 0
        % flow upward: inflow at bottom (z=0)
        Tghost = T_bc;
        d1(1) = (T(1) - Tghost) / dz;
        d1(2:N) = (T(2:N) - T(1:N-1)) / dz;
    elseif v < 0
        % flow downward: inflow at top (z=H)
        Tghost = T_bc;
        d1(N) = (Tghost - T(N)) / dz;   % consistent with backward difference direction
        d1(1:N-1) = (T(2:N) - T(1:N-1)) / dz;
    else
        d1(:) = 0;
    end
end

function dTdt_free = free_convection_single_volume( ...
    T, Tin, idx_in, dz, dt_mix, mdot, c_w, A, rho_w)
% Free convection module adapted to single water volume with an internal injection point.
%
% Rules:
% - If Tin > T(idx_in): buoyant hot water rises -> mix upwards (increasing z)
% - If Tin < T(idx_in): colder water sinks -> mix downwards (decreasing z)
%
% This follows your structure: determine mixing zone, compute I_mix, build geom_factor, add Qdot term.

    N = numel(T);
    dTdt_free = zeros(N,1);

    T_in_local = Tin;
    T_at_in    = T(idx_in);

    % Decide direction
    if T_in_local > T_at_in
        % upward mixing
        z_mix_idx = idx_in;
        while (T_in_local > T(z_mix_idx) && z_mix_idx < N)
            z_mix_idx = z_mix_idx + 1;
        end
        ids_mix = idx_in:z_mix_idx;
        T_mix   = T(ids_mix);
        I_mix   = sum((T_in_local - T_mix) * dz);  % [K*m]
    else
        % downward mixing
        z_mix_idx = idx_in;
        while (T_in_local < T(z_mix_idx) && z_mix_idx > 1)
            z_mix_idx = z_mix_idx - 1;
        end
        ids_mix = idx_in:-1:z_mix_idx;
        T_mix   = T(ids_mix);
        I_mix   = sum((T_in_local - T_mix) * dz * (-1)); % [K*m]
    end

    % Energy flow rate associated with injection temperature difference at inlet cell
    Qdot = mdot * c_w * (T_in_local - T_at_in); % [W]

    dz_mix = max((numel(ids_mix)-1) * dz, eps); % [m]
    z_dist = (ids_mix - z_mix_idx).' * dz;       % [m]
    geom_factor = 2 * z_dist / (dz_mix^2);       % [1/m]

    dTdt_free(ids_mix) = (1/dt_mix) * (geom_factor * I_mix + (T_in_local - T_mix)) ...
                       - Qdot * geom_factor / (A * rho_w * c_w);
end
