%% Free convection (EXTRAPOLATED zmix) + optional diffusion toy example
clear; clc; close all;

%%% -----------------------------
%  Toggle diffusion on/off
%%% -----------------------------
enableDiffusion = false;   % <-- switch here

%% Spatial grid
z  = linspace(0,20,4000)';          % [m]
dz = z(2) - z(1);                   % [m]

%% Initial temperature profile
T0 = 50 + z;                        % linear reference profile

%% Dirichlet boundary conditions
T_left  = T0(1);
T_right = T0(end);

%% -----------------------------
%  Water material parameters (same as original script)
%% -----------------------------
rho_w    = 997.6;                   % [kg/m^3]
cp_w     = 4.177;                   % [kJ/(kg*K)]
lambda_w = 0.6081;                  % [W/(m*K)]
alpha    = lambda_w / (rho_w * cp_w);

%% Time steps to visualize
dt_list = [0,1,2,3,4,5,6]*7200;      % [s]

%% Plot
figure; tiledlayout(1,2,'Padding','compact'); grid on;

%% ==============================================================
%  CASE A: Tin BELOW minimum temperature  (cold inflow)
% ==============================================================
nexttile; hold on; grid on;
title('Extrapolated FK: Tin < Tmin(Tmix)');

for dt = dt_list

    T = T0;

    % ------------------------------------------------------------
    % (1) Extrapolated free convection operator
    % ------------------------------------------------------------

    ids = (z >= 0 & z <= 10);        % downward mixing region
    z_in = 10;                       % [m]
    z_mix = 0;                       % [m]

    Tmix = T(ids);
    zloc = z(ids);

    Tin = min(Tmix) - 10;            % <<< extrapolated COLD inlet

    % Integral over actual mixing zone
    I_num = abs(sum((Tin - Tmix)) * dz);

    % Mixing power (same toy value as original)
    Qdot_mix = 577;                 % [W]
    A_hws    = 300;                  % [m^2]
    rho      = rho_w;
    cp       = cp_w * 1e3;           % [J/(kg*K)]

    dt_mix = dt;

    % Energy term
    E = abs(Qdot_mix * dt_mix) / (rho * cp * A_hws);

    % Linear regression slope over physical zone
    mix_height = abs(z_mix - z_in);
    dTdz = (Tmix(end) - Tmix(1)) / mix_height;

    % Virtual mixing height
    dz_star = abs((Tin - Tmix(end)) / dTdz);

    % Geometric integral
    z0 = 0;
    z1 = mix_height;
    z1s = z1 + dz_star;

    G = 0.5*(z1^2 - z0^2) - z1s*(z1 - z0);

    % New slope
    a = (E - I_num) / G;

    % Apply extrapolated FK update
    z_phys = (zloc - z_mix);
    T(ids) = Tin + a * (z_phys - dz_star);

    % ------------------------------------------------------------
    % (2) Optional diffusion
    % ------------------------------------------------------------
    if enableDiffusion
        T = diffuse_dirichlet_explicit(T, alpha, dz, dt, T_left, T_right);
    end

    plot(T, z, 'LineWidth', 1.5);
end

plot(T0, z, 'k--', 'LineWidth', 1.5);
xlabel('T'); ylabel('z');

%% ==============================================================
%  CASE B: Tin ABOVE maximum temperature (hot inflow)
% ==============================================================
nexttile; hold on; grid on;
title('Extrapolated FK: Tin > Tmax(Tmix)');

for dt = dt_list

    T = T0;

    % ------------------------------------------------------------
    % (1) Extrapolated free convection operator
    % ------------------------------------------------------------

    ids = (z >= 10 & z <= 20);        % downward mixing region
    z_in = 10;                       % [m]
    z_mix = 20;                       % [m]

    Tmix = T(ids);
    zloc = z(ids);

    Tin = max(Tmix) + 10;            % <<< extrapolated HOT inlet

    I_num = abs(sum((Tin - Tmix)) * dz);

    Qdot_mix = 500000;                 % [W]
    A_hws    = 300;                  % [m^2]
    rho      = rho_w;
    cp       = cp_w * 1e3;           % [J/(kg*K)]

    dt_mix = dt;
    E = abs(Qdot_mix * dt_mix) / (rho * cp * A_hws);

    mix_height = abs(z_mix - z_in);
    dTdz = (Tmix(end) - Tmix(1)) / mix_height;

    dz_star = abs((Tin - Tmix(end)) / dTdz);

    z0 = 0;
    z1 = mix_height;
    z1s = z1 + dz_star;

    G = 0.5*(z1^2 - z0^2) - z1s*(z1 - z0);
    a = (E - I_num) / G;
    disp(a);
    z_phys = (zloc - z_mix);
    T(ids) = Tin + a * (z_phys - dz_star);

    if enableDiffusion
        T = diffuse_dirichlet_explicit(T, alpha, dz, dt, T_left, T_right);
    end

    plot(T, z, 'LineWidth', 1.5);
end

plot(T0, z, 'k--', 'LineWidth', 1.5);
xlabel('T'); ylabel('z');

%% ==============================================================
% Local diffusion function (unchanged)
% ==============================================================
function T = diffuse_dirichlet_explicit(T, alpha, dz, dt_total, T_left, T_right)

    if dt_total <= 0
        T(1)   = T_left;
        T(end) = T_right;
        return;
    end

    dt_stable = 0.45 * dz^2 / (2*alpha);
    nSteps = max(1, ceil(dt_total / dt_stable));
    dt = dt_total / nSteps;

    N = numel(T);
    for k = 1:nSteps
        Tnew = T;
        Tzz = (T(3:N) - 2*T(2:N-1) + T(1:N-2)) / dz^2;
        Tnew(2:N-1) = T(2:N-1) + alpha * dt * Tzz;
        Tnew(1)   = T_left;
        Tnew(end) = T_right;
        T = Tnew;
    end
end
