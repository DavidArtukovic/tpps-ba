%% test_energy_inflow_2vol_bypass_soil.m
% -------------------------------------------------------------------------
%  SUMMARY:
%   Energy consistency test for a 1D vertical TPPS-like system consisting
%   of two hydraulically connected water volumes (bypass), an insulation
%   layer above, and a soil layer below.
%
%  PURPOSE:
%   The script compares two numerical modeling strategies for handling
%   thermal interfaces and advective energy transport:
%
%   Method A ("contact temperature approach"):
%     - Water–solid interfaces are enforced via a contact temperature.
%     - No explicit interfacial heat flux terms are added.
%     - Equivalent to the legacy TPPS implementation using interface nodes.
%
%   Method B ("explicit flux approach"):
%     - Water–solid interfaces are modeled via explicit heat fluxes.
%     - Inflow and outflow are treated as energetic boundary terms.
%     - No artificial contact temperature is imposed.
%
%  PHYSICAL SETUP:
%   - Upper water volume: inflow at the top, advection downward to bypass.
%   - Lower water volume: inflow from bypass, advection downward to outlet.
%   - Bypass is assumed lossless and purely advective.
%   - Insulation and soil exchange heat only via conduction.
%   - All internal interfaces (water–water, water–solid) are adiabatic
%     unless explicitly modeled via heat fluxes.
%
%  OBJECTIVE:
%   - Verify global energy conservation:
%       E(t) = E(0) + Ein(t) - Eout(t)
%   - Identify systematic differences between interface modeling strategies.
%   - Detect numerical artifacts caused by interface discretization.
%
%  REMARKS:
%   - Spatial discretization: uniform FD grid.
%   - Time integration: explicit Euler.
%   - Advection: first-order upwind.
%   - Diffusion: second-order central differences.
%
% -------------------------------------------------------------------------

clear; clc;

%%% ------------------------
% 1. Geometry & discretization
%%% ------------------------
use_diffusion = true;   % true: with diffusion, false: pure advection

H_w1 = 1.0;      % upper water height [m]
H_w2 = 1.0;      % lower water height [m]
H_i  = 1.0;      % insulation height [m] (above upper water)
H_s  = 1.0;      % soil height [m] (below lower water)

dz  = 0.005;

z_w1 = (0:dz:H_w1).';
z_w2 = (0:dz:H_w2).';
z_i  = (0:dz:H_i).';
z_s  = (0:dz:H_s).';

Nw1 = length(z_w1);
Nw2 = length(z_w2);
Ni  = length(z_i);
Ns  = length(z_s);

%%% ------------------------
% 2. Bypass location
%%% ------------------------
z_bypass = 0.2;                              % [m] in each volume (local coordinate)
idx_b1   = round(z_bypass/dz) + 1;           % in upper water
idx_b2   = round(0.5/dz) + 1;           % in lower water

assert(idx_b1 >= 2 && idx_b1 <= Nw1-1, "Bypass index in water1 must be interior.");
assert(idx_b2 >= 2 && idx_b2 <= Nw2-1, "Bypass index in water2 must be interior.");

%%% ------------------------
% 3. Material parameters
%%% ------------------------
% Water
rho_cp_w      = 1000*4200;                   % [J/(m^3 K)]
alpha_w       = 0.6/(1000*4200);             % [m^2/s]
thermal_eff_w = 1587.5;                      % [J/(m^2 K s^0.5)] (given)

% Insulation
rho_cp_i      = 170*850;
alpha_i       = 0.08/(170*850);              % insulation diffusion always on
thermal_eff_i = 107.5;

% Soil
lambda_s      = 2.1;                         % [W/mK]
rho_s         = 2400;                        % [kg/m^3]
cp_s          = 1000;                        % [J/kgK]
rho_cp_s      = rho_s*cp_s;                  % [J/(m^3 K)]
alpha_s       = lambda_s/(rho_s*cp_s);       % [m^2/s]
thermal_eff_s = sqrt(lambda_s*rho_s*cp_s);   % [J/(m^2 K s^0.5)]

A = 255;                                     % area factor (as in your script)

if ~use_diffusion
    alpha_w = 0.0;
end

%%% ------------------------
% 4. Flow & temperatures
%%% ------------------------
flow = 5e-5;        % positive: increasing index direction (downwards)
T_init = 45.0;
T_in   = 80.0;

%%% ------------------------
% 5.  Time discretization
%%% ------------------------
dt = 0.1;
t_end = 12*3600;
t = 0:dt:t_end;
Nt = numel(t);

%%% ------------------------
% 6. Temperature matrices
%%% ------------------------
Tw1_hist_A = zeros(Nw1, Nt);
Tw2_hist_A = zeros(Nw2, Nt);
Ti_hist_A  = zeros(Ni,  Nt);
Ts_hist_A  = zeros(Ns,  Nt);

Tw1_hist_B = zeros(Nw1, Nt);
Tw2_hist_B = zeros(Nw2, Nt);
Ti_hist_B  = zeros(Ni,  Nt);
Ts_hist_B  = zeros(Ns,  Nt);

%%% ------------------------
% 7. Initial conditions
%%% ------------------------
T_w1_0 = T_init * ones(Nw1,1);
T_w2_0 = T_init * ones(Nw2,1);
T_i0   = T_init * ones(Ni,1);
T_s0   = T_init * ones(Ns,1);

% Initial total energy in system
E0 = A * ( ...
    rho_cp_i * trapz(z_i,  T_i0)  + ...
    rho_cp_w * trapz(z_w1, T_w1_0)+ ...
    rho_cp_w * trapz(z_w2, T_w2_0)+ ...
    rho_cp_s * trapz(z_s,  T_s0) );

% Ground truth energy: E(t) = E0 + inflow - outflow (bypass is internal)
E_gt_old = zeros(size(t)); E_gt_old(1) = E0;
E_gt_new = zeros(size(t)); E_gt_new(1) = E0;

CFL = abs(flow) * dt / dz;
fprintf("CFL = %.3f\n", CFL);
assert(CFL < 0.8, "Advection likely unstable (CFL too high). Reduce dt or flow.");

%%% ========================
% METHOD A: old (contact temp)
%%% ========================
T_w1 = T_w1_0;
T_w2 = T_w2_0;
T_i  = T_i0;
T_s  = T_s0;

E_old = zeros(size(t));

for n = 1:Nt
    %%% ------------------------
    % ground truth energy balance (system boundary only)
    %%% ------------------------
    Ein  = abs(flow) * A * rho_cp_w * T_in       * dt;      % inflow at top
    Eout = abs(flow) * A * rho_cp_w * T_w2(end)  * dt;      % outflow at bottom of lower water
    if n > 1
        E_gt_old(n) = E_gt_old(n-1) + (Ein - Eout);
    end

    % --- top contact temperature (water1 <-> insulation, with T_in reference) ---
    T_contact_top = (thermal_eff_w*T_in + thermal_eff_i*T_i(2)) / ...
                    (thermal_eff_w + thermal_eff_i);
    T_i(1)  = T_contact_top;
    T_w1(1) = T_contact_top;

    % --- bottom contact temperature (water2 <-> soil, symmetric contact) ---
    % uses adjacent interior nodes as "bulk" references
    T_contact_bot = (thermal_eff_w*T_w2(end-1) + thermal_eff_s*T_s(2)) / ...
                    (thermal_eff_w + thermal_eff_s);
    T_w2(end) = T_contact_bot;
    T_s(1)    = T_contact_bot;

    % --- bypass coupling (lossless): inlet temperature into lower = outlet temp from upper ---
    T_bypass = T_w1(idx_b1);

    %%% ------------------------
    % Upper water volume: advection only from top down to bypass, diffusion everywhere
    %%% ------------------------

    dTw1 = zeros(Nw1,1);

    % cell 2..idx_b1: flowing segment
    for i = 2:idx_b1
        adv  = -flow * (T_w1(i) - T_w1(i-1)) / dz;
        diff = alpha_w * (T_w1(i+1) - 2*T_w1(i) + T_w1(i-1)) / dz^2;
        dTw1(i) = adv + diff;
    end

    % below bypass: no flow, only diffusion
    for i = idx_b1+1:Nw1-1
        diff = alpha_w * (T_w1(i+1) - 2*T_w1(i) + T_w1(i-1)) / dz^2;
        dTw1(i) = diff;
    end

    % bottom of upper water: adiabatic for diffusion
    diff_bottom1 = alpha_w * (2*(T_w1(Nw1-1) - T_w1(Nw1))) / dz^2; % Neumann(0)
    dTw1(Nw1) = diff_bottom1;   % no advection here (no flow)

    T_w1 = T_w1 + dt*dTw1;
    
    %%% ------------------------
    % Lower water volume: advection only from bypass downwards, diffusion everywhere
    %%% ------------------------

    dTw2 = zeros(Nw2,1);

    % top of lower water: adiabatic diffusion (no flow above bypass)
    dTw2(1) = alpha_w * (2*(T_w2(2) - T_w2(1))) / dz^2;

    % above bypass: no flow, only diffusion
    for i = 2:idx_b2-1
        diff = alpha_w * (T_w2(i+1) - 2*T_w2(i) + T_w2(i-1)) / dz^2;
        dTw2(i) = diff;
    end

    % bypass inlet cell: enforce advective inlet from T_bypass (boundary-like)
    dTw2(idx_b2) = -(flow/dz) * (T_w2(idx_b2) - T_bypass) + ...
                   alpha_w * (T_w2(idx_b2+1) - 2*T_w2(idx_b2) + T_w2(idx_b2-1)) / dz^2;

    % flowing interior below bypass
    for i = idx_b2+1:Nw2-1
        adv  = -flow * (T_w2(i) - T_w2(i-1)) / dz;
        diff = alpha_w * (T_w2(i+1) - 2*T_w2(i) + T_w2(i-1)) / dz^2;
        dTw2(i) = adv + diff;
    end

    % bottom cell already overwritten by contact; keep derivative zero
    dTw2(Nw2) = 0;

    T_w2 = T_w2 + dt*dTw2;

    %%% ------------------------
    % Insulation diffusion (always on): adiabatic at top
    %%% ------------------------

    dTi = zeros(Ni,1);
    dTi(1) = 0; % contact-enforced
    for j = 2:Ni-1
        dTi(j) = alpha_i * (T_i(j+1) - 2*T_i(j) + T_i(j-1)) / dz^2;
    end
    dTi(Ni) = alpha_i * (T_i(Ni-1) - 2*T_i(Ni) + T_i(Ni-1)) / dz^2; % adiabatic at top
    T_i = T_i + dt*dTi;

    %%% ------------------------
    % Soil diffusion (always on): adiabatic at bottom
    %%% ------------------------

    dTs = zeros(Ns,1);
    dTs(1) = 0; % contact-enforced
    for j = 2:Ns-1
        dTs(j) = alpha_s * (T_s(j+1) - 2*T_s(j) + T_s(j-1)) / dz^2;
    end
    dTs(Ns) = alpha_s * (T_s(Ns-1) - 2*T_s(Ns) + T_s(Ns-1)) / dz^2; % adiabatic at bottom
    T_s = T_s + dt*dTs;

    % store
    Tw1_hist_A(:,n) = T_w1;
    Tw2_hist_A(:,n) = T_w2;
    Ti_hist_A(:,n)  = T_i;
    Ts_hist_A(:,n)  = T_s;

    % total energy
    E_old(n) = A * ( ...
        rho_cp_i * trapz(z_i,  T_i)  + ...
        rho_cp_w * trapz(z_w1, T_w1) + ...
        rho_cp_w * trapz(z_w2, T_w2) + ...
        rho_cp_s * trapz(z_s,  T_s) );
end

%%% ========================
% METHOD B: new (explicit inflow + explicit interface fluxes)
%%% ========================
T_w1 = T_w1_0;
T_w2 = T_w2_0;
T_i  = T_i0;
T_s  = T_s0;

E_new = zeros(size(t));

for n = 1:Nt

    %%% ------------------------
    % ground truth energy balance (system boundary only)
    %%% ------------------------
    Ein  = abs(flow) * A * rho_cp_w * T_in      * dt;
    Eout = abs(flow) * A * rho_cp_w * T_w2(end) * dt;
    if n > 1
        E_gt_new(n) = E_gt_new(n-1) + (Ein - Eout);
    end

    % --- bypass coupling (lossless): inlet temperature into lower = outlet temp from upper ---
    T_bypass = T_w1(idx_b1);

    %%% ------------------------
    % Top interface (water1 <-> insulation)
    %%% ------------------------

    h_int_wi = 1 / (1/thermal_eff_w + 1/thermal_eff_i);
    q_wi     = h_int_wi * (T_w1(1) - T_i(1));     % +: water -> insulation

    %%% ------------------------
    % Bottom interface (water2 <-> soil)
    %%% ------------------------

    h_int_ws = 1 / (1/thermal_eff_w + 1/thermal_eff_s);
    q_ws     = h_int_ws * (T_w2(end) - T_s(1));   % +: water -> soil

    %%% ------------------------
    % Upper water volume (flow only to bypass)
    %%% ------------------------

    dTw1 = zeros(Nw1,1);

    % top water cell: inflow + diffusion (Neumann) + interface heat flux
    dTw1(1) = -2*(flow/dz) * (T_w1(1) - T_in) ...
              + alpha_w * (2*(T_w1(2) - T_w1(1))) / dz^2 ...
              - 2*q_wi / (rho_cp_w * dz);

    % flowing segment up to bypass
    for i = 2:idx_b1
        adv  = -flow * (T_w1(i) - T_w1(i-1)) / dz;
        diff = alpha_w * (T_w1(i+1) - 2*T_w1(i) + T_w1(i-1)) / dz^2;
        dTw1(i) = adv + diff;
    end

    % below bypass: no flow, only diffusion
    for i = idx_b1+1:Nw1-1
        diff = alpha_w * (T_w1(i+1) - 2*T_w1(i) + T_w1(i-1)) / dz^2;
        dTw1(i) = diff;
    end

    % bottom of upper water: adiabatic diffusion
    dTw1(Nw1) = alpha_w * (2*(T_w1(Nw1-1) - T_w1(Nw1))) / dz^2;

    T_w1 = T_w1 + dt*dTw1;

    %%% ------------------------
    % Lower water volume (flow only from bypass downwards)
    %%% ------------------------

    dTw2 = zeros(Nw2,1);

    % top of lower water: adiabatic diffusion (no flow above bypass)
    dTw2(1) = alpha_w * (2*(T_w2(2) - T_w2(1))) / dz^2;

    % above bypass: no flow, only diffusion
    for i = 2:idx_b2-1
        diff = alpha_w * (T_w2(i+1) - 2*T_w2(i) + T_w2(i-1)) / dz^2;
        dTw2(i) = diff;
    end

    % bypass inlet cell: inflow from T_bypass + diffusion
    dTw2(idx_b2) = -(flow/dz) * (T_w2(idx_b2) - T_bypass) ...
                   + alpha_w * (T_w2(idx_b2+1) - 2*T_w2(idx_b2) + T_w2(idx_b2-1)) / dz^2;

    % flowing interior
    for i = idx_b2+1:Nw2-1
        adv  = -flow * (T_w2(i) - T_w2(i-1)) / dz;
        diff = alpha_w * (T_w2(i+1) - 2*T_w2(i) + T_w2(i-1)) / dz^2;
        dTw2(i) = adv + diff;
    end

    % bottom water cell: outflow advection + diffusion to satisfy Neumann? + interface heat flux to soil
    adv_bottom2  = -2*flow * (T_w2(Nw2) - T_w2(Nw2-1)) / dz;
    diff_bottom2 = alpha_w * (2*(T_w2(Nw2-1) - T_w2(Nw2))) / dz^2;
    dTw2(Nw2)    = adv_bottom2 + diff_bottom2 - 2*q_ws/(rho_cp_w * dz);

    T_w2 = T_w2 + dt*dTw2;

    %%% ------------------------
    % Insulation diffusion (always on) + interface flux from water
    %%% ------------------------

    dTi = zeros(Ni,1);
    dTi(1) = alpha_i * 2*(T_i(2) - T_i(1)) / dz^2 + 2*q_wi / (rho_cp_i * dz);
    for j = 2:Ni-1
        dTi(j) = alpha_i * (T_i(j+1) - 2*T_i(j) + T_i(j-1)) / dz^2;
    end
    dTi(Ni) = alpha_i * (T_i(Ni-1) - 2*T_i(Ni) + T_i(Ni-1)) / dz^2;
    T_i = T_i + dt*dTi;

    %%% ------------------------
    % Soil diffusion (always on) + interface flux from water
    %%% ------------------------

    dTs = zeros(Ns,1);
    dTs(1) = alpha_s * 2*(T_s(2) - T_s(1)) / dz^2 + 2*q_ws / (rho_cp_s * dz);
    for j = 2:Ns-1
        dTs(j) = alpha_s * (T_s(j+1) - 2*T_s(j) + T_s(j-1)) / dz^2;
    end
    dTs(Ns) = alpha_s * (T_s(Ns-1) - 2*T_s(Ns) + T_s(Ns-1)) / dz^2; % adiabatic at bottom
    T_s = T_s + dt*dTs;

    % store
    Tw1_hist_B(:,n) = T_w1;
    Tw2_hist_B(:,n) = T_w2;
    Ti_hist_B(:,n)  = T_i;
    Ts_hist_B(:,n)  = T_s;

    % total energy
    E_new(n) = A * ( ...
        rho_cp_i * trapz(z_i,  T_i)  + ...
        rho_cp_w * trapz(z_w1, T_w1) + ...
        rho_cp_w * trapz(z_w2, T_w2) + ...
        rho_cp_s * trapz(z_s,  T_s) );
end

%%% ------------------------
% 8. Plot energy comparison
%%% ------------------------
figure;
plot(t/3600, E_old, 'r-', 'LineWidth', 1.5); hold on;
plot(t/3600, E_new, 'b-', 'LineWidth', 1.5); hold on;
plot(t/3600, E_gt_new, 'k--', 'LineWidth', 2);
xlabel('Time [h]'); ylabel('Total energy [arb.]');
legend('Method A','Method B','Ground truth','Location','best');
grid on;

figure; hold on;
plot(t/3600, (E_old - E_gt_old)./E_gt_old, 'r');
plot(t/3600, (E_new - E_gt_new)./E_gt_new, 'b');
yline(0,'k:');
xlabel('Time [h]'); ylabel('relative Energy error');
legend('Method A','Method B','Location','best');
grid on;

%%% ------------------------
% 90. Plot final temperature profile (stacked)
%%% ------------------------
Tw1A_end = Tw1_hist_A(:, end);
Tw2A_end = Tw2_hist_A(:, end);
TiA_end  = Ti_hist_A(:,  end);
TsA_end  = Ts_hist_A(:,  end);

Tw1B_end = Tw1_hist_B(:, end);
Tw2B_end = Tw2_hist_B(:, end);
TiB_end  = Ti_hist_B(:,  end);
TsB_end  = Ts_hist_B(:,  end);

% coordinate for plotting:
% insulation above interface: negative
z_i_plot  = -flipud(z_i);

% upper water: 0..H_w1
z_w1_plot = z_w1;

% lower water: H_w1..H_w1+H_w2
z_w2_plot = H_w1 + z_w2;

% soil: below lower water
z_s_plot  = H_w1 + H_w2 + z_s;

z_all = [z_i_plot; z_w1_plot; z_w2_plot; z_s_plot];

T_A = [flipud(TiA_end); Tw1A_end; Tw2A_end; TsA_end];
T_B = [flipud(TiB_end); Tw1B_end; Tw2B_end; TsB_end];

figure; hold on;
plot(T_A, z_all, '-',  'LineWidth', 1.8);
plot(T_B, z_all, '--', 'LineWidth', 1.8);

yline(0, 'k:', 'LineWidth', 2);                 % insulation-water interface
yline(1, 'k:', 'LineWidth', 2);                 % water-water interface
yline(2, 'k:', 'LineWidth', 2);                 % water-soil interface

yline(H_w1, 'k:', 'LineWidth', 1);              % between water1 & water2 (adiabatic separation)
yline(H_w1+H_w2, 'k:', 'LineWidth', 1);         % water2-soil interface

xline(T_in, 'k:', 'LineWidth', 1);              % inlet temperature reference (optional)

xlabel('Temperature [°C]');
ylabel('Height [m]');
title('Final stacked temperature profile');
legend('Method A','Method B','Location','best');
grid on;

% --- textual annotations for adiabatic interfaces ---
text(min(T_A)+5, H_w1, ...
     '\textit{adiabatic interface (water--water)}', ...
     'Interpreter','latex', ...
     'VerticalAlignment','bottom');

text(min(T_A)+5, -H_i, ...
     '\textit{adiabatic boundary (top insulation)}', ...
     'Interpreter','latex', ...
     'VerticalAlignment','top');

text(min(T_A)+5, H_w1+H_w2+H_s, ...
     '\textit{adiabatic boundary (soil)}', ...
     'Interpreter','latex', ...
     'VerticalAlignment','bottom');