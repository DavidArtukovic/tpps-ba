%% test_energy_inflow.m
clear; clc;

%% ------------------------
% Geometry & discretization
%% ------------------------

use_diffusion = true;   % true: with diffusion, false: pure advection

H_w = 1.0;      % water height [m]
H_i = 1.0;      % insulation height [m]
dz  = 0.005;

z_w = (0:dz:H_w).';
z_i = (0:dz:H_i).';

Nw = length(z_w);
Ni = length(z_i);

%% ------------------------
% Material parameters
%% ------------------------
rho_cp_w = 1000*4200;   % water [J/(m3 K)]
rho_cp_i = 170*850;
thermal_eff_w = 1587.5;
thermal_eff_i = 107.5;

A = 255;

alpha_w = 0.6/(1000*4200);     % water thermal diffusivity [m^2/s]
alpha_i = 0.08/(170*850);      % insulation thermal diffusivity [m^2/s]  (always on)


if ~use_diffusion
    alpha_w = 0.0;
end

%% ------------------------
% Flow & temperatures
%% ------------------------
flow = 5e-5;        % inflow from top
T_init = 45.0;
T_in   = 80.0;

%% ------------------------
% Time discretization
%% ------------------------
dt = 1;
t_end = 4*3600;
t = 0:dt:t_end;


%% ------------------------
% Time Stamps for Plotting
%% ------------------------

Nt = numel(t);

Tw_hist_A = zeros(Nw, Nt);
Tw_hist_B = zeros(Nw, Nt);

Ti_hist_A = zeros(Ni, Nt);
Ti_hist_B = zeros(Ni, Nt);

%% ------------------------
% Initial conditions
%% ------------------------
T_w0 = T_init * ones(Nw,1);
T_i0 = T_init * ones(Ni,1);


% initial total energy in system
E0 = A * ( rho_cp_w * trapz(z_w, T_w0) + rho_cp_i * trapz(z_i, T_i0) );

% ground truth energy (E0 + inflow - outflow)
E_gt_old = zeros(size(t));
E_gt_old(1) = E0;

E_gt_new = zeros(size(t));
E_gt_new(1) = E0;

CFL = abs(flow) * dt / dz;
fprintf("CFL = %.3f\n", CFL);
assert(CFL < 0.8, "Advection likely unstable (CFL too high). Reduce dt or flow.");
%% ========================
% METHOD A: old (contact temp)
%% ========================
T_w = T_w0;
T_i = T_i0;
E_old = zeros(size(t));

for n = 1:length(t)
    % --- energy ---
    % ground truth energy balance
    Ein  = abs(flow) * A*rho_cp_w * T_in * dt;          % inflowing energy
    Eout = abs(flow) * A*rho_cp_w * T_w(end) * dt;      % outflow at bottom
    if n>1
        E_gt_old(n) = E_gt_old(n-1) + sum(Ein - Eout);
    end

    % --- contact temperature ---
    T_contact = (thermal_eff_w*T_in + thermal_eff_i*T_i(2)) / ...
                (thermal_eff_w + thermal_eff_i);

    % overwrite top node via contact
    T_i(1) = T_contact;
    T_w(1) = T_contact;

    % --- advection (upwind, no diffusion) ---
    dTw = zeros(Nw,1);

    for i = 2:Nw-1
        adv  = -flow * (T_w(i) - T_w(i-1)) / dz;
        diff = alpha_w * (T_w(i+1) - 2*T_w(i) + T_w(i-1)) / dz^2;
        dTw(i) = adv + diff;
    end
    % --- bottom water boundary (adiabatic for diffusion) ---
    % dT/dz = 0  -> ghost node: T_w(Nw+1) = T_w(Nw-1)
    adv_bottom  = -flow * (T_w(Nw) - T_w(Nw-1)) / dz;            % keep upwind advection
    diff_bottom = alpha_w * (2*(T_w(Nw-1) - T_w(Nw))) / dz^2;    % Neumann(0) diffusion

    dTw(Nw) = adv_bottom + diff_bottom;
    T_w =T_w + dt*dTw;
    Tw_hist_A(:, n) = T_w;

    % --- insulation diffusion (always on) ---
    dTi = zeros(Ni,1);
    dTi(1) = 0;
    % interior nodes
    for j = 2:Ni-1
        dTi(j) = alpha_i * (T_i(j+1) - 2*T_i(j) + T_i(j-1)) / dz^2;
    end
    
    % adiabatic at insulation top: dT/dz = 0  -> ghost node T(Ni+1)=T(Ni-1)
    dTi(Ni) = alpha_i * (T_i(Ni-1) - 2*T_i(Ni) + T_i(Ni-1)) / dz^2;
    T_i = T_i + dt*dTi;
    
    Ti_hist_A(:, n) = T_i;
    E_old(n) = A * ( rho_cp_w * trapz(z_w, T_w) + rho_cp_i * trapz(z_i, T_i) );

end

%% ========================
% METHOD B: new (explicit inflow)
%% ========================

T_w = T_w0;
T_i = T_i0;
E_new = zeros(size(t));

for n = 1:length(t)
    % --- energy ---
    % ground truth energy balance
    Ein  = abs(flow) * A*rho_cp_w * T_in * dt;          % inflowing energy
    Eout = abs(flow) * A*rho_cp_w * T_w(end) * dt;      % outflow at bottom
    if n>1
        E_gt_new(n) = E_gt_new(n-1) + sum(Ein - Eout);
    end
    
    % effective heat transfer coefficient [W/(m^2 K)]
    h_int = 1 / (1/thermal_eff_w + 1/thermal_eff_i);
    
    % heat flux water -> insulation [W/m^2]
    q_int = h_int * (T_w(1) - T_i(1));

    % --- advection ---
    dTw = zeros(Nw,1);
    dTw(1) = -2*(flow/dz) * (T_w(1) - T_in) ...
       + alpha_w * (2*(T_w(2) - T_w(1))) / dz^2 ...  % Neumann(0) diffusion
       - 2*q_int / (rho_cp_w * dz);

    for i = 2:Nw-1
        adv  = -flow * (T_w(i) - T_w(i-1)) / dz;
        diff = alpha_w * (T_w(i+1) - 2*T_w(i) + T_w(i-1)) / dz^2;
        dTw(i) = adv + diff;
    end
    % --- bottom water boundary (adiabatic for diffusion) ---
    % dT/dz = 0  -> ghost node: T_w(Nw+1) = T_w(Nw-1)
    adv_bottom  = -flow * (T_w(Nw) - T_w(Nw-1)) / dz;            % keep upwind advection
    diff_bottom = alpha_w * (2*(T_w(Nw-1) - T_w(Nw))) / dz^2;    % Neumann(0) diffusion

    dTw(Nw) = adv_bottom + diff_bottom;
    
    T_w =T_w+ dt*dTw;
    Tw_hist_B(:, n) = T_w;

    % --- insulation diffusion (always on) ---
    dTi = zeros(Ni,1);

    dTi(1) = ...
    alpha_i * 2*(T_i(2) - T_i(1)) / dz^2 ...   % axial diffusion into insulation
    + 2*q_int / (rho_cp_i * dz);               % interface heat flux

    % interior nodes
    for j = 2:Ni-1
        dTi(j) = alpha_i * (T_i(j+1) - 2*T_i(j) + T_i(j-1)) / dz^2;
    end
    
    % adiabatic at insulation top: dT/dz = 0  -> ghost node T(Ni+1)=T(Ni-1)
    dTi(Ni) = alpha_i * (T_i(Ni-1) - 2*T_i(Ni) + T_i(Ni-1)) / dz^2;
    T_i = T_i + dt*dTi;
    
    Ti_hist_B(:, n) = T_i;
    E_new(n) = A * ( rho_cp_w * trapz(z_w, T_w) + rho_cp_i * trapz(z_i, T_i) );

end

%% ------------------------
% Plot energy comparison
%% ------------------------
figure;
plot(t/3600, E_old, 'r-', 'LineWidth', 1.5); 
hold on;
plot(t/3600, E_new, 'b-',  'LineWidth', 1.5);
hold on;
plot(t/3600, E_gt_old,  'k--',  'LineWidth', 2);
hold on;

xlabel('Time [h]');
ylabel('Total energy [arb.]');
legend('Old method','New method','Ground truth','Location','best');
grid on;


figure; hold on;
plot(t/3600, (E_old - E_gt_old)./E_gt_old, 'r');
plot(t/3600, (E_new - E_gt_new)./E_gt_new, 'b');
yline(0,'k:');
xlabel('Time [h]');
ylabel('relative Energy error');
legend('Method A','Method B','Location','best');
grid on;

%%

T_wA_end = Tw_hist_A(:, end);
T_wB_end = Tw_hist_B(:, end);

T_iA_end = Ti_hist_A(:,end);   % oder Ti_hist_A(:,end), falls du das auch speicherst
T_iB_end = Ti_hist_B(:,end);

figure; hold on;

z_i_plot = -flipud(z_i);   % insulation: negative
z_w_plot =  z_w;           % water: positive

z_all = [z_i_plot; z_w_plot];

T_A = [flipud(T_iA_end); T_wA_end];
T_B = [flipud(T_iB_end); T_wB_end];

plot(T_A, z_all, '-',  'LineWidth', 1.8);
plot(T_B, z_all, '--', 'LineWidth', 1.8);

% interface
yline(0, 'k:', 'LineWidth', 1);

xlabel('Temperature [°C]');
ylabel('Height [m]');
title('Temperature profile at final time step');
legend('Method A','Method B','Location','best');
grid on;

