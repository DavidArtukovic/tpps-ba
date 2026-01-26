%% test_total_energy_balance.m
% ---------------------------------------------------------------
% PURPOSE:
%   Test global energy consistency (1st law) of the TPPS model.
%   Compares:
%     - Change of total internal system energy
%     - Cumulative supplied and removed flow energy
%
%   The balance defines an upper energetic bound:
%       E_in - E_out >= Delta E_sys
%
%   Valid for comparison of:
%     - MATLAB (no FK)
%     - MATLAB (FK v15)
%     - MATLAB (FK v16)
%
% AUTHOR:
%   BA TPPS – Energy consistency check
% ---------------------------------------------------------------

clc
clear
close all

%%%------------------------------------------%%%
% 01. Load paths and data
%%%------------------------------------------%%%

run(fullfile('..','..','..', 'configs', 'paths_local.m'));

DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_OWN  = fullfile(DATA_BASE, 'scenario1_freeConv');
DATA_SCEN1      = fullfile(DATA_BASE, 'scenario1');

% load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));

load(fullfile(DATA_SCEN1, 'Init_d18_h18_time8.mat'));   % Geometry, material values and initial values

data_no_fk = load(fullfile(DATA_SCEN1_BASE, ...
    'd18_h18_Res_Matlab_d18_18.mat'));

% data_fk_v15 = load(fullfile(DATA_SCEN1_OWN, ...
%     '260122_d18_h18_Res_Matlab_FK_v15.mat'));

data_fk_v16 = load(fullfile(DATA_SCEN1_OWN, ...
    '260124_d18_h18_Res_Matlab_FK_v16.mat'));
%%
%%%------------------------------------------%%%
% 02. Physical parameters
%%%------------------------------------------%%%

% Geometric radii and areas
d_ST   = d(1);           % Diameter in m
r_ST   = d_ST / 2;       % Radius in m

A(1) = d(1)^2 * pi/4;    % Area of the cylindrical storage

dt   = 900;              % timestep [s]
rho  = SW(1,2);          % water density [kg/m^3]
cp   = SW(1,3);          % heat capacity [J/(kg K)]
Acs  = A(1);             % storage cross section [m^2]

Tin  = 80;               % supply temperature [°C]
Tout = 45;               % return temperature [°C]

%%%------------------------------------------%%%
% 03. Helper function handle
%%%------------------------------------------%%%

compute_balance = @(Res_900) compute_energy_balance_run( ...
    Res_900, rho, cp, Acs, Tin, Tout, dt);

%%%------------------------------------------%%%
% 04. Compute balances
%%%------------------------------------------%%%

[E_sys_noFK, E_flow_noFK] = compute_balance(data_no_fk.Res_900_d18_18);
% [E_sys_v15,  E_flow_v15 ] = compute_balance(data_fk_v15.Res_900_d18_18_FK);
[E_sys_v16,  E_flow_v16 ] = compute_balance(data_fk_v16.Res_900_d18_18_FK);
%%
%%%------------------------------------------%%%
% 05. Plot results
%%%------------------------------------------%%%

t_h = (0:length(E_sys_noFK)-1) * dt / 3600;

figure('Color','w'); hold on; grid on;
Nplot = 960;
plot(t_h(1:Nplot), E_sys_noFK(1:Nplot), 'k-',  'LineWidth',1.4)
% plot(t_h(1:Nplot), E_sys_v15(1:Nplot),  'g-',  'LineWidth',1.4)
plot(t_h(1:Nplot), E_sys_v16(1:Nplot),  'b-',  'LineWidth',1.4)

plot(t_h(1:Nplot), E_flow_noFK(1:Nplot),'k--', 'LineWidth',1.2)
% plot(t_h(1:Nplot), E_flow_v15(1:Nplot), 'g--', 'LineWidth',1.2)
plot(t_h(1:Nplot), E_flow_v16(1:Nplot), 'b--', 'LineWidth',1.2)

xlabel('Time [h]')
ylabel('Energy [J]')
title('Total energy balance – system energy vs. flow bound')

legend( ...
    'ΔE_{sys} no FK', ...
    'ΔE_{sys} FK v16', ...
    'E_{in}-E_{out} no FK', ...
    'E_{in}-E_{out} FK v16', ...
    'Location','best')

%%%------------------------------------------%%%
% 06. Residual diagnostics
%%%------------------------------------------%%%

res_noFK = E_flow_noFK - E_sys_noFK;
res_v15  = E_flow_v15  - E_sys_v15;
res_v16  = E_flow_v16  - E_sys_v16;

fprintf('\n===== Energy balance residuals =====\n');
fprintf('no FK : min = %.3e J | max = %.3e J\n', min(res_noFK), max(res_noFK));
% fprintf('FK v15: min = %.3e J | max = %.3e J\n', min(res_v15 ), max(res_v15 ));
fprintf('FK v16: min = %.3e J | max = %.3e J\n', min(res_v16 ), max(res_v16 ));

%% ============================================================= %%
%                        LOCAL FUNCTION
%% ============================================================= %%

function [E_sys, E_flow] = compute_energy_balance_run( ...
    Res_900, rho, cp, A, Tin, Tout, dt)
% COMPUTE_ENERGY_BALANCE_RUN
%   Computes:
%     - ΔE_sys : total internal energy change
%     - E_flow : cumulative flow energy (upper bound)
%
%   Assumes:
%     - Res_900(1:6,:) store absolute energies [J]
%     - Res_900(9,:)   charging velocity  [m/s]
%     - Res_900(10,:)  discharging velocity [m/s]

    % --- System energy change ---
    E_abs = sum(Res_900(1:6,:), 1);
    E_sys = E_abs - E_abs(1);

    % --- Signed flow velocity ---
    v = Res_900(9,2:end) + Res_900(10,2:end);

    % --- Cumulative flow energy ---
    dE = -rho * cp .* (v * A) .* (Tin - Tout) * dt;
    E_flow = [0, cumsum(dE)];
end
