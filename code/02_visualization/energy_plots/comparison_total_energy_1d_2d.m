%% comparison_total_energy_1d_2d.m
% ---------------------------------------------------------------
% PURPOSE:
%   Compare the global energy balance of the legacy 1D TPPS model
%   and the newer 2D-piston TPPS model.
%
% DESCRIPTION:
%   - Loads one legacy 1D result file and one newer ResOut-based 2D result file.
%   - Maps both file formats to one common internal data structure.
%   - Computes:
%       (1) change of total internal system energy
%       (2) cumulative net flow enthalpy using the actual outlet temperature
%       (3) cumulative simple reference enthalpy using constant 80/45 °C
%   - Plots both model variants against their corresponding flow-based
%     ground-truth curves.
%   - Prints residual diagnostics for the exact and simplified balances.
%
% ENERGY BALANCE:
%   Charging  (v < 0):
%       hot water enters at the top with Tin_hot = 80 °C
%       water leaves at the bottom with T_bottom
%
%   Discharging (v > 0):
%       cold water enters at the bottom with Tin_cold = 45 °C
%       water leaves at the top with T_top
%
%   The cumulative flow enthalpy is compared against:
%       Delta E_sys = E_sys(t) - E_sys(0)
%
% NOTES:
%   - The script assumes that stored water temperatures are oriented such
%     that row 1 = top node and row end = bottom node.
%   - If the legacy file uses the opposite orientation, adapt the loader
%     function 'unpack_legacy_run' accordingly.
%
% ---------------------------------------------------------------

clc
clear
close all

%%%------------------------------------------%%%
% 01. Load paths and data
%%%------------------------------------------%%%

run(fullfile('..','..','..', 'configs', 'paths_local.m'));

DATA_INIT = fullfile(DATA_BASE, 'init');
DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_FK   = fullfile(DATA_BASE, 'scenario1_freeConv');
RESULTS_VIS_BASE  = fullfile(RESULTS_BASE, "02_visualizations");
RESULTS_ENERGY = fullfile(RESULTS_VIS_BASE, "01_energy");

if ~exist(RESULTS_ENERGY, 'dir')
    mkdir(RESULTS_ENERGY);
end

load(fullfile(DATA_INIT, '20260330_d18_hp18.0_gap0.5_2D_chunk009.mat'));   % Geometry, material values and initial values

% --- Legacy 1D result file
data_1d = load(fullfile(DATA_SCEN1_BASE, ...
    'd18_h18_Res_Matlab_d18_18.mat'));

% --- Newer 2D result file (ResOut structure)
data_2d = load(fullfile(DATA_SCEN1_FK, ...
    '260331_d18_h18_Res_Matlab_FK_2d_v12.mat'));

export_png_dpi = 600;
%%
%%%------------------------------------------%%%
% 02. Physical parameters
%%%------------------------------------------%%%

% Cross section of storage [m^2]
Acs = InitOut.geom.d(1)^2 * pi/4;

% Water properties
rho = InitOut.param.SW(1,2);              % density [kg/m^3]
cp  = InitOut.param.SW(1,3);              % heat capacity [J/(kg K)]

% Time step
dt = 900;                   % [s]

% Inlet temperatures [°C]
Tin_hot  = 80;
Tin_cold = 45;

%%
%%%------------------------------------------%%%
% 03. Map data files to unified structure
%%%------------------------------------------%%%

run_1d = unpack_legacy_run(data_1d, 'MATLAB 1D');
run_2d = unpack_resout_run(data_2d, 'MATLAB 2D piston');

%%
%%%------------------------------------------%%%
% 04. Compute energy balances
%%%------------------------------------------%%%

[E_sys_1d, E_flow_true_1d, E_flow_simple_1d] = ...
    compute_energy_balance_true( ...
        run_1d.Res_900, run_1d.T_W_900, ...
        rho, cp, Acs, Tin_hot, Tin_cold, dt);

[E_sys_2d, E_flow_true_2d, E_flow_simple_2d] = ...
    compute_energy_balance_true( ...
        run_2d.Res_900, run_2d.T_W_900, ...
        rho, cp, Acs, Tin_hot, Tin_cold, dt);

%%
%%%------------------------------------------%%%
% 05. Time axis
%%%------------------------------------------%%%

Nt_1d = length(E_sys_1d);
Nt_2d = length(E_sys_2d);

t_h_1d = (0:Nt_1d-1) * dt / 3600;
t_h_2d = (0:Nt_2d-1) * dt / 3600;

% convert to days
t_d_1d = t_h_1d / 24;
t_d_2d = t_h_2d / 24;

Nplot_1d = Nt_1d;
Nplot_2d = Nt_2d;

%%
%%%------------------------------------------%%%
% 06. Plot: total system energy vs. flow enthalpy
%%%------------------------------------------%%%

% Academic / muted colors
c_base = [0.15 0.15 0.15];   % dark gray
c_cmsl = [0.00 0.45 0.70];   % blue
c_ext  = [0.80 0.40 0.00];   % orangen

width_cm  = 14;
height_cm = 9;

fig = figure( ...
    'Color', 'w', ...
    'Name', 'MATLAB energy comparison', ...
    'Units', 'centimeters', ...
    'Position', [2 2 width_cm height_cm]);

set(fig, 'Units', 'centimeters');
set(fig, 'Position', [2 2 width_cm height_cm]);

set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [width_cm height_cm]);
set(fig, 'PaperPosition', [0 0 width_cm height_cm]);
set(fig, 'Renderer', 'painters');

ax = gca;
set(ax,'FontName','Times New Roman');
set(ax,'FontSize',9);
set(ax,'LineWidth',0.5);

set(gca,'GridAlpha',0.10);             % lighter grid
set(gca,'MinorGridAlpha',0.05);

hold on
grid on
box on

% --- 1D model
plot(t_d_1d(1:Nplot_1d), E_sys_1d(1:Nplot_1d), ...
    '-','Color', c_base, 'LineWidth', 0.7);

plot(t_d_1d(1:Nplot_1d), E_flow_true_1d(1:Nplot_1d), ...
    '--', 'Color', c_base, 'LineWidth', 0.7);

% --- 2D model
plot(t_d_2d(1:Nplot_2d), E_sys_2d(1:Nplot_2d), ...
    '-', 'Color', c_ext, 'LineWidth', 0.7);

plot(t_d_2d(1:Nplot_2d), E_flow_true_2d(1:Nplot_2d), ...
    '--', 'Color', c_ext,  'LineWidth', 0.7);


xlabel('Time [days]')
ylabel('Energy change [GJ]')
% title('Total energy balance: internal energy vs. cumulative flow enthalpy')

legend( ...
    '\DeltaE_{sys} - MATLAB Baseline', ...
    'E_{flow} - MATLAB Baseline', ...
    '\DeltaE_{sys} - MATLAB Extended', ...
    'E_{flow} - MATLAB Extended', ...
    'Location', 'south',...
    'LineWidth', 0.3);

%%%------------------------------------------%%%
% 07. Export
%%%------------------------------------------%%%
version = 'v2';
dateTag = datestr(now, 'yyyymmdd');
baseName = sprintf('%s_energy_balance_%s', dateTag, version);

out_png = fullfile(RESULTS_ENERGY, [baseName '.png']);
out_pdf = fullfile(RESULTS_ENERGY, [baseName '.pdf']);

exportgraphics(fig, out_png, 'Resolution', export_png_dpi);
exportgraphics(fig, out_pdf, 'ContentType', 'vector', 'Resolution', export_png_dpi);

disp(['Saved PNG: ' out_png]);
disp(['Saved PDF: ' out_pdf]);
%%
%%%------------------------------------------%%%
% 08. Plot: residuals
%%%------------------------------------------%%%

res_true_1d   = E_flow_true_1d   - E_sys_1d;
res_simple_1d = E_flow_simple_1d - E_sys_1d;

res_true_2d   = E_flow_true_2d   - E_sys_2d;
res_simple_2d = E_flow_simple_2d - E_sys_2d;


%%%------------------------------------------%%%
% 08b. Residual in percent
%%%------------------------------------------%%%

E_ref_1d = max(abs(E_flow_true_1d));
E_ref_2d = max(abs(E_flow_true_2d));

res_pct_1d = 100 * res_true_1d / E_ref_1d;
res_pct_2d = 100 * res_true_2d / E_ref_2d;

width_cm  = 14;
height_cm = 9;

fig_res = figure( ...
    'Color','w', ...
    'Name','Energy balance residuals', ...
    'Units','centimeters', ...
    'Position',[2 2 width_cm height_cm]);

set(fig_res, 'Renderer', 'painters');
set(fig_res, 'PaperUnits', 'centimeters');
set(fig_res, 'PaperSize', [width_cm height_cm]);
set(fig_res, 'PaperPosition', [0 0 width_cm height_cm]);

ax = gca;
set(ax,'FontName', 'Times New Roman');
set(ax,'FontSize',9);
set(ax,'LineWidth',0.5);

set(gca,'GridAlpha',0.10)             % lighter grid
set(gca,'MinorGridAlpha',0.05)

hold on
grid on
box on

plot(t_d_1d(1:Nplot_1d), res_pct_1d(1:Nplot_1d), ...
    '-', 'Color', c_base,  'LineWidth', 0.7);

plot(t_d_2d(1:Nplot_2d), res_pct_2d(1:Nplot_2d), ...
    '-','Color', c_ext,  'LineWidth', 0.7);

xlabel('Time [days]')
ylabel('Energy residual [%]')
% title('Energy balance residual: flow enthalpy minus internal energy change')

legend( ...
    'E_{flow} - \DeltaE_{sys}  | MATLAB Baseline', ...
    'E_{flow} - \DeltaE_{sys}  | MATLAB Extended', ...
    'Location', 'northwest',...
    'LineWidth', 0.3);

%%%------------------------------------------%%%
% 09. Export Residual
%%%------------------------------------------%%%
version = 'v2';
dateTag = datestr(now, 'yyyymmdd');
baseName = sprintf('%s_energy_balance_residual_%s', dateTag, version);

out_png = fullfile(RESULTS_ENERGY, [baseName '.png']);
out_pdf = fullfile(RESULTS_ENERGY, [baseName '.pdf']);

exportgraphics(fig_res, out_png, 'Resolution', export_png_dpi);
exportgraphics(fig_res, out_pdf, 'ContentType', 'vector', 'Resolution', export_png_dpi);

disp(['Saved PNG: ' out_png]);
disp(['Saved PDF: ' out_pdf]);

%%
%%%------------------------------------------%%%
% 10. Residual diagnostics
%%%------------------------------------------%%%

fprintf('\n============================================================\n');
fprintf('Energy balance diagnostics\n');
fprintf('============================================================\n');

fprintf('\n--- MATLAB 1D ---\n');
fprintf('True residual   : min = %.3e J | max = %.3e J | end = %.3e J\n', ...
    min(res_true_1d), max(res_true_1d), res_true_1d(end));
fprintf('Simple residual : min = %.3e J | max = %.3e J | end = %.3e J\n', ...
    min(res_simple_1d), max(res_simple_1d), res_simple_1d(end));

fprintf('\n--- MATLAB 2D piston ---\n');
fprintf('True residual   : min = %.3e J | max = %.3e J | end = %.3e J\n', ...
    min(res_true_2d), max(res_true_2d), res_true_2d(end));
fprintf('Simple residual : min = %.3e J | max = %.3e J | end = %.3e J\n', ...
    min(res_simple_2d), max(res_simple_2d), res_simple_2d(end));

%%
%%% ============================================================ %%%
%                           LOCAL FUNCTIONS
%%% ============================================================ %%%

function runData = unpack_legacy_run(dataStruct, label)
% UNPACK_LEGACY_RUN
% Map legacy 1D result file to a unified internal structure.
%
% INPUT:
%   dataStruct - loaded .mat structure with legacy field names
%   label      - descriptive label for plotting
%
% OUTPUT:
%   runData.Res_900 - quarter-hour result matrix
%   runData.T_W_900 - water temperature field [top-to-bottom x time]
%   runData.label   - plot label

    runData = struct();

    runData.label   = label;
    runData.Res_900 = dataStruct.Res_900_d18_18;
    runData.T_W_900 = dataStruct.Res_Wasser_d18_18;

    % If required, flip the legacy water field here:
    % runData.T_W_900 = flipud(runData.T_W_900);
end

function runData = unpack_resout_run(dataStruct, label)
% UNPACK_RESOUT_RUN
% Map newer ResOut-based result file to a unified internal structure.
%
% INPUT:
%   dataStruct - loaded .mat structure with ResOut field
%   label      - descriptive label for plotting
%
% OUTPUT:
%   runData.Res_900 - quarter-hour result matrix
%   runData.T_W_900 - water temperature field [top-to-bottom x time]
%   runData.label   - plot label

    runData = struct();

    runData.label   = label;
    runData.Res_900 = dataStruct.ResOut.res.series_900;
    runData.T_W_900 = double(dataStruct.ResOut.temperature.water);
end

function [E_sys, E_flow_true, E_flow_simple] = ...
    compute_energy_balance_true(Res_900, T_W_900, rho, cp, A, Tin_hot, Tin_cold, dt)
% COMPUTE_ENERGY_BALANCE_TRUE
% ---------------------------------------------------------------
% Computes:
%   - change of total internal system energy
%   - cumulative flow enthalpy using actual outlet temperatures
%   - cumulative simple reference enthalpy using constant 80/45 °C
%
% INPUT:
%   Res_900   - result matrix
%   T_W_900   - water temperature field [top-to-bottom x time]
%   rho, cp   - water material properties
%   A         - storage cross section [m^2]
%   Tin_hot   - hot inlet temperature during charging [°C]
%   Tin_cold  - cold inlet temperature during discharging [°C]
%   dt        - time step [s]
%
% OUTPUT:
%   E_sys         - change of total internal energy [J]
%   E_flow_true   - cumulative flow enthalpy with actual outlet temperature [J]
%   E_flow_simple - cumulative simple reference enthalpy [J]
%
% CONVENTION:
%   v < 0 : charging
%           inflow at top    = Tin_hot
%           outflow at bottom= T_bottom
%
%   v > 0 : discharging
%           inflow at bottom = Tin_cold
%           outflow at top   = T_top

    % --- Change of total internal system energy ---
    E_abs = sum(Res_900(1:6,:), 1);
    E_sys = (E_abs - E_abs(1))*1e-9;

    % --- Signed flow velocity in storage ---
    v = Res_900(9,2:end) + Res_900(10,2:end);
    
    % --- Stored water temperatures ---
    T_top    = T_W_900(1,   :);
    T_bottom = T_W_900(end, :);

    % Match time dimension to flow vector
    n_steps = length(v);

    if size(T_W_900, 2) < n_steps
        error('T_W_900 has fewer time columns than required by Res_900.');
    end

    T_top    = T_top(1:n_steps);
    T_bottom = T_bottom(1:n_steps);

    % --- Stepwise flow enthalpy increments ---
    dE_true   = zeros(1, n_steps);
    dE_simple = zeros(1, n_steps);

    for k = 1:n_steps

        v_k = v(k);

        if abs(v_k) < eps
            continue
        end

        m_dot = rho * abs(v_k) * A;   % [kg/s]
        if v_k < 0
            % Charging:
            % hot water enters at top, water leaves at bottom
            dE_true(k)   = m_dot * cp * (Tin_hot - T_bottom(k)) * dt;
            dE_simple(k) = m_dot * cp * (Tin_hot - Tin_cold ) * dt;

        elseif v_k > 0
            % Discharging:
            % cold water enters at bottom, water leaves at top
            dE_true(k)   = -m_dot * cp * (T_top(k) - Tin_cold) * dt;
            dE_simple(k) = -m_dot * cp * (Tin_hot - Tin_cold) * dt;
        end
    end

    E_flow_true   = [0, cumsum(dE_true)]*1e-9;
    E_flow_simple = [0, cumsum(dE_simple)]*1e-9;
end