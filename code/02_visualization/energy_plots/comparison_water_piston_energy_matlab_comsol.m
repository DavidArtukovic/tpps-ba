%% comparison_water_piston_energy_matlab_comsol.m
% ---------------------------------------------------------------
% PURPOSE:
%   Compare the thermal energy evolution of the water domain and
%   the piston domain for three model variants:
%       (1) MATLAB 1D baseline
%       (2) COMSOL reference
%       (3) MATLAB 2D piston extended model
%
% DESCRIPTION:
%   - Loads MATLAB and COMSOL scenario results.
%   - Extracts water and piston temperature fields for all models.
%   - Computes the total thermal energy change relative to t = 0
%     separately for:
%         (a) water
%         (b) piston
%   - Uses a common total domain volume for each physical domain
%     to make the comparison independent of discretization.
%   - Plots:
%         solid lines   -> water energy change
%         dashed lines  -> piston energy change
%     with:
%         left y-axis   -> water energy [GJ]
%         right y-axis  -> piston energy [GJ]
%
% NOTES:
%   - The comparison is based on total domain energies, not local
%     or moving sub-volume energies.
%   - All energies are plotted as changes relative to the first
%     stored time step.
%   - Different spatial discretizations are handled by assigning
%     each temperature node an equivalent cell volume such that
%     the total physical domain volume is preserved.
% ---------------------------------------------------------------

clc
clear
close all

%%%------------------------------------------%%%
% 01. Load paths and data
%%%------------------------------------------%%%

run(fullfile('..','..','..', 'configs', 'paths_local.m'));

DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_NOFK = fullfile(DATA_BASE, 'scenario1');
DATA_SCEN1_FK   = fullfile(DATA_BASE, 'scenario1_freeConv');
DATA_INIT       = fullfile(DATA_BASE, 'init');

% Piston position information 
load(fullfile(DATA_BASE, 'scenario1/SzenarioComsol.mat'));

% --- Init / geometry data
load(fullfile(DATA_INIT, '20260330_d18_hp18.0_gap0.5_2D_chunk009.mat'));

% --- MATLAB 1D baseline
data_1 = load(fullfile(DATA_SCEN1_NOFK, ...
    'd18_h18_Res_Matlab_d18_18.mat'));

% --- COMSOL reference
data_2 = load(fullfile(DATA_SCEN1_BASE, ...
    'T1_1818_900.mat'));

% --- MATLAB 2D piston extended
data_3 = load(fullfile(DATA_SCEN1_FK, ...
    '260331_d18_h18_Res_Matlab_FK_2d_v12.mat'));

%%
%%%------------------------------------------%%%
% 02. Basic settings and material properties
%%%------------------------------------------%%%

dz_M = 0.005;                     % MATLAB axial resolution [m]
dz_C = 0.05;                      % COMSOL axial resolution [m]
dt   = 900;                       % 15 min [s]

% --- Time step count (use common length for all models)
n_steps = 1921;

% --- Time axis
t_days = (0:n_steps-1) * dt / 86400;

% --- Material properties
rho_w = InitOut.param.SW(1,2);    % water density [kg/m^3]
cp_w  = InitOut.param.SW(1,3);    % water heat capacity [J/(kg K)]

rho_p = InitOut.param.SW(2,2);    % piston density [kg/m^3]
cp_p  = InitOut.param.SW(2,3);    % piston heat capacity [J/(kg K)]

%%
%%%------------------------------------------%%%
% 03. Extract temperature from Comsol
%%%------------------------------------------%%%

T_W_C = double(data_2.W(:,1:n_steps));
T_P_C = double(data_2.P(:,1:n_steps));

% First and last node in Comsol have NaN temperatues
% set neighbouring temperature as replacement
% some columns in the middle have also single missing values
% interpolate them.
T_W_C_interp = fillmissing(T_W_C, 'linear', 1, 'EndValues', 'nearest');
%%
%%%------------------------------------------%%%
% 04. Domain volumes
%%%------------------------------------------%%%

% --- Storage geometry
r_ST   = InitOut.geom.d(1) / 2;             % storage radius [m]
gap    = InitOut.meta.r_gap;                 % assumed ring-gap thickness [m]
r_pist = r_ST - gap;                        % piston radius [m]
h_pist = 18.0;                              % piston height [m]

h_kupp = 1; % height of dome

% Piston lower edge from scenario mapping
z_piston_ts = t_900(2,1:n_steps-1) * (-19) + (1 - t_900(2,1:n_steps-1)) * (-35);
z_piston_ts = [-19 z_piston_ts];
%%
%%%------------------------------------------%%%
% 05. Compute thermal energy change
%%%------------------------------------------%%%

% --- Water energy change [GJ]
dE_W_1 = (data_1.Res_900_d18_18(2,1:end)-data_1.Res_900_d18_18(2,1))*1e-9;

% --- COMSOL water energy [GJ] ---
dE_W_C = compute_energy_uniform( ...
    T_W_C_interp, z_piston_ts, ...
    rho_w, cp_w, dz_C, ...
    r_ST, r_pist, h_kupp, h_pist);

dE_W_3 = (data_3.ResOut.res.series_900(2,1:n_steps)...
        -data_3.ResOut.res.series_900(2,1))*1e-9;

% --- Piston energy change [GJ]
dE_P_1 = (data_1.Res_900_d18_18(3,1:end)-data_1.Res_900_d18_18(3,1))*1e-9;

% --- COMSOL piston energy [GJ] ---
dE_P_C = compute_energy_piston_2d( ...
    T_P_C, ...
    rho_p, cp_p, dz_C, ...
    r_pist);

dE_P_3 = (data_3.ResOut.res.series_900(3,1:n_steps)...
        -data_3.ResOut.res.series_900(3,1))*1e-9;
%%
%%%------------------------------------------%%%
% 06. Plot
%%%------------------------------------------%%%

% Academic / muted colors
c_base = [0.10 0.10 0.10];   % dark gray
c_cmsl = [0.00 0.45 0.74];   % muted blue
c_ext  = [0 0.39 0];         % muted green

figure('Color','w', ...
    'Name','Water and piston energy comparison: MATLAB vs COMSOL', ...
    'Position',[100 100 1300 650]);

hold on
box on
grid on

ax = gca;

set(gcf, 'Color', 'w');
set(ax, 'Color', 'w');

ax.LineWidth = 1.0;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 12;

ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridAlpha = 0.18;
ax.MinorGridAlpha = 0.10;
ax.GridColor = [0 0 0];

ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

% ax.TickDir = 'in';
% ax.TickLength = [0.015 0.015];

ax.Layer = 'top';

% --- Left axis: water energy
yyaxis left
ax = gca;
ax.YColor = 'k';   % left y-axis black

hW1 = plot(t_days, dE_W_1, '-', ...
    'LineWidth', 1.8, 'Color', c_base);
hold on

hWC = plot(t_days, dE_W_C, '-', ...
    'LineWidth', 1.6, 'Color', c_cmsl);

hW3 = plot(t_days, dE_W_3, '-', ...
    'LineWidth', 1.6, 'Color', c_ext);

ylabel('\Delta E_{water} [GJ]')

% --- Right axis: piston energy
yyaxis right
ax.YColor = 'k';   % right y-axis black
ax.XColor = 'k';   % x-axis black

hP1 = plot(t_days, dE_P_1, '--', ...
    'LineWidth', 1.8, 'Color', c_base);
hold on

hPC = plot(t_days, dE_P_C, '--', ...
    'LineWidth', 1.6, 'Color', c_cmsl);

hP3 = plot(t_days, dE_P_3, '--', ...
    'LineWidth', 1.6, 'Color', c_ext);

ylabel('\Delta E_{piston} [GJ]')

xlabel('Time [days]')

legend([hW1 hWC hW3 hP1 hPC hP3], ...
    'Water - MATLAB 1D baseline', ...
    'Water - COMSOL', ...
    'Water - MATLAB 2D extended', ...
    'Piston - MATLAB 1D baseline', ...
    'Piston - COMSOL', ...
    'Piston - MATLAB 2D extended', ...
    'Location', 'south');

lgd.Box = 'off';
% No title -> use LaTeX caption instead

%%
%%%------------------------------------------%%%
% 07. Diagnostics
%%%------------------------------------------%%%

fprintf('\n============================================================\n');
fprintf('Partial energy comparison diagnostics\n');
fprintf('============================================================\n');

fprintf('\n--- Final water energy change [GJ] ---\n');
fprintf('MATLAB 1D baseline : %.6f\n', dE_W_1(end));
fprintf('COMSOL             : %.6f\n', dE_W_C(end));
fprintf('MATLAB 2D extended : %.6f\n', dE_W_3(end));

fprintf('\n--- Final piston energy change [GJ] ---\n');
fprintf('MATLAB 1D baseline : %.6f\n', dE_P_1(end));
fprintf('COMSOL             : %.6f\n', dE_P_C(end));
fprintf('MATLAB 2D extended : %.6f\n', dE_P_3(end));

%%
%%% ============================================================ %%%
%                           LOCAL FUNCTIONS
%%% ============================================================ %%%

function dE_GJ = compute_energy_uniform(T_W, z_piston_ts, rho, cp, dz, ...
    r_ST, r_pist, h_kupp, h_pist)
% COMPUTE_ENERGY_UNIFORM
% ---------------------------------------------------------------
% Compute the COMSOL water energy change relative to the first time
% step using variable node volumes:
%   - top dome
%   - bottom dome
%   - full cylinder cross section outside piston height
%   - ring-gap cross section at piston height
%
% INPUT:
%   T_W         - COMSOL water temperatures [Nz x Nt]
%   z_piston_ts - piston lower edge in physical coordinates [1 x Nt]
%   rho, cp     - water properties
%   dz          - COMSOL axial spacing [m]
%   r_ST        - storage radius [m]
%   r_pist      - piston radius [m]
%   h_kupp      - dome height [m]
%   h_pist      - piston height [m]
%
% OUTPUT:
%   dE_GJ       - water energy change relative to t = 0 [GJ]
%
% NOTES:
%   - Uses the original COMSOL 38 m axis:
%         0 ... 1 m   : top dome
%         1 ... 37 m  : cylindrical part
%         37 ... 38 m : bottom dome
%   - The piston occupies an 18 m interval inside the cylindrical part.
%   - On piston height, the water cross section is reduced to the ring gap.

    T_W = double(T_W);

    [n_z, n_t] = size(T_W);
    z_raw = (0:n_z-1)' * dz;    % 0 ... 38 m from top to bottom

    % --- Constant cross sections ---
    A_full = pi * r_ST^2;
    A_gap  = pi * (r_ST^2 - r_pist^2);

    E_abs = zeros(1, n_t);

    for k = 1:n_t

        % --- Build cross section vector for current piston position ---
        A_vec = zeros(n_z, 1);

        % Convert physical piston position to original COMSOL raw axis.
        % In the cylindrical middle region the mapping is:
        %   z_raw = 1 - z_phys
        z_pist_low_raw = 1 - z_piston_ts(k);         % lower piston edge
        z_pist_top_raw = z_pist_low_raw - h_pist;    % upper piston edge

        for i = 1:n_z

            zi = z_raw(i);

            if zi < h_kupp
                % --- Top dome: paraboloid-like cap ---
                % x = distance from top tip downward, 0 ... h_kupp
                x = zi;
                A_vec(i) = dome_area_paraboloid(x, r_ST, h_kupp);

            elseif zi > (38 - h_kupp)
                % --- Bottom dome: mirrored paraboloid-like cap ---
                % x = distance from bottom tip upward, 0 ... h_kupp
                x = 38 - zi;
                A_vec(i) = dome_area_paraboloid(x, r_ST, h_kupp);

            else
                % --- Cylindrical part ---
                if zi >= z_pist_top_raw && zi <= z_pist_low_raw
                    A_vec(i) = A_gap;
                else
                    A_vec(i) = A_full;
                end
            end
        end

        % --- Node volumes and energy ---
        dV_vec = A_vec * dz;
        E_abs(k) = rho * cp * sum(T_W(:,k) .* dV_vec);
    end

    dE_GJ = (E_abs - E_abs(1)) * 1e-9;
end

function dE_GJ = compute_energy_piston_2d(T_P, rho, cp, dz, r_pist)
% COMPUTE_ENERGY_PISTON_2D
% ---------------------------------------------------------------
% Compute the COMSOL piston energy change relative to the first time
% step using a constant piston cross section.
%
% INPUT:
%   T_P    - COMSOL piston temperatures [Nz x Nt]
%   rho, cp- piston material properties
%   dz     - COMSOL axial spacing [m]
%   r_pist - piston radius [m]
%
% OUTPUT:
%   dE_GJ  - piston energy change relative to t = 0 [GJ]

    T_P = double(T_P);

    A_pist = pi * r_pist^2;
    dV = A_pist * dz;

    E_abs = rho * cp * dV * sum(T_P, 1);
    dE_GJ = (E_abs - E_abs(1)) * 1e-9;
end

function A = dome_area_paraboloid(x, r, h)
% DOME_AREA_PARABOLOID
% ---------------------------------------------------------------
% Axial cross-sectional area of a dome with volume
%   V = (2/3) * pi * r^2 * h
% using a paraboloid-like profile.
%
% INPUT:
%   x - axial coordinate from dome tip toward cylinder [m]
%   r - storage radius at cylinder junction [m]
%   h - dome height [m]
%
% OUTPUT:
%   A - local cross-sectional area [m^2]

    xi = max(0, min(x / h, 1));
    A = pi * r^2 * (2 * xi - xi^2);
end