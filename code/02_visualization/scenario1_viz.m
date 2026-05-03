clc
clear
close

%%%------------------------------------------%%%
% 01. Load paths and data
%%%------------------------------------------%%%

run(fullfile('..','..', 'configs', 'paths_local.m'));

RESULTS_VIS_BASE  = fullfile(RESULTS_BASE, "02_visualizations");

DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_INIT = fullfile(DATA_BASE, 'init');

load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));
load(fullfile(DATA_BASE, 'scenario1/SzenarioComsol.mat'));

load(fullfile(DATA_INIT, '20260330_d18_hp18.0_gap0.5_2D_chunk009.mat'));       % Actual init with modified flux

Lflow = 4.9236e-05; % m/s
Dflow = 3.7395e-05; % m/s

L_mflow = Lflow*pi*9^2*1000;
D_mflow = Dflow*pi*9^2*1000;
%%
%%%------------------------------------------%%%
% Scenario plot
%%%------------------------------------------%%%

t_h = (0:size(t_900,2)-1) * 0.25;   % 15 min steps in hours
t_days = t_h / 24;
% Example: convert operation flag to mass flow
op_flag = t_900(1,:);      % 1 = charging, -1 = discharging, 0 = idle
piston_flag = t_900(2,:);  % 1 = top, 0 = bottom

mdot = zeros(size(op_flag));

% Example assignment:
% mflow_charge must be known from your scenario definition
% mflow_discharge can be set analogously if available
mdot(op_flag == 1)  = L_mflow;
mdot(op_flag == -1) = -D_mflow;
mdot(op_flag == 0)  = 0;

% Piston top position
z_piston_top = (piston_flag-1) *16;
%%

figure('Color','w', 'Name', 'Operating scenario');

fig = gcf;

% Size in centimeters
fig.Units = 'centimeters';
fig.Position = [5 5 14 8];   % [x y width height]

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');


%%%------------------------------------------%%%
% 02 Mass flow
%%%------------------------------------------%%%
ax1 = nexttile;
stairs(ax1, t_days, mdot, 'k', 'LineWidth',1.0)
grid(ax1, 'on');
ax1.Box = 'off';
set(ax1, 'FontSize',10);
ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';
ax1.TickDir = 'out';
ax1.XTickLabel = [];
ax1.GridAlpha = 0.15;
ylabel(ax1, '$\dot{m}$ [kg\,s$^{-1}$]','Interpreter','latex')
ylim(ax1,[-13 17]);


%%%------------------------------------------%%%
% 03 Piston position
%%%------------------------------------------%%%
ax2 = nexttile;
stairs(ax2, t_days, z_piston_top, 'k', 'LineWidth',1.0)
grid(ax2, 'on');
set(ax2, 'FontSize',10);
ax2.Box = 'off';
ax2.TickDir = 'out';
ax2.XAxisLocation = 'bottom';
ax2.YAxisLocation = 'left';
ax2.GridAlpha = 0.15;
ylabel(ax2,'$z_{\mathrm{piston,top}}$ [m]','Interpreter','latex')
xlabel(ax2, 'Time [days]')
ylim(ax2,[-17 1])


%%%------------------------------------------%%%
% 04 Export
%%%------------------------------------------%%%
version = 'v4';
dateTag = datestr(now, 'yyyymmdd');
export_png_dpi = 900;
baseName = sprintf('%s_scenario1_operating_schema_%s', dateTag, version);

out_pdf = fullfile(RESULTS_VIS_BASE, [baseName '.pdf']);

exportgraphics(fig, out_pdf, 'ContentType', 'vector', 'Resolution', export_png_dpi);

disp(['Saved PDF: ' out_pdf]);