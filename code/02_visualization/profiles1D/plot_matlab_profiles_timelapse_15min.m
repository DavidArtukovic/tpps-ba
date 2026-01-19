clc
clear
close all
%%%------------------------------------------%%%
%  plot_matlab_profiles_timelapse_15min.m
%
%  SUMMARY:
%   Time-lapse visualization of 1D temperature profiles (15 min resolution)
%   for the same MATLAB result files as in plot_matlab_profiles.m.
%
%  NOTES:
%   - Uses the physical z-axis convention: top at +0.665 m, bottom at -36.665 m.
%   - Visualizes piston position as a semi-transparent rectangle.
%   - Exports an MP4 video into RESULTS_BASE/02_visualizations/timelapse.
%%%------------------------------------------%%%

%%%------------------------------------------%%%
% 01. Load paths and data
%%%------------------------------------------%%%

run(fullfile('..','..','..', 'configs', 'paths_local.m'));

RESULTS_VIS_BASE  = fullfile(RESULTS_BASE, "02_visualizations");
RESULTS_TIMELAPSE = fullfile(RESULTS_VIS_BASE, "timelapse");
if ~exist(RESULTS_TIMELAPSE, "dir")
    mkdir(RESULTS_TIMELAPSE);
end

DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_OWN  = fullfile(DATA_BASE, 'scenario1_freeConv');

% Scenario/time information
load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));

% Piston information (only for plotting the piston)
load(fullfile(DATA_BASE, 'scenario1/SzenarioComsol.mat'));

% MATLAB results (same as plot_matlab_profiles.m)
data_no_fk = load(fullfile(DATA_SCEN1_BASE, 'd18_h18_Res_Matlab_d18_18.mat'));
data_fk_v11 = load(fullfile(DATA_SCEN1_OWN, '260116_d18_h18_Res_Matlab_FK_v11.mat'));
data_fk_v13 = load(fullfile(DATA_SCEN1_OWN, '260117_d18_h18_Res_Matlab_FK_v13.mat'));

%%%------------------------------------------%%%
% 02. Spatial grid (MATLAB only)
%%%------------------------------------------%%%

dz_M = 0.005;  % MATLAB resolution
Nz   = size(data_no_fk.Res_System_d18_18,1) - 400;

z_M_star = flip((0:Nz-1)' * dz_M);
% Physical coordinate: top at +0.665 m, bottom at -36.665 m
z_M = -1*(37.33 - z_M_star) + 0.665;

% Inlet index (from your static script)
inlet_idx = 4288;
z_inlet   = z_M(inlet_idx);

%%%------------------------------------------%%%
% 03. Extract full 15-min resolution
%%%------------------------------------------%%%

T_noFK = data_no_fk.Res_System_d18_18(401:end,:);
T_fk11 = data_fk_v11.Res_System_d18_18_FK(401:end,:);
T_fk13 = data_fk_v13.Res_System_d18_18_FK(401:end,:);

n_steps = size(T_noFK,2);

%%%------------------------------------------%%%
% 04. Time axis (15 min per step)
%%%------------------------------------------%%%

dt_min  = 15;
t_hours = (0:n_steps-1) * dt_min / 60;

%%%------------------------------------------%%%
% 05. Piston position per time step
%%%------------------------------------------%%%
% Physical piston position (lower edge) from your mapping
% z_piston = t_900(2,i)*(-19) + (1-t_900(2,i))*(-36)
% Make sure indices are valid.
if size(t_900,2) < n_steps
    warning('t_900 has fewer time steps (%d) than temperature data (%d). Using min length.', size(t_900,2), n_steps);
    n_steps = min(size(t_900,2), n_steps);
    T_noFK = T_noFK(:,1:n_steps);
    T_fk11 = T_fk11(:,1:n_steps);
    T_fk13 = T_fk13(:,1:n_steps);
    t_hours = t_hours(1:n_steps);
end

z_piston_ts = t_900(2,1:n_steps) * (-19) + (1 - t_900(2,1:n_steps)) * (-36);

%%%------------------------------------------%%%
% 06. Model container
%%%------------------------------------------%%%

models(1).name  = 'MATLAB (no FK)';
models(1).T     = T_noFK;
models(1).z     = z_M;
models(1).style = '-';
models(1).color = [0.8500 0.3250 0.0980];

models(2).name  = 'MATLAB (FK v11)';
models(2).T     = T_fk11;
models(2).z     = z_M;
models(2).style = '--';
models(2).color = [0.4660 0.6740 0.1880];

models(3).name  = 'MATLAB (FK v13)';
models(3).T     = T_fk13;
models(3).z     = z_M;
models(3).style = '--';
models(3).color = [0.4940 0.1840 0.5560];

%%%------------------------------------------%%%
% 07. Video writer
%%%------------------------------------------%%%

videoname = fullfile(RESULTS_TIMELAPSE, 'timelapse_1D_temperature_profiles_15min_matlab_only.mp4');

v = VideoWriter(videoname,'MPEG-4');
% Keep it readable; increase if you want a faster video
v.FrameRate = 10;
open(v);

%%%------------------------------------------%%%
% 08. Figure initialization
%%%------------------------------------------%%%

figure('Color','w','Name','1D temperature profiles – time-lapse (15 min)');
hold on
grid on

for mm = 1:numel(models)
    h(mm) = plot(models(mm).T(:,1), models(mm).z, ...
        models(mm).style, ...
        'Color', models(mm).color, ...
        'LineWidth',1.5);
end

% Inlet line (constant)
h_inlet = yline(z_inlet, 'k--', 'LineWidth', 1);

xlabel('Temperature [$^{\circ}$C]','Interpreter','Latex')
ylabel('System depth [m]','Interpreter','Latex')
legend({models.name},'Interpreter','Latex','Location','best')

xlim([40 80])
ylim([min(z_M) max(z_M)])

% Piston rectangle (initialized)
h_piston = 18;   % [m]
z_piston = z_piston_ts(1);

h_rect = rectangle('Position',[40 z_piston 40 h_piston], ...
                   'FaceColor',[1 1 1], 'FaceAlpha',0.12, ...
                   'EdgeColor',[0 0 0]);

% Put rectangle behind lines
uistack(h_rect,'bottom');

% Title handle
h_title = title('', 'Interpreter','Latex');

%%%------------------------------------------%%%
% 09. Animation loop
%%%------------------------------------------%%%

for i = 1:n_steps

    for mm = 1:numel(models)
        set(h(mm), 'XData', models(mm).T(:,i), ...
                   'YData', models(mm).z);
    end

    z_piston = z_piston_ts(i);
    set(h_rect, 'Position', [40 z_piston 40 h_piston]);

    h_title.String = sprintf('1D temperature profiles (15 min) - t = %.2f h', t_hours(i));

    drawnow
    writeVideo(v, getframe(gcf));
end

close(v);

disp(['Saved video: ' videoname]);
