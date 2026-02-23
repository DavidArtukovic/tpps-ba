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


%%
DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_NOFK  = fullfile(DATA_BASE, 'scenario1');
DATA_SCEN1_FK  = fullfile(DATA_BASE, 'scenario1_freeConv');


% Scenario/time information
load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));

% Piston information (only for plotting the piston)
load(fullfile(DATA_BASE, 'scenario1/SzenarioComsol.mat'));

% MATLAB results (same as plot_matlab_profiles.m)
data_1 = load(fullfile(DATA_SCEN1_NOFK, 'd18_h18_Res_Matlab_d18_18.mat'));
data_2 = load(fullfile(DATA_SCEN1_FK, '260207_d18_h18_Res_Matlab_FK_v18.mat'));
data_3 = load(fullfile(DATA_SCEN1_FK, '260219_d18_h18_Res_Matlab_FK_v20.mat'));

%%
%%%------------------------------------------%%%
% 02. Spatial grid
%%%------------------------------------------%%%

dz_M = 0.005;  % MATLAB resolution
Nz   = size(data_1.Res_System_d18_18,1) - 400;

z_M_star = flip((0:Nz-1)' * dz_M);
% Physical coordinate: top at +0.665 m, bottom at -36.665 m
z_M = -1*(37.33 - z_M_star) + 0.665;

% --- geometric bypass positions (relative to z = 0) ---
z_ps_up_geom  = -1 - 16/2;
z_ps_low_geom = -3*1 - 0.05 - 16;

% --- corresponding indices on z_M grid ---
[~, idx_bypass_upper] = min(abs(z_M - z_ps_up_geom));
[~, idx_bypass_lower] = min(abs(z_M - z_ps_low_geom));


%%%------------------------------------------%%%
% 03. Extract full 15-min resolution
%%%------------------------------------------%%%

T_1 = data_1.Res_System_d18_18(401:end,:);
T_2 = data_2.ResOut.temperature.system(401:end,:);
T_3 = data_3.ResOut.temperature.system(401:end,:);

n_steps = size(T_1, 2);
%%
%-----------------------------------------
% 2.1 Extract Ring-Gap Info
%------------------------------------------

% Part for Ringgap
blockLen = 389;
startIdx = data_1.Res_900_d18_18(8,2:n_steps+1)-1;   % [1 x Nt]
% --- build row index matrix ---
rowOffsets = (0:blockLen-1)';                 % [389 x 1]
rows = startIdx(:)' + rowOffsets;         % [389 x Nt_reduced]
colsMat = repmat(1:n_steps, blockLen, 1);       % [389 x Nt_reduced]

%%% No FK %%%
A        = data_1.Res_Wasser_d18_18(:,1:n_steps);     % [Nz x Nt]
linInd = sub2ind(size(A), rows, colsMat);
T_ring_1 = A(linInd);                % [389 x Nt_reduced]

%%% FK %%%
A        = data_2.ResOut.temperature.water(:,1:n_steps);     % [Nz x Nt]
linInd = sub2ind(size(A), rows, colsMat);
T_ring_2 = A(linInd);                % [389 x Nt_reduced]

%%% FK %%%
A        = data_3.ResOut.temperature.water(:,1:n_steps);     % [Nz x Nt]
linInd = sub2ind(size(A), rows, colsMat);
T_ring_3 = A(linInd);                % [389 x Nt_reduced]

clear A linInd startIdx rowOffsets rows

%%
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
    T_1 = T_noFK(:,1:n_steps);
    T_2 = T_fk(:,1:n_steps);
    t_hours = t_hours(1:n_steps);
end

z_piston_ts = t_900(2,1:n_steps) * (-19) + (1 - t_900(2,1:n_steps)) * (-36);
%%
%%%------------------------------------------%%%
% 06. Model container
%%%------------------------------------------%%%

models(1).name   = 'no FK';
models(1).T_sys  = T_1;
models(1).T_ring = T_ring_1;
models(1).z      = z_M;
models(1).style  = '-';
models(1).color  = [0 0 0];          % black

models(2).name   = 'FK ';
models(2).T_sys  = T_2;
models(2).T_ring = T_ring_2;
models(2).z      = z_M;
models(2).style  = '-';
models(2).color  = [0.0 0.7 0.0];    % green

models(3).name   = 'FK - adapted diffusion + radial coupling';
models(3).T_sys  = T_3;
models(3).T_ring = T_ring_3;
models(3).z      = z_M;
models(3).style  = '-';
models(3).color  = [0.6 0.0 0.0];    % blue

z_ring_norm = linspace(0,1,blockLen)';   % normalized vertical coordinate

%%%------------------------------------------%%%
% 07. Video writer
%%%------------------------------------------%%%
dateTag = datestr(now,'yyyymmdd');
videoname = fullfile(RESULTS_TIMELAPSE, ...
    [dateTag '_timelapse_1D_temperature_profiles_15min_matlab_only_v2.mp4']);

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

% --- right-side markers (short horizontal segments) ---
x_mark = [78 80];

h_ps_up  = plot(x_mark, [NaN NaN], 'k-', 'LineWidth',3);
h_ps_low = plot(x_mark, [NaN NaN], 'k-', 'LineWidth',3);
h_mem    = plot(x_mark, [NaN NaN], 'k-', 'LineWidth',3);

x_txt = 80.3;   % slightly right of marker lines

h_txt_up  = text(x_txt, NaN, 'upper bypass inlet', 'FontSize',9, 'Color','k', 'VerticalAlignment','middle');
h_txt_mem = text(x_txt, NaN, 'membrane',           'FontSize',9, 'Color','k', 'VerticalAlignment','middle');
h_txt_low = text(x_txt, NaN, 'lower bypass inlet', 'FontSize',9, 'Color','k', 'VerticalAlignment','middle');


for mm = 1:numel(models)
    % Water-Piston system
    h_sys(mm) = plot(models(mm).T_sys(:,1), models(mm).z, ...
    models(mm).style, ...
    'Color', models(mm).color, ...
    'LineWidth',1.5);

    % --- Ring-gap handles for model 2 and 3 ---
    if mm >= 2
        h_ring(mm) = plot(nan, nan, '--', ...
            'Color', models(mm).color, ...
            'LineWidth',1.5);
    end
end



xlabel('Temperature [$^{\circ}$C]','Interpreter','Latex')
ylabel('System depth [m]','Interpreter','Latex')
% Build legend entries (system + ring-gap)
legend_handles = [];
legend_labels  = {};

for mm = 1:numel(models)
    legend_handles(end+1) = h_sys(mm);
    legend_labels{end+1}  = [models(mm).name ' – system'];

    if mm >= 2
        legend_handles(end+1) = h_ring(mm);
        legend_labels{end+1}  = [models(mm).name ' – ring gap'];
    end
end


legend(legend_handles, legend_labels, 'Interpreter','none', 'Location','best')

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

% FK highlight rectangle (whole volume)
h_fk_rect = rectangle('Position',[40 0 40 0], ...
    'EdgeColor','none', ...
    'LineWidth',3, ...
    'FaceColor','none');

uistack(h_fk_rect,'bottom');

% Title handle
h_title = title('', 'Interpreter','Latex');
set(h_ps_up ,  'YData', [z_M(idx_bypass_upper) z_M(idx_bypass_upper)]);
set(h_ps_low, 'YData', [z_M(idx_bypass_lower) z_M(idx_bypass_lower)]);
set(h_txt_up ,  'Position', [x_txt z_M(idx_bypass_upper) 0]);
set(h_txt_low,  'Position', [x_txt z_M(idx_bypass_lower) 0]);

%%%------------------------------------------%%%
% 09. Animation loop
%%%------------------------------------------%%%

% set(h_sys(2), 'Visible', 'off');

for i = 1:n_steps

    SoC = t_900(2,i);   % state of charge in [0,1]

    z_membrane = -2*1 - 16 * (1 - SoC/2);
    [~, idx_membrane] = min(abs(z_M - z_membrane));

    set(h_mem, 'YData', [z_M(idx_membrane) z_M(idx_membrane)], 'Color', [1 0.85 0]);
    set(h_txt_mem, 'Position', [x_txt z_M(idx_membrane) 0]);

    z_piston = z_piston_ts(i);
    z_ring   = flip(z_piston + z_ring_norm * h_piston);

    % --- FK visualization (whole volume) ---
    fk_code = data_3.ResOut.fk.logging_code(i);
    
    if fk_code == 0
        set(h_fk_rect,'EdgeColor','none');
    else
        % Color: warm in = red, cold in = blue
        if abs(fk_code) == 2
            fk_color = [0.85 0.1 0.1];   % red
        else
            fk_color = [0.1 0.3 0.85];   % blue
        end
    
        if fk_code > 0
            % --- upper water volume ---
            y0 = z_M(idx_bypass_upper);
            height_fk  = abs(z_M(idx_bypass_upper) -max(z_M));
        else
            % --- lower water volume ---
            y0 = min(z_M);
            height_fk  = z_piston - min(z_M);
        end
    
        set(h_fk_rect, ...
            'Position',[40 y0 40 height_fk], ...
            'EdgeColor',fk_color);
    end


    % Load models   
    for mm = 1:numel(models)
        % System profile
        set(h_sys(mm), 'XData', models(mm).T_sys(:,i), ...
                   'YData', models(mm).z);

    end
    % Ring-gap profile mapped to piston height
    for mm = 2:numel(models)
        set(h_ring(mm), ...
            'XData', models(mm).T_ring(:,i), ...
            'YData', z_ring);
    end

    set(h_rect, 'Position', [40 z_piston 40 h_piston]);

    h_title.String = sprintf('1D temperature profiles (15 min) - t = %.2f h', t_hours(i));

    drawnow
    writeVideo(v, getframe(gcf));
end

close(v);

disp(['Saved video: ' videoname]);
