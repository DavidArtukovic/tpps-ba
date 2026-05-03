clc
clear
close all
%%%------------------------------------------%%%
%  timelapse_matlab_comsol_15min_2d.m
%
%  SUMMARY:
%   Timelapse visualization of matlab/comsol 1D temperature profiles (15 min resolution)
%   and 2D piston temperature field from Matlab, side-by-side.
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
RESULTS_TIMELAPSE = fullfile(RESULTS_VIS_BASE, "04_timelapse");

DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_NOFK  = fullfile(DATA_BASE, 'scenario1');
DATA_SCEN1_FK  = fullfile(DATA_BASE, 'scenario1_freeConv');
DATA_INIT = fullfile(DATA_BASE, 'init');

% Load init results and scenario files
load(fullfile(DATA_INIT, '20260330_d18_hp18.0_gap0.5_2D_chunk009.mat')); %2D init file
% Scenario/time information
load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));

% Piston information (only for plotting the piston)
load(fullfile(DATA_BASE, 'scenario1/SzenarioComsol.mat'));

% MATLAB and COMSOL results (same as plot_matlab_profiles.m)
data_1 = load(fullfile(DATA_SCEN1_NOFK, 'd18_h18_Res_Matlab_d18_18.mat')); % Matlab 1D no FK
data_2 = load(fullfile(DATA_SCEN1_BASE, 'T1_1818_900.mat')); % COMSOL
data_3 = load(fullfile(DATA_SCEN1_FK, '260330_d18_h18_Res_Matlab_FK_2d_v10.mat')); % Matlab 2D piston with FK


%%%------------------------------------------%%%
% 02. Spatial grid
%%%------------------------------------------%%%

dz_M = 0.005;  % MATLAB resolution
Nz   = size(data_1.Res_System_d18_18,1) - 400; % remove isolation nodes

z_M_star = flip((0:Nz-1)' * dz_M);
% Physical coordinate: top at +0.665 m, bottom at -36.665 m
z_M = -1*(37.33 - z_M_star) + 0.665;

% geometric bypass positions (relative to z = 0)
% based on Dominic#s sheet
z_ps_up_geom  = -1 - 16/2;
z_ps_low_geom = -3*1 - 0.05 - 16;

% corresponding indices on z_M grid
[~, idx_bypass_upper] = min(abs(z_M - z_ps_up_geom));
[~, idx_bypass_lower] = min(abs(z_M - z_ps_low_geom));

%%%------------------------------------------%%%
% 03. Extract full 15-min resolution
%%%------------------------------------------%%%

T_1 = data_1.Res_System_d18_18(401:end,:);
T_3 = data_3.ResOut.temperature.system(401:end,:);
T_3_p = data_3.ResOut.temperature.piston; % 2d piston

% COMSOL: drop first stored time point (t = 0)
W_C_raw = data_2.W(:,2:end);
P_C_raw = data_2.P(:,2:end);


n_steps = 375;

% truncate to n_steps
T_1     = T_1(:,1:n_steps);
T_3     = T_3(:,1:n_steps);
T_3_p   = T_3_p(:,1:n_steps);
W_C_raw = W_C_raw(:,1:n_steps);
P_C_raw = P_C_raw(:,1:n_steps);

%-----------------------------------------
% 3.1 Extract Ring-Gap Temperatures for Matlab models
%------------------------------------------

% Part for Ringgap
blockLen = 389; % hard coded ring-gap size
startIdx = data_1.Res_900_d18_18(8,2:n_steps+1)-1;   % [1 x Nt]

% --- build row index matrix ---
rowOffsets = (0:blockLen-1)';                   % [389 x 1]
rows = startIdx(:)' + rowOffsets;               % [389 x Nt_reduced]
colsMat = repmat(1:n_steps, blockLen, 1);       % [389 x Nt_reduced]

%%% No FK %%%
A        = data_1.Res_Wasser_d18_18(:,1:n_steps);     % [Nz x Nt]
linInd = sub2ind(size(A), rows, colsMat);
T_ring_1 = A(linInd);                                 % [389 x Nt_reduced]

%%% FK %%%
A        = data_3.ResOut.temperature.water(:,1:n_steps);     % [Nz x Nt]
linInd = sub2ind(size(A), rows, colsMat);
T_ring_3 = A(linInd);                                        % [389 x Nt_reduced]

clear A linInd startIdx rowOffsets rows

%%%------------------------------------------%%%
% 04. Time axis (15 min per step)
%%%------------------------------------------%%%

dt_min  = 15;
t_hours = (1:n_steps) * dt_min / 60;

%%%------------------------------------------%%%
% 04.1 Annotation events for presentation video
%%%------------------------------------------%%%

event_step = [];
event_text = {};
event_xTipB = [];
event_yTipB = [];
event_xTipE = [];
event_yTipE = [];
event_xTxt = [];
event_yTxt = [];

% 80 °C inflow at top: n_step 24 ... 44
steps = 24:44;
event_step = [event_step, steps];
event_text = [event_text, repmat({'80°C water inflow at top'}, 1, numel(steps))];
event_xTipB = [event_xTipB, 78 * ones(1,numel(steps))];
event_yTipB = [event_yTipB, 0.3 * ones(1,numel(steps))];
event_xTipE = [event_xTipE, 75.5 * ones(1,numel(steps))];
event_yTipE = [event_yTipE, -3.3 * ones(1,numel(steps))];
event_xTxt = [event_xTxt, 72.5 * ones(1,numel(steps))];
event_yTxt = [event_yTxt, -3.7 * ones(1,numel(steps))];

% faster mixing in Comsol: n_step 46 ... 60
steps = 46:60;
xTip = 55 + 0.8*(0:numel(steps)-1);
event_step = [event_step, steps];
event_text = [event_text, repmat({'faster mixing in Comsol'}, 1, numel(steps))];
event_xTipB = [event_xTipB, xTip];
event_yTipB = [event_yTipB, -3 * ones(1,numel(steps))];
event_xTipE = [event_xTipE, xTip + 3];
event_yTipE = [event_yTipE, -3 * ones(1,numel(steps))];
event_xTxt = [event_xTxt, xTip + 3];
event_yTxt = [event_yTxt, -3 * ones(1,numel(steps))];

% warm water reaches ring-gap in Matlab: n_step 56 ... 72
steps = 58:72;
xTip = [52 53 54 56 58 60 61 61 61 61 61 61 61 61 61];
yTip = [-1 -1 -1 -1 -1.5 -1.5 -2 -2.0 -2.5 -3 -3.5 -3.5 -4 -4.5 -4.7];
event_step = [event_step, steps];
event_text = [event_text, repmat({'warm water reaches ring-gap in Matlab'}, 1, numel(steps))];
event_xTipB = [event_xTipB, xTip];
event_yTipB = [event_yTipB, yTip];
event_xTipE = [event_xTipE, xTip + 2];
event_yTipE = [event_yTipE, yTip-5];
event_xTxt = [event_xTxt, xTip -1];
event_yTxt = [event_yTxt, yTip-7];


% warm water reaches lower ring-gap via bypass: n_step 128 ... 144
steps = 128:144;

xTip = [50 50.3 50.6 50.9 51.2 51.5 51.8 52.1 52.4 52.7 53 53.3 53.6 53.9 54.2 54.5 54.8];
yTip = -17 * ones(1,numel(steps));

event_step = [event_step, steps];
event_text = [event_text, repmat({'warm water reaches lower ring-gap via bypass'}, 1, numel(steps))];
event_xTipB = [event_xTipB, xTip];
event_yTipB = [event_yTipB, yTip];
event_xTipE = [event_xTipE, xTip + 3];
event_yTipE = [event_yTipE, yTip];
event_xTxt = [event_xTxt, xTip + 3.5];
event_yTxt = [event_yTxt, yTip];


% extended discrete FK update in Matlab: n_step 217 ... 252 & 312 ... 339

steps = [217:252, 312:339];

xTip = [ ...
    62 62.3 62.6 62.9 63.2 63.5 63.8 64.1 64.4 64.7 ...
    65 65.3 65.6 65.9 66.2 66.5 66.8 67.1 67.4 67.7 ...
    68 68.15 68.25 68.35 68.45 68.55 68.65 68.75 68.90 69.05 ...
    69.20 69.35 69.50 69.65 69.80 70 ...
    69 69.15 69.25 69.35 69.45 69.55 69.65 69.75 69.85 69.95 ...
    70.05 70.15 70.25 70.35 70.45 70.55 70.65 70.75 70.85 70.95 ...
    71.15 71.35 71.55 71.75 71.95 72.15 72.35 72.55];

event_step = [event_step, steps];
event_text = [event_text, repmat({'discrete FK update\newline in Matlab'}, 1, numel(steps))];

event_xTipB = [event_xTipB, xTip];
event_yTipB = [event_yTipB, -15 * ones(1,numel(steps))];

event_xTipE = [event_xTipE, xTip + 2];
event_yTipE = [event_yTipE, -15 * ones(1,numel(steps))];

event_xTxt = [event_xTxt, xTip + 2.5];
event_yTxt = [event_yTxt, -15 * ones(1,numel(steps))];


% inverse thermocline at piston bottom: n_step 272 ... 300
steps = 272:300;

event_step = [event_step, steps];
event_text = [event_text, repmat({'inverse thermocline\newline at piston bottom'}, 1, numel(steps))];
event_xTipB = [event_xTipB, 70.5 * ones(1,numel(steps))];
event_yTipB = [event_yTipB, -18.9 * ones(1,numel(steps))];
event_xTipE = [event_xTipE, 73 * ones(1,numel(steps))];
event_yTipE = [event_yTipE, -17.6 * ones(1,numel(steps))];
event_xTxt = [event_xTxt, 73.5 * ones(1,numel(steps))];
event_yTxt = [event_yTxt, -16.5 * ones(1,numel(steps))];



% piston top surface heats up: n_step 80 ... 100
steps = 80:100;

event_step_2d = [];
event_text_2d = {};
event_xTipB_2d = [];
event_yTipB_2d = [];
event_xTipE_2d = [];
event_yTipE_2d = [];
event_xTxt_2d = [];
event_yTxt_2d = [];

event_step_2d = [event_step_2d, steps];
event_text_2d = [event_text_2d, repmat({'piston top surface\newline heats up'}, 1, numel(steps))];
event_xTipB_2d = [event_xTipB_2d, 7.0  * ones(1,numel(steps))];
event_yTipB_2d = [event_yTipB_2d, 17.7 * ones(1,numel(steps))];
event_xTipE_2d = [event_xTipE_2d, 6.5  * ones(1,numel(steps))];
event_yTipE_2d = [event_yTipE_2d, 16.5 * ones(1,numel(steps))];
event_xTxt_2d = [event_xTxt_2d, 6.0  * ones(1,numel(steps))];
event_yTxt_2d = [event_yTxt_2d, 16.2 * ones(1,numel(steps))];



% piston bottom 65°C vs water 70°C: n_step 272 ... 300
steps = 272:300;

event_step_2d = [event_step_2d, steps];
event_text_2d = [event_text_2d, repmat({'piston bottom 65°C vs water 70°C'}, 1, numel(steps))];
event_xTipB_2d = [event_xTipB_2d, 4.0 * ones(1,numel(steps))];
event_yTipB_2d = [event_yTipB_2d, 0.3 * ones(1,numel(steps))];
event_xTipE_2d = [event_xTipE_2d, 4.0 * ones(1,numel(steps))];
event_yTipE_2d = [event_yTipE_2d, 1.0 * ones(1,numel(steps))];
event_xTxt_2d = [event_xTxt_2d, 3.4 * ones(1,numel(steps))];
event_yTxt_2d = [event_yTxt_2d, 1.7 * ones(1,numel(steps))];

%%%------------------------------------------%%%
% 05. COMSOL vertical coordinate projected to MATLAB geometry
%%%------------------------------------------%%%

dz_C = 0.05;                                        % COMSOL resolution in m
zC_raw = (0:size(data_2.W,1)-1)' * dz_C;            % 0 ... 38 m from top to bottom

zC_proj = zeros(size(zC_raw));

% Top dead zone: 2.0 m -> 1.665 m
mask_top = zC_raw <= 2.0;
zC_proj(mask_top) = zC_raw(mask_top) * (1.665 / 2.0);

% Middle part: unchanged
mask_mid = zC_raw > 2.0 & zC_raw < 36.0;
zC_proj(mask_mid) = 1.665 + (zC_raw(mask_mid) - 2.0);

% Bottom dead zone: 2.0 m -> 1.665 m
mask_bot = zC_raw >= 36.0;
zC_proj(mask_bot) = 1.665 + 34.0 + (zC_raw(mask_bot) - 36.0) * (1.665 / 2.0);

% Convert to physical z-axis used in plot:
% top = +0.665 m, bottom = -36.665 m
z_C = 0.665 - zC_proj;


%%% COMSOL %%%
% Keep original COMSOL grid and map via interpolation onto projected z-grid
T_ring_comsol = nan(blockLen, n_steps);


%%%------------------------------------------%%%
% 06. COMSOL -> reconstruct global system profile
%%%------------------------------------------%%%

h_piston    = 18;                         % [m]
z_ring_norm = linspace(0,1,blockLen)';    % normalized vertical coordinate

T_sys_comsol = nan(Nz, n_steps);

% Piston lower edge from scenario mapping
z_piston_ts = t_900(2,1:n_steps) * (-19) + (1 - t_900(2,1:n_steps)) * (-35);

% COMSOL piston coordinate (18 m, no projection needed)
zP_rel = (0:size(P_C_raw,1)-1)' * dz_C;     % 0 ... 18 m from piston top to bottom

for i = 1:(n_steps)

    z_piston = z_piston_ts(i);              % lower piston edge
    z_piston_top = z_piston + 18.0;         % upper piston edge

    % --- projected COMSOL water on MATLAB z-grid ---
    T_Wi = interp1(z_C, W_C_raw(:,i), z_M, 'linear', 'extrap');

    % --- piston mean temperature placed into current piston interval ---
    z_Pi = z_piston_top - zP_rel;           % physical piston coordinates, top -> bottom
    T_Pi = interp1(z_Pi, P_C_raw(:,i), z_M, 'linear', 'extrap');

    % start with projected water everywhere
    T_sys_comsol(:,i) = T_Wi;

    % overwrite piston region
    mask_piston = (z_M <= z_piston_top) & (z_M >= z_piston);
    T_sys_comsol(mask_piston,i) = T_Pi(mask_piston);

    % ring-gap temperatures along current ring-gap height
    z_ring = flip(z_piston + z_ring_norm * h_piston);
    T_ring_comsol(:,i) = interp1(z_C, W_C_raw(:,i), z_ring, 'linear', 'extrap');
end

%%
%%%------------------------------------------%%%
% 07. Model container
%%%------------------------------------------%%%

models(1).name   = 'MATLAB baseline';
models(1).T_sys  = T_1;
models(1).T_ring = T_ring_1;
models(1).z      = z_M;
models(1).style  = '-';
models(1).color  = [0 0 0];          % black

models(2).name   = 'COMSOL';
models(2).T_sys  = T_sys_comsol;
models(2).T_ring = T_ring_comsol;
models(2).z      = z_M;
models(2).style  = '-';
models(2).color  = [0.0 0.7 0.0];    % green

models(3).name   = 'MATLAB extended';
models(3).T_sys  = T_3;
models(3).T_ring = T_ring_3;
models(3).T_piston = flip(T_3_p);
models(3).z      = z_M;
models(3).style  = '-';
models(3).color  = [0.6 0.0 0.0];    % blue

z_ring_norm = linspace(0,1,blockLen)';   % normalized vertical coordinate

%%%------------------------------------------%%%
% 08. Video writer
%%%------------------------------------------%%%
dateTag = datestr(now,'yyyymmdd');
version = 'v1';
part = 'part1';
videoname = fullfile(RESULTS_TIMELAPSE, ...
    [dateTag '_' version '_timelapse_1D+2D_temperature_profiles_15min_matlab_comsol_part1.mp4']);

v = VideoWriter(videoname,'MPEG-4');
% Keep it readable; increase if you want a faster video
v.FrameRate = 3;
open(v);

%%%------------------------------------------%%%
% 09. Figure initialization
%%%------------------------------------------%%%

figure('Color','w','Name','1D system + 2D piston temperature profiles – time-lapse (15 min)');

tlo = tiledlayout(1,2,'TileSpacing','loose','Padding','compact');

%%%------------------------------------------%%%
% 10. LEFT: 1D temperature profile
%%%------------------------------------------%%%

ax1 = nexttile;
hold(ax1,'on')
grid(ax1,'on')

% right-side markers (short horizontal segments) for upper/lower bypass and membrane
x_mark = [79.0 80.5];

h_ps_up  = plot(ax1, x_mark, [NaN NaN], 'k-', 'LineWidth',1.5);
h_ps_low = plot(ax1, x_mark, [NaN NaN], 'k-', 'LineWidth',1.5);
h_mem    = plot(ax1, x_mark, [NaN NaN], 'k-', 'LineWidth',1.5);

x_txt = 80.8;   % slightly right of marker lines

h_txt_up  = text(ax1, x_txt, NaN, 'upper bypass', ...
    'FontSize',9, 'Color','k', 'VerticalAlignment','middle');

h_txt_mem = text(ax1, x_txt, NaN, 'membrane', ...
    'FontSize',9, 'Color','k', 'VerticalAlignment','middle');

h_txt_low = text(ax1, x_txt, NaN, 'lower bypass', ...
    'FontSize',9, 'Color','k', 'VerticalAlignment','middle');

for mm = 1:numel(models)

    % Water-Piston system
    h_sys(mm) = plot(ax1, models(mm).T_sys(:,1), models(mm).z, ...
        models(mm).style, ...
        'Color', models(mm).color, ...
        'LineWidth',1.5);

    % Ring-gap handles for model 2 and 3
    if mm >= 2
        h_ring(mm) = plot(ax1, nan, nan, '--', ...
            'Color', models(mm).color, ...
            'LineWidth',1.5);
    end
end

xlabel(ax1,'Temperature [$^{\circ}$C]','Interpreter','Latex')
ylabel(ax1,'System depth [m]','Interpreter','Latex')

% Build legend entries (system + ring-gap)
legend_handles = [];
legend_labels  = {};

for mm = 1:numel(models)

    legend_handles(end+1) = h_sys(mm);
    legend_labels{end+1}  = [models(mm).name ' – system'];

    % handles for ring-gap only for models 2 and 3 (comsol+,matlab)
    if mm >= 2
        legend_handles(end+1) = h_ring(mm);
        legend_labels{end+1}  = [models(mm).name ' – ring gap'];
    end
end

legend(ax1, legend_handles, legend_labels, ...
    'Interpreter','none', 'Location','southeast','HandleVisibility','off')

% set limits
xlim(ax1,[40 81])
ylim(ax1,[min(z_M) max(z_M)])

% Piston rectangle (initialized)
h_piston = 18;   % [m]
z_piston = z_piston_ts(1);

h_piston_top = plot(ax1, [40 81], [z_piston + h_piston, z_piston + h_piston], ...
    ':', 'Color', [0 0 1], 'LineWidth', 1.5);

h_piston_bottom = plot(ax1, [40 81], [z_piston, z_piston], ...
    ':', 'Color', [0 0 1], 'LineWidth', 1.5);


legend(ax1, [legend_handles h_piston_top], [legend_labels {'upper/lower piston end'}], ...
    'Interpreter','none', 'Location','southeast', 'HandleVisibility','off')

% Title handle
h_title = title(ax1,'','Interpreter','Latex');

set(h_ps_up ,  'YData', [z_M(idx_bypass_upper) z_M(idx_bypass_upper)]);
set(h_ps_low,  'YData', [z_M(idx_bypass_lower) z_M(idx_bypass_lower)]);
set(h_txt_up , 'Position', [x_txt z_M(idx_bypass_upper) 0]);
set(h_txt_low, 'Position', [x_txt z_M(idx_bypass_lower) 0]);

% Event annotation handles
max_events_per_frame = 2;

h_event_arrow = gobjects(max_events_per_frame,1);
h_event_text  = gobjects(max_events_per_frame,1);

for kk = 1:max_events_per_frame

    h_event_arrow(kk) = quiver(ax1, NaN, NaN, NaN, NaN, 0, ...
        'Color', [0.1 0.1 0.1], ...
        'LineWidth', 1.5, ...
        'MaxHeadSize', 0.8,'HandleVisibility','off');

    h_event_text(kk) = text(ax1, NaN, NaN, '', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', [0.1 0.1 0.1], ...
        'BackgroundColor', 'None', ...
        'Margin', 3, ...
        'VerticalAlignment', 'middle');
end


%%%------------------------------------------%%%
% 11. RIGHT: 2D radial piston temperature field for Matlab model
%%%------------------------------------------%%%

ax2 = nexttile;

Nz_p = InitOut.grid.Nz(3);
Nr_p = InitOut.grid.Nz(14);

r_vec = InitOut.grid.z_RP(1,:);        % radial node positions
z_vis = linspace(0, h_piston, Nz_p);   % piston height coordinate

% ----------------------------------------------------------
% Uniform radial visualization grid
% ----------------------------------------------------------
dr_vis = 0.005;
r_vis  = 0:dr_vis:max(r_vec);

% ----------------------------------------------------------
% Precompute radial cell assignment indices
% Each r_vis gets temperature of its corresponding ring cell
% ----------------------------------------------------------
idx_r = zeros(size(r_vis));

for k = 1:length(r_vis)
    
    ind = find(r_vis(k) <= r_vec, 1, 'first');
    
    if isempty(ind)
        ind = Nr_p;   % safety for outermost radius
    end
    
    idx_r(k) = ind;
end

% ----------------------------------------------------------
% Initial field
% ----------------------------------------------------------
block   = models(3).T_piston(:,1);
T_field = reshape(block, Nz_p, Nr_p);

% Undo storage flip if necessary
% T_field = rot90(T_field,2);

% Assign ring temperatures piecewise (no interpolation)
T_vis = T_field(:, idx_r);

% ----------------------------------------------------------
% Plot
% ----------------------------------------------------------
h_img = imagesc(ax2, r_vis, z_vis, T_vis);
set(ax2,'YDir','normal')

xlabel(ax2,'Radius [m]')
ylabel(ax2,'Piston height [m]')
title(ax2,'2D piston temperature - MATLAB extended')

colormap(ax2,'turbo')
colorbar(ax2)

clim(ax2,[40 80])
xlim(ax2,[0 max(r_vec)])
axis(ax2,'tight')

% right-side markers (short horizontal segments) for upper/lower bypass and membrane
x_mark_r = [8.0 8.5];

hold(ax2,'on')
h_mem_r    = plot(ax2, x_mark_r, [NaN NaN], 'k-', 'LineWidth',2);


% Event annotation handles for 2D
max_events_per_frame_2d = 2;

h_event_arrow_2d = gobjects(max_events_per_frame_2d,1);
h_event_text_2d  = gobjects(max_events_per_frame_2d,1);

for kk = 1:max_events_per_frame_2d

    h_event_arrow_2d(kk) = quiver(ax2, NaN, NaN, NaN, NaN, 0, ...
        'Color',[0.1 0.1 0.1], ...
        'LineWidth',1.5, ...
        'MaxHeadSize',0.8,'HandleVisibility','off');

    h_event_text_2d(kk) = text(ax2, NaN, NaN, '', ...
        'FontSize',10, ...
        'FontWeight','bold', ...
        'BackgroundColor','none', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top');
end

%%%------------------------------------------%%%
% 12. Animation loop
%%%------------------------------------------%%%

for i = 1:n_steps

    SoC = t_900(2,i);

    % set position of membrane marker based on SoC
    z_membrane = -2*1 - 16 * (1 - SoC/2);
    [~, idx_membrane] = min(abs(z_M - z_membrane));

    set(h_mem,'YData',[z_M(idx_membrane) z_M(idx_membrane)], ...
        'Color',[1 0.85 0]);

    set(h_txt_mem,'Position',[x_txt z_M(idx_membrane) 0]);

    z_piston = z_piston_ts(i);
    z_ring   = flip(z_piston + z_ring_norm * h_piston);

    % set membrane positions for variable piston position
    y_mem_loc = z_M(idx_membrane)     - z_piston;

    % Update 1D system profiles
    for mm = 1:numel(models)
        set(h_sys(mm), ...
            'XData', models(mm).T_sys(:,i), ...
            'YData', models(mm).z);
    end

    % Update ring-gap profiles
    for mm = 2:numel(models)
        set(h_ring(mm), ...
            'XData', models(mm).T_ring(:,i), ...
            'YData', z_ring);
    end

    set(h_piston_top, ...
    'XData', [40 81], ...
    'YData', [z_piston + h_piston, z_piston + h_piston]);

    set(h_piston_bottom, ...
        'XData', [40 81], ...
        'YData', [z_piston, z_piston]);

    % Update 2D piston field
    block = models(3).T_piston(:,i);
    T_field = reshape(block, Nz_p, Nr_p);
    T_vis = T_field(:, idx_r);
    set(h_img,'CData',T_vis);

    h_title.String = sprintf( ...
        '1D system + 2D piston temperature profiles (15 min) - t = %.2f h', ...
        t_hours(i));

    % set membrane dynamically based on current piston position
    set(h_mem_r   , 'YData', [y_mem_loc y_mem_loc], 'Color', [1 0.85 0]);

    % ----------------------------------------------------------
    % Update event annotations
    % ----------------------------------------------------------
    idx_events = find(event_step == i);

    % hide all annotations first
    for kk = 1:max_events_per_frame
        set(h_event_arrow(kk), ...
            'XData', NaN, 'YData', NaN, ...
            'UData', NaN, 'VData', NaN);

        set(h_event_text(kk), ...
            'Position', [NaN NaN 0], ...
            'String', '');
    end

    % show active annotations
    for kk = 1:min(numel(idx_events), max_events_per_frame)

        jj = idx_events(kk);

        % arrow from defined begin point to defined end point
        set(h_event_arrow(kk), ...
            'XData', event_xTipE(jj), ...
            'YData', event_yTipE(jj), ...
            'UData', event_xTipB(jj) - event_xTipE(jj), ...
            'VData', event_yTipB(jj) - event_yTipE(jj));

        set(h_event_text(kk), ...
            'Position', [event_xTxt(jj), event_yTxt(jj), 0], ...
            'String', event_text{jj});
    end


    drawnow
    writeVideo(v, getframe(gcf));


    % ----------------------------------------------------------
    % Update 2D event annotations
    % ----------------------------------------------------------
    idx_events_2d = find(event_step_2d == i);
    
    % hide
    for kk = 1:max_events_per_frame_2d
        set(h_event_arrow_2d(kk), ...
            'XData',NaN,'YData',NaN,'UData',NaN,'VData',NaN);
    
        set(h_event_text_2d(kk), ...
            'Position',[NaN NaN 0], ...
            'String','');
    end
    
    % show
    for kk = 1:min(numel(idx_events_2d), max_events_per_frame_2d)
    
        jj = idx_events_2d(kk);
    
        set(h_event_arrow_2d(kk), ...
        'XData', event_xTipE_2d(jj), ...
        'YData', event_yTipE_2d(jj), ...
        'UData', event_xTipB_2d(jj) - event_xTipE_2d(jj), ...
        'VData', event_yTipB_2d(jj) - event_yTipE_2d(jj));
    
        set(h_event_text_2d(kk), ...
            'Position',[event_xTxt_2d(jj), event_yTxt_2d(jj), 0], ...
            'String', event_text_2d{jj});
    end
end

close(v);

disp(['Saved video: ' videoname]);