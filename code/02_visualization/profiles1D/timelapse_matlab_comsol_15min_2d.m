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
RESULTS_TIMELAPSE = fullfile(RESULTS_VIS_BASE, "timelapse");

DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_NOFK  = fullfile(DATA_BASE, 'scenario1');
DATA_SCEN1_FK  = fullfile(DATA_BASE, 'scenario1_freeConv');
DATA_INIT = fullfile(DATA_BASE, 'init');

% Load init results and scenario files
load(fullfile(DATA_INIT, '20260324_d18_hp16.0_gap0.5_2D_chunk009_v1.mat')); %2D init file
% Scenario/time information
load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));

% Piston information (only for plotting the piston)
load(fullfile(DATA_BASE, 'scenario1/SzenarioComsol.mat'));

% MATLAB and COMSOL results (same as plot_matlab_profiles.m)
data_1 = load(fullfile(DATA_SCEN1_NOFK, 'd18_h18_Res_Matlab_d18_18.mat')); % Matlab 1D no FK
data_2 = load(fullfile(DATA_SCEN1_BASE, 'T1_1818_900.mat')); % COMSOL
data_3 = load(fullfile(DATA_SCEN1_FK, '260324_d18_h18_Res_Matlab_FK_2d_v6.mat')); % Matlab 2D piston with FK


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

% need to ensure all variables have the same number of time steps for the animation loop
n_steps = min([ ...
    size(T_1,2), ...
    size(T_3,2), ...
    size(T_3_p,2), ...
    size(W_C_raw,2), ...
    size(P_C_raw,2), ...
    size(t_900,2)]);

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

for i = 1:n_steps

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

models(1).name   = 'no FK 1D';
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

models(3).name   = 'FK 2D piston';
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
videoname = fullfile(RESULTS_TIMELAPSE, ...
    [dateTag '_timelapse_1D+2D_temperature_profiles_15min_matlab_comsol.mp4']);

v = VideoWriter(videoname,'MPEG-4');
% Keep it readable; increase if you want a faster video
v.FrameRate = 8;
open(v);

%%%------------------------------------------%%%
% 09. Figure initialization
%%%------------------------------------------%%%

figure('Color','w','Name','1D system + 2D piston temperature profiles – time-lapse (15 min)');

tlo = tiledlayout(1,2,'TileSpacing','loose','Padding','compact');

%% ==========================================================
% LEFT: 1D temperature profile (unchanged visual appearance)
% ==========================================================

ax1 = nexttile;
hold(ax1,'on')
grid(ax1,'on')

% right-side markers (short horizontal segments) for upper/lower bypass and membrane
x_mark = [78.8 80.3];

h_ps_up  = plot(ax1, x_mark, [NaN NaN], 'k-', 'LineWidth',2);
h_ps_low = plot(ax1, x_mark, [NaN NaN], 'k-', 'LineWidth',2);
h_mem    = plot(ax1, x_mark, [NaN NaN], 'k-', 'LineWidth',2);

x_txt = 80.3;   % slightly right of marker lines

h_txt_up  = text(ax1, x_txt, NaN, 'upper inlet', ...
    'FontSize',9, 'Color','k', 'VerticalAlignment','middle');

h_txt_mem = text(ax1, x_txt, NaN, 'membrane', ...
    'FontSize',9, 'Color','k', 'VerticalAlignment','middle');

h_txt_low = text(ax1, x_txt, NaN, 'lower inlet', ...
    'FontSize',9, 'Color','k', 'VerticalAlignment','middle');

for mm = 1:numel(models)

    % Water-Piston system
    h_sys(mm) = plot(ax1, models(mm).T_sys(:,1), models(mm).z, ...
        models(mm).style, ...
        'Color', models(mm).color, ...
        'LineWidth',1.5);

    % --- Ring-gap handles for model 2 and 3 ---
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

    if mm >= 2
        legend_handles(end+1) = h_ring(mm);
        legend_labels{end+1}  = [models(mm).name ' – ring gap'];
    end
end

legend(ax1, legend_handles, legend_labels, ...
    'Interpreter','none', 'Location','southwest')

xlim(ax1,[40 80])
ylim(ax1,[min(z_M) max(z_M)])

% Piston rectangle (initialized)
h_piston = 18;   % [m]
z_piston = z_piston_ts(1);

h_rect = rectangle(ax1, ...
    'Position',[40 z_piston 40 h_piston], ...
    'FaceColor',[1 1 1], ...
    'FaceAlpha',0.12, ...
    'EdgeColor',[0 0 0]);

uistack(h_rect,'bottom');

% FK highlight rectangle (whole volume)
h_fk_rect = rectangle(ax1, ...
    'Position',[40 0 40 0], ...
    'EdgeColor','none', ...
    'LineWidth',3, ...
    'FaceColor','none');

uistack(h_fk_rect,'bottom');

% Title handle
h_title = title(ax1,'','Interpreter','Latex');

set(h_ps_up ,  'YData', [z_M(idx_bypass_upper) z_M(idx_bypass_upper)]);
set(h_ps_low,  'YData', [z_M(idx_bypass_lower) z_M(idx_bypass_lower)]);
set(h_txt_up , 'Position', [x_txt z_M(idx_bypass_upper) 0]);
set(h_txt_low, 'Position', [x_txt z_M(idx_bypass_lower) 0]);

%% ==========================================================
% RIGHT: 2D radial piston temperature field (no interpolation)
% ==========================================================

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
title(ax2,'2D piston temperature')

colormap(ax2,'turbo')
colorbar(ax2)

clim(ax2,[40 80])
xlim(ax2,[0 max(r_vec)])
axis(ax2,'tight')

% right-side markers (short horizontal segments) for upper/lower bypass and membrane
x_mark_r = [8.2 8.5];

hold(ax2,'on')
h_ps_up_r  = plot(ax2, x_mark_r, [NaN NaN], 'k-', 'LineWidth',2);
h_ps_low_r = plot(ax2, x_mark_r, [NaN NaN], 'k-', 'LineWidth',2);
h_mem_r    = plot(ax2, x_mark_r, [NaN NaN], 'k-', 'LineWidth',2);

%%%------------------------------------------%%%
% 10. Animation loop
%%%------------------------------------------%%%

for i = 1:n_steps

    SoC = t_900(2,i);

    z_membrane = -2*1 - 16 * (1 - SoC/2);
    [~, idx_membrane] = min(abs(z_M - z_membrane));

    set(h_mem,'YData',[z_M(idx_membrane) z_M(idx_membrane)], ...
        'Color',[1 0.85 0]);

    set(h_txt_mem,'Position',[x_txt z_M(idx_membrane) 0]);

    z_piston = z_piston_ts(i);
    z_ring   = flip(z_piston + z_ring_norm * h_piston);

    % set inlets and membrane positions for variable piston position
    y_up_loc  = z_M(idx_bypass_upper) - z_piston;
    y_low_loc = z_M(idx_bypass_lower) - z_piston+0.05;
    y_mem_loc = z_M(idx_membrane)     - z_piston;


    % --- FK visualization ---
    fk_code = data_3.ResOut.fk.logging_code(i);

    if fk_code == 0
        set(h_fk_rect,'EdgeColor','none');
    else

        if abs(fk_code) == 2
            fk_color = [0.85 0.1 0.1];
        else
            fk_color = [0.1 0.3 0.85];
        end

        if fk_code > 0
            y0 = z_M(idx_bypass_upper);
            height_fk = abs(z_M(idx_bypass_upper) - max(z_M));
        else
            y0 = min(z_M);
            height_fk = z_piston - min(z_M);
        end

        set(h_fk_rect, ...
            'Position',[40 y0 40 height_fk], ...
            'EdgeColor',fk_color);
    end

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

    set(h_rect,'Position',[40 z_piston 40 h_piston]);

    % Update 2D piston field
    block = models(3).T_piston(:,i);
    T_field = reshape(block, Nz_p, Nr_p);
    T_vis = T_field(:, idx_r);
    set(h_img,'CData',T_vis);

    h_title.String = sprintf( ...
        '1D system + 2D piston temperature profiles (15 min) - t = %.2f h', ...
        t_hours(i));

    
    set(h_ps_up_r , 'YData', [y_up_loc y_up_loc]);
    set(h_ps_low_r, 'YData', [y_low_loc y_low_loc]);
    set(h_mem_r   , 'YData', [y_mem_loc y_mem_loc], 'Color', [1 0.85 0]);

    drawnow
    % writeVideo(v, getframe(gcf));
end

close(v);

disp(['Saved video: ' videoname]);