clc
clear
close all
%%%------------------------------------------%%%
%  plot_snapshots_matlab_comsol_2d_piston.m
%
%  SUMMARY:
%   Generate six academic-style snapshot plots for the comparison of
%   MATLAB / COMSOL 1D system temperature profiles and the 2D piston
%   temperature field from the MATLAB FK-2D model.
%
%  NOTES:
%   - Each snapshot is shown as a pair:
%       left  = 1D system temperature profile
%       right = 2D piston temperature field
%   - The six snapshots are arranged in a compact 3x4 tiled layout.
%   - The figure is designed to remain readable in a thesis context.
%%%------------------------------------------%%%

%%%------------------------------------------%%%
% 01. User settings
%%%------------------------------------------%%%

snapshot_times_h = [0.25, 10, 59, 149, 275, 450];

c_base = [0.15 0.15 0.15];   % dark gray
c_cmsl = [0.00 0.45 0.70];   % blue
c_ext  = [0.80 0.40 0.00];   % orange

c_mem = [0 0 0];   % blue
c_inlet  = [0.49 0.18 0.56];   % purple

lw_main     = 0.65;
lw_inlet    = 0.7;
lw_mark     = 0.55;

fs_axis   = 7;
fs_label  = 5.5;
fs_title  = 7;
fs_legend = 10;

T_lim = [40 80];
z_lim = [-36.665 0.665];

export_png_dpi = 1200;

%%%------------------------------------------%%%
% 02. Load paths and data
%%%------------------------------------------%%%

run(fullfile('..','..','..', 'configs', 'paths_local.m'));

RESULTS_VIS_BASE  = fullfile(RESULTS_BASE, "02_visualizations");
RESULTS_SNAPSHOTS = fullfile(RESULTS_VIS_BASE, "03_snapshots");

if ~exist(RESULTS_SNAPSHOTS, 'dir')
    mkdir(RESULTS_SNAPSHOTS);
end

DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_NOFK = fullfile(DATA_BASE, 'scenario1');
DATA_SCEN1_FK   = fullfile(DATA_BASE, 'scenario1_freeConv');
DATA_INIT       = fullfile(DATA_BASE, 'init');

load(fullfile(DATA_INIT, '20260330_d18_hp18.0_gap0.5_2D_chunk009.mat'));
load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));
load(fullfile(DATA_BASE, 'scenario1/SzenarioComsol.mat'));

data_1 = load(fullfile(DATA_SCEN1_NOFK, 'd18_h18_Res_Matlab_d18_18.mat'));
data_2 = load(fullfile(DATA_SCEN1_BASE, 'T1_1818_900.mat'));
data_3 = load(fullfile(DATA_SCEN1_FK, '260331_d18_h18_Res_Matlab_FK_2d_v12.mat'));

%%%------------------------------------------%%%
% 03. Spatial grids and geometry
%%%------------------------------------------%%%

dz_M = 0.005;
Nz   = size(data_1.Res_System_d18_18,1) - 400;

z_M_star = flip((0:Nz-1)' * dz_M);
z_M      = -1 * (37.33 - z_M_star) + 0.665;

z_ps_up_geom  = -1 - 16/2;
z_ps_low_geom = -3*1 - 0.05 - 16;

[~, idx_bypass_upper] = min(abs(z_M - z_ps_up_geom));
[~, idx_bypass_lower] = min(abs(z_M - z_ps_low_geom));

%%%------------------------------------------%%%
% 04. Extract full 15-min resolution
%%%------------------------------------------%%%

T_1   = data_1.Res_System_d18_18(401:end,:);
T_3   = data_3.ResOut.temperature.system(401:end,:);
T_3_p = data_3.ResOut.temperature.piston;

W_C_raw = data_2.W(:,2:end);
P_C_raw = data_2.P(:,2:end);

n_steps = min([ ...
    size(T_1,2), ...
    size(T_3,2), ...
    size(T_3_p,2), ...
    size(W_C_raw,2), ...
    size(P_C_raw,2), ...
    size(t_900,2)]);

T_1     = T_1(:,1:n_steps);
T_3     = T_3(:,1:n_steps);
T_3_p   = T_3_p(:,1:n_steps);
W_C_raw = W_C_raw(:,1:n_steps);
P_C_raw = P_C_raw(:,1:n_steps);

%%%------------------------------------------%%%
% 05. Extract ring-gap temperatures for MATLAB models
%%%------------------------------------------%%%

blockLen = 389;
startIdx = data_1.Res_900_d18_18(8,2:n_steps+1) - 1;

rowOffsets = (0:blockLen-1)';
rows       = startIdx(:)' + rowOffsets;
colsMat    = repmat(1:n_steps, blockLen, 1);

A = data_1.Res_Wasser_d18_18(:,1:n_steps);
linInd = sub2ind(size(A), rows, colsMat);
T_ring_1 = A(linInd);

A = data_3.ResOut.temperature.water(:,1:n_steps);
linInd = sub2ind(size(A), rows, colsMat);
T_ring_3 = A(linInd);

clear A linInd rows colsMat rowOffsets startIdx

%%%------------------------------------------%%%
% 06. Time axis
%%%------------------------------------------%%%

dt_min  = 15;
t_hours = (1:n_steps) * dt_min / 60;

%%%------------------------------------------%%%
% 07. COMSOL vertical coordinate projected to MATLAB geometry
%%%------------------------------------------%%%

dz_C  = 0.05;
zC_raw = (0:size(data_2.W,1)-1)' * dz_C;

zC_proj = zeros(size(zC_raw));

mask_top = zC_raw <= 2.0;
zC_proj(mask_top) = zC_raw(mask_top) * (1.665 / 2.0);

mask_mid = zC_raw > 2.0 & zC_raw < 36.0;
zC_proj(mask_mid) = 1.665 + (zC_raw(mask_mid) - 2.0);

mask_bot = zC_raw >= 36.0;
zC_proj(mask_bot) = 1.665 + 34.0 + (zC_raw(mask_bot) - 36.0) * (1.665 / 2.0);

z_C = 0.665 - zC_proj;

%%%------------------------------------------%%%
% 08. Reconstruct COMSOL global system profile
%%%------------------------------------------%%%

h_piston    = 18.0;
z_ring_norm = linspace(0,1,blockLen)';

T_ring_comsol = nan(blockLen, n_steps);
T_sys_comsol  = nan(Nz, n_steps);

z_piston_ts = t_900(2,1:n_steps) * (-19) + (1 - t_900(2,1:n_steps)) * (-35);
zP_rel      = (0:size(P_C_raw,1)-1)' * dz_C;

for i = 1:n_steps

    z_piston_i     = z_piston_ts(i);
    z_piston_top_i = z_piston_i + h_piston;

    T_Wi = interp1(z_C, W_C_raw(:,i), z_M, 'linear', 'extrap');

    z_Pi = z_piston_top_i - zP_rel;
    T_Pi = interp1(z_Pi, P_C_raw(:,i), z_M, 'linear', 'extrap');

    T_sys_comsol(:,i) = T_Wi;

    mask_piston = (z_M <= z_piston_top_i) & (z_M >= z_piston_i);
    T_sys_comsol(mask_piston,i) = T_Pi(mask_piston);

    z_ring_i = flip(z_piston_i + z_ring_norm * h_piston);
    T_ring_comsol(:,i) = interp1(z_C, W_C_raw(:,i), z_ring_i, 'linear', 'extrap');
end

%%%------------------------------------------%%%
% 09. Model container
%%%------------------------------------------%%%

models(1).name    = '1D MATLAB';
models(1).T_sys   = T_1;
models(1).T_ring  = T_ring_1;
models(1).z       = z_M;
models(1).style   = '-';
models(1).color   = c_base;

models(2).name    = 'COMSOL';
models(2).T_sys   = T_sys_comsol;
models(2).T_ring  = T_ring_comsol;
models(2).z       = z_M;
models(2).style   = '-';
models(2).color   = c_cmsl;

models(3).name    = '2D MATLAB + FK';
models(3).T_sys   = T_3;
models(3).T_ring  = T_ring_3;
models(3).T_piston = flip(T_3_p);
models(3).z       = z_M;
models(3).style   = '-';
models(3).color   = c_ext;

%%%------------------------------------------%%%
% 10. 2D piston visualization grid
%%%------------------------------------------%%%

Nz_p = InitOut.grid.Nz(3);
Nr_p = InitOut.grid.Nz(14);

r_vec = InitOut.grid.z_RP(1,:);
z_vis = linspace(0, h_piston, Nz_p);

dr_vis = 0.005;
r_vis  = 0:dr_vis:max(r_vec);

idx_r = zeros(size(r_vis));

for k = 1:length(r_vis)

    ind = find(r_vis(k) <= r_vec, 1, 'first');

    if isempty(ind)
        ind = Nr_p;
    end

    idx_r(k) = ind;
end

%%%------------------------------------------%%%
% 11. Snapshot indices
%%%------------------------------------------%%%

snap_idx = zeros(size(snapshot_times_h));

for k = 1:numel(snapshot_times_h)
    [~, snap_idx(k)] = min(abs(t_hours - snapshot_times_h(k)));
end

snapshot_times_h = t_hours(snap_idx);

%%%------------------------------------------%%%
% 12. Figure setup
%%%------------------------------------------%%%

width_cm  = 17.0;
height_cm = 17.0;

fig = figure( ...
    'Color', 'w', ...
    'Name', 'MATLAB / COMSOL snapshot comparison', ...
    'Units', 'centimeters', ...
    'Position', [2 2 width_cm height_cm]);

set(fig, 'Renderer', 'painters');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [width_cm height_cm]);
set(fig, 'PaperPosition', [0 0 width_cm height_cm]);
set(fig, 'DefaultAxesFontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman');

tlo = tiledlayout(fig, 3, 4, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

colormap(fig, turbo);

legend_handles = gobjects(0);
legend_labels  = {};

panel_labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

%%%------------------------------------------%%%
% 13. Plot snapshots
%%%------------------------------------------%%%

for s = 1:numel(snap_idx)

    i = snap_idx(s);

    z_piston_i     = z_piston_ts(i);
    z_piston_top_i = z_piston_i + h_piston;
    z_ring_i       = flip(z_piston_i + z_ring_norm * h_piston);

    SoC = t_900(2,i);
    z_membrane = -2*1 - 16 * (1 - SoC/2);
    [~, idx_membrane] = min(abs(z_M - z_membrane));

    y_mem_local = z_M(idx_membrane) - z_piston_i;

    %%%--------------------------------------%%%
    % LEFT PANEL: 1D system profile
    %%%--------------------------------------%%%

    axL = nexttile(tlo, 2*(s-1) + 1);
    hold(axL, 'on')
    box(axL, 'on')
    grid(axL, 'on')

    % Model system lines
    h1 = plot(axL, models(1).T_sys(:,i), models(1).z, ...
        '-', 'Color', models(1).color, 'LineWidth', lw_main);

    h2 = plot(axL, models(2).T_sys(:,i), models(2).z, ...
        '-', 'Color', models(2).color, 'LineWidth', lw_main);

    h3 = plot(axL, models(3).T_sys(:,i), models(3).z, ...
        '-', 'Color', models(3).color, 'LineWidth', lw_main);

    % Ring-gap lines
    h4 = plot(axL, models(2).T_ring(:,i), z_ring_i, ...
        ':', 'Color', models(2).color, 'LineWidth', lw_main);

    h5 = plot(axL, models(3).T_ring(:,i), z_ring_i, ...
        ':', 'Color', models(3).color, 'LineWidth', lw_main);

    % Piston end lines
    h6 = plot(axL, T_lim, [z_piston_top_i z_piston_top_i], ...
        '--', 'Color', [0.25 0.25 0.25], 'LineWidth', lw_mark);

    plot(axL, T_lim, [z_piston_i z_piston_i], ...
        '--', 'Color', [0.25 0.25 0.25], 'LineWidth', lw_mark);

    % Interface markers (short segments at right edge)
    x_mark = [78.0 80.0];

    plot(axL, x_mark, [z_M(idx_bypass_upper) z_M(idx_bypass_upper)], ...
        '-', 'Color', c_inlet, 'LineWidth', lw_inlet);

    plot(axL, x_mark, [z_M(idx_membrane) z_M(idx_membrane)], ...
        '-', 'Color', c_mem, 'LineWidth', lw_inlet);

    plot(axL, x_mark, [z_M(idx_bypass_lower) z_M(idx_bypass_lower)], ...
        '-', 'Color', c_inlet, 'LineWidth', lw_inlet);

    % set y label only left
    if mod(s,2) == 0
        x_txt = 80.2;
        text(axL, x_txt, z_M(idx_bypass_upper), 'upper inlet', ...
            'FontSize', fs_label, ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left', ...
            'Clipping', 'off');
        
        text(axL, x_txt, z_M(idx_membrane), 'membrane', ...
            'FontSize', fs_label, ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left', ...
            'Clipping', 'off');
        
        text(axL, x_txt, z_M(idx_bypass_lower), 'lower inlet', ...
            'FontSize', fs_label, ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left', ...
            'Clipping', 'off');
    end


    xlim(axL, T_lim);
    ylim(axL, z_lim);
    axL.XTick = [40 50 60 70 80];

    axL.FontSize = fs_axis;
    axL.GridAlpha = 0.15;
    axL.LineWidth = 0.25;
    axL.Layer = 'top';

    % set y label only left
    if mod(s,2) == 1
        ylabel(axL, 'System depth [m]')
    else
        ylabel(axL, '')
    end

    % set x label only at bottom
    if s >= 5
        xlabel(axL, 'Temperature [°C]')
    else
        xlabel(axL, '')
        axL.XTickLabel = [];
    end


    time_str = format_snapshot_time(snapshot_times_h(s));

    title(axL, sprintf('%s  %s', panel_labels{s}, time_str), ...
        'FontWeight', 'normal', 'FontSize', fs_title)

    if s == 1
        legend_handles = [h1 h3 h2 h5 h4 h6];
        legend_labels  = { ...
            'MATLAB baseline – system', ...
            'MATLAB extended – system', ...
            'COMSOL – system', ...
            'MATLAB extended – ring gap', ...
            'COMSOL – ring gap', ...
            'upper/lower piston end'}; 
    
        h_box = patch(NaN, NaN, [0.8 0.8 0.8], ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', [0.5 0.5 0.5]);

        legend_handles = [legend_handles h_box];
        legend_labels  = [legend_labels {'highlighted region'}];
    end
%--------------------------------------%
% Highlight regions for discussion
%--------------------------------------%

x_rect = 41;
w_rect = 34;   % width across most of temperature axis

switch s

    % Plot 2 → s = 2 → 0.665 to -7 m
    case 2
        y_top = 0.665;
        y_bot = -7;
        x_rect = 43;
        w_rect = 38;   
        rectangle(axL, ...
        'Position', [x_rect, y_bot, w_rect, y_top - y_bot], ...
        'FaceColor', [0.8 0.8 0.8], ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', [0.5 0.5 0.5], ...
        'LineStyle', '-', ...
        'LineWidth', 0.3);

    % Plot 3 → -10 to -20 m
    case 3
        y_top = -10;
        y_bot = -19.6;
        x_rect = 59;
        w_rect = 11; 
        rectangle(axL, ...
        'Position', [x_rect, y_bot, w_rect, y_top - y_bot], ...
        'FaceColor', [0.8 0.8 0.8], ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', [0.5 0.5 0.5], ...
        'LineStyle', '-', ...
        'LineWidth', 0.3);


    % Plot 4 → -18.5 to -19.5 m
    case 4
        y_top = -17.5;
        y_bot = -22.5;
        x_rect = 66;
        w_rect = 10; 
        rectangle(axL, ...
        'Position', [x_rect, y_bot, w_rect, y_top - y_bot], ...
        'FaceColor', [0.8 0.8 0.8], ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', [0.5 0.5 0.5], ...
        'LineStyle', '-', ...
        'LineWidth', 0.3);

    % Plot 6 → two regions
    case 6
        % region 1: -8.5 to -10
        y_top = -8.5;
        y_bot = -9.8;
        x_rect = 53;
        w_rect = 10; 
        rectangle(axL, ...
        'Position', [x_rect, y_bot, w_rect, y_top - y_bot], ...
        'FaceColor', [0.8 0.8 0.8], ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', [0.5 0.5 0.5], ...
        'LineStyle', '-', ...
        'LineWidth', 0.3);

        % region 2: -18.7 to -19.3
        y_top = -17;
        y_bot = -19.8;
        x_rect = 51;
        w_rect = 19; 
        rectangle(axL, ...
        'Position', [x_rect, y_bot, w_rect, y_top - y_bot], ...
        'FaceColor', [0.8 0.8 0.8], ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', [0.5 0.5 0.5], ...
        'LineStyle', '-', ...
        'LineWidth', 0.3);
end
    %%%--------------------------------------%%%
    % RIGHT PANEL: 2D piston field
    %%%--------------------------------------%%%

    axR = nexttile(tlo, 2*(s-1) + 2);
    hold(axR, 'on')
    box(axR, 'on')

    block   = models(3).T_piston(:,i);
    T_field = reshape(block, Nz_p, Nr_p);
    T_vis   = T_field(:, idx_r);

    imagesc(axR, r_vis, z_vis, T_vis)
    set(axR, 'YDir', 'normal')


    % Membrane marker in piston-local coordinates
    plot(axR, [max(r_vec)-0.5 max(r_vec)], [y_mem_local y_mem_local], ...
        '-', 'Color', [0 0 0], 'LineWidth', lw_mark)

    xlim(axR, [0 max(r_vec)])
    ylim(axR, [0 h_piston])
    clim(axR, T_lim)

    axR.FontSize = fs_axis;
    axR.LineWidth = 0.25;
    axR.Layer = 'top';
    grid(axR, 'off');


    % Y label only for left column of piston plots
    if mod(s,2) == 1
        ylabel(axR, 'Piston height [m]')
    else
        ylabel(axR, '')
        axR.YTickLabel = [];
    end


    % X label only bottom row
    if s >= 5
        xlabel(axR, 'Radius [m]')
    else
        xlabel(axR, '')
        axR.XTickLabel = [];
    end

    axL.Box = 'on';
    axR.Box = 'on';
    axL.Layer = 'top';
    axR.Layer = 'top';


    % set y label only left
    if s==2 || s==4
        x_txt = 8.6;

        text(axR, x_txt, 9, 'membrane', ...
            'FontSize', fs_axis, ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left', ...
            'Clipping', 'off');
    end
end

%%%------------------------------------------%%%
% 14. Shared colorbar and legend
%%%------------------------------------------%%%

cb = colorbar;
cb.Location = 'eastoutside';
cb.Label.String = 'Temperature [°C]';
cb.FontSize = fs_axis;
cb.Limits = T_lim;

lgd = legend(legend_handles, legend_labels, ...
    'NumColumns', 4, ...
    'Location', 'southoutside', ...
    'Interpreter', 'none', ...
    'Box', 'off');

lgd.Layout.Tile = 'south';
lgd.FontSize = fs_axis;

%%%------------------------------------------%%%
% 15. Export
%%%------------------------------------------%%%
version = 'v6';
dateTag = datestr(now, 'yyyymmdd');
baseName = sprintf('%s_snapshots_matlab_comsol_2d_piston_%s', dateTag, version);

out_png = fullfile(RESULTS_SNAPSHOTS, [baseName '.png']);
out_pdf = fullfile(RESULTS_SNAPSHOTS, [baseName '.pdf']);

exportgraphics(fig, out_png, 'Resolution', export_png_dpi);
exportgraphics(fig, out_pdf, 'ContentType', 'vector', 'Resolution', export_png_dpi);

disp(['Saved PNG: ' out_png]);
disp(['Saved PDF: ' out_pdf]);
%% Local Function

function time_str = format_snapshot_time(t_h)
% Format snapshot time for panel titles
if t_h <1
    time_str = sprintf('t = %.2f h', 0);
elseif t_h < 24
    time_str = sprintf('t = %.2f h', t_h);
else
    n_days = floor(t_h / 24);
    n_hours = round(t_h - 24 * n_days);

    if n_hours == 24
        n_days = n_days + 1;
        n_hours = 0;
    end

    if n_hours == 0
        time_str = sprintf('t = %d days', n_days);
    else
        time_str = sprintf('t = %d days, %dh', n_days, n_hours);
    end
end
end