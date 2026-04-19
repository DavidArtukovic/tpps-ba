clc
clear
close all
%%%------------------------------------------%%%
%  timelapse_matlab_15min_2d.m
%
%  SUMMARY:
%   Timelapse visualization of matlab 1D temperature profiles (15 min resolution)
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


%%
DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_NOFK  = fullfile(DATA_BASE, 'scenario1');
DATA_SCEN1_FK  = fullfile(DATA_BASE, 'scenario1_freeConv');
DATA_INIT = fullfile(DATA_BASE, 'init');

% Load init and scenario files
load(fullfile(DATA_INIT, '20260324_d18_hp16.0_gap0.5_2D_chunk009_v1.mat'));       % Synthetic data for initial temperature fields
% Scenario/time information
load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));

% Piston information (only for plotting the piston)
load(fullfile(DATA_BASE, 'scenario1/SzenarioComsol.mat'));

% MATLAB results (same as plot_matlab_profiles.m)
data_3 = load(fullfile(DATA_SCEN1_FK, '260419_d56_h112_Res_Matlab_FK_2d_v11.mat'));

%%
%%%------------------------------------------%%%
% 02. Spatial grid
%%%------------------------------------------%%%

dz_M = 0.005;  % MATLAB resolution
Nz   = size(data_3.ResOut.temperature.system, 1) - 400;

z_M_star = flip((0:Nz-1)' * dz_M);
% Physical coordinate: top at +0.665 m, bottom at -36.665 m
z_M = -1*(226 - z_M_star);

% --- geometric bypass positions (relative to z = 0) ---
z_ps_up_geom  = -1 - 112/2;
z_ps_low_geom = -1 - 112;

% --- corresponding indices on z_M grid ---
[~, idx_bypass_upper] = min(abs(z_M - z_ps_up_geom));
[~, idx_bypass_lower] = min(abs(z_M - z_ps_low_geom));


%%%------------------------------------------%%%
% 03. Extract full 15-min resolution
%%%------------------------------------------%%%

T_3 = data_3.ResOut.temperature.system(401:end,:);
T_3_p = data_3.ResOut.temperature.piston; % 2d piston

n_steps = 60;
%%
%-----------------------------------------
% 2.1 Extract Ring-Gap Info
%------------------------------------------

% Part for Ringgap
blockLen = 1573;
startIdx = 22601 * ones(1, n_steps);   % [1 x Nt]

% --- build row index matrix ---
rowOffsets = (0:blockLen-1)';                   % [389 x 1]
rows = startIdx(:)' + rowOffsets;               % [389 x Nt_reduced]
colsMat = repmat(1:n_steps, blockLen, 1);       % [389 x Nt_reduced]

%%% FK %%%
A        = data_3.ResOut.temperature.water(:,1:n_steps);     % [Nz x Nt]
linInd = sub2ind(size(A), rows, colsMat);
T_ring_3 = A(linInd);                                        % [389 x Nt_reduced]

% clear A linInd startIdx rowOffsets rows

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


z_piston_ts = t_900(2,1:n_steps) * (-113) + (1 - t_900(2,1:n_steps)) * (-225);
%%
%%%------------------------------------------%%%
% 06. Model container
%%%------------------------------------------%%%

models(1).name   = 'no FK 1D';
models(1).T_sys  = T_3;
models(1).T_ring = T_ring_3;
models(1).T_piston = flip(T_3_p);
models(1).z      = z_M;
models(1).style  = '-';
models(1).color  = [0 0 0];          % black

z_ring_norm = linspace(0,1,blockLen)';   % normalized vertical coordinate

%%%------------------------------------------%%%
% 07. Video writer
%%%------------------------------------------%%%
dateTag = datestr(now,'yyyymmdd');
videoname = fullfile(RESULTS_TIMELAPSE, ...
    [dateTag '_timelapse_1D+2D_temperature_profiles_15min_matlab_only.mp4']);

v = VideoWriter(videoname,'MPEG-4');
% Keep it readable; increase if you want a faster video
v.FrameRate = 10;
open(v);
%%%------------------------------------------%%%
% 08. Figure initialization
%%%------------------------------------------%%%

figure('Color','w','Name','1D temperature profiles – time-lapse (15 min)');

tlo = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

%% ==========================================================
% LEFT: 1D temperature profile (unchanged visual appearance)
% ==========================================================

ax1 = nexttile;
hold(ax1,'on')
grid(ax1,'on')

% --- right-side markers (short horizontal segments) ---
x_mark = [78 80];

h_ps_up  = plot(ax1, x_mark, [NaN NaN], 'k-', 'LineWidth',3);
h_ps_low = plot(ax1, x_mark, [NaN NaN], 'k-', 'LineWidth',3);
h_mem    = plot(ax1, x_mark, [NaN NaN], 'k-', 'LineWidth',3);

x_txt = 80.3;   % slightly right of marker lines

h_txt_up  = text(ax1, x_txt, NaN, 'upper bypass inlet', ...
    'FontSize',9, 'Color','k', 'VerticalAlignment','middle');

h_txt_mem = text(ax1, x_txt, NaN, 'membrane', ...
    'FontSize',9, 'Color','k', 'VerticalAlignment','middle');

h_txt_low = text(ax1, x_txt, NaN, 'lower bypass inlet', ...
    'FontSize',9, 'Color','k', 'VerticalAlignment','middle');

for mm = 1:numel(models)

    % Water-Piston system
    h_sys(mm) = plot(ax1, models(mm).T_sys(:,1), models(mm).z, ...
        models(mm).style, ...
        'Color', models(mm).color, ...
        'LineWidth',1.5);


    h_ring(mm) = plot(ax1, nan, nan, '--', ...
        'Color', models(mm).color, ...
        'LineWidth',1.5);
end

xlabel(ax1,'Temperature [$^{\circ}$C]','Interpreter','Latex')
ylabel(ax1,'System depth [m]','Interpreter','Latex')

legend(ax1, [h_sys(1), h_ring(1)], ...
    {'large system – system', 'large system – ring gap'}, ...
    'Interpreter','none', 'Location','southwest')

xlim(ax1,[10 80])
ylim(ax1,[min(z_M) max(z_M)])

% Piston rectangle (initialized)
h_piston = 112;   % [m]
z_piston = z_piston_ts(1);

h_rect = rectangle(ax1, ...
    'Position',[10 z_piston 70 h_piston], ...
    'FaceColor',[1 1 1], ...
    'FaceAlpha',0.12, ...
    'EdgeColor',[0 0 0]);

uistack(h_rect,'bottom');

% FK highlight rectangle (whole volume)
% h_fk_rect = rectangle(ax1, ...
%     'Position',[10 0 10 0], ...
%     'EdgeColor','none', ...
%     'LineWidth',3, ...
%     'FaceColor','none');

% uistack(h_fk_rect,'bottom');

% Title handle
h_title = title(ax1,'','Interpreter','Latex');

set(h_ps_up ,  'YData', [z_M(idx_bypass_upper) z_M(idx_bypass_upper)]);
set(h_ps_low,  'YData', [z_M(idx_bypass_lower) z_M(idx_bypass_lower)]);
set(h_txt_up , 'Position', [x_txt z_M(idx_bypass_upper) 0]);
set(h_txt_low,'Position', [x_txt z_M(idx_bypass_lower) 0]);

%% ==========================================================
% RIGHT: 2D radial piston temperature field (no interpolation)
% ==========================================================

ax2 = nexttile;
%%
Nz_p = 22401;
Nr_p = 12;

R = 56/2;          % radius [m]
Nr = 12;           % number of radial cells

a = R / sum((1:Nr).^2);
dr = a * (1:Nr).^2;     % cell widths from outside to inside

sum(dr)                  % should be 28

r_edge = [0, cumsum(flip(dr))];


r_vec = InitOut.grid.z_RP(1,:);        % radial node positions
z_vis = linspace(0, h_piston, Nz_p);   % piston height coordinate
%%
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
block   = models(1).T_piston(:,1);
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

clim(ax2,[10 80])
xlim(ax2,[0 max(r_vec)])
axis(ax2,'tight')
%%%------------------------------------------%%%
% 09. Animation loop
%%%------------------------------------------%%%

for i = 1:n_steps

    SoC = t_900(2,i);

    z_membrane = -2*1 - 112 * (1 - SoC/2);
    [~, idx_membrane] = min(abs(z_M - z_membrane));

    set(h_mem,'YData',[z_M(idx_membrane) z_M(idx_membrane)], ...
        'Color',[1 0.85 0]);

    set(h_txt_mem,'Position',[x_txt z_M(idx_membrane) 0]);

    z_piston = z_piston_ts(i);
    z_ring   = flip(z_piston + z_ring_norm * h_piston);

    % --- FK visualization ---
    fk_code = data_3.ResOut.fk.logging_code(i);

    % if fk_code == 0
    %     set(h_fk_rect,'EdgeColor','none');
    % else
    % 
    %     if abs(fk_code) == 2
    %         fk_color = [0.85 0.1 0.1];
    %     else
    %         fk_color = [0.1 0.3 0.85];
    %     end
    % 
    %     if fk_code > 0
    %         y0 = z_M(idx_bypass_upper);
    %         height_fk = abs(z_M(idx_bypass_upper) - max(z_M));
    %     else
    %         y0 = min(z_M);
    %         height_fk = z_piston - min(z_M);
    %     end
    % 
    %     % set(h_fk_rect, ...
    %     %     'Position',[10 y0 70 height_fk], ...
    %     %     'EdgeColor',fk_color);
    % end

    % Update 1D system profiles
    for mm = 1:numel(models)
        set(h_sys(mm), ...
            'XData', models(mm).T_sys(:,i), ...
            'YData', models(mm).z);
    end

    % Update ring-gap profiles
    for mm = 1:numel(models)
        set(h_ring(mm), ...
            'XData', models(mm).T_ring(:,i), ...
            'YData', z_ring);
    end

    % set(h_rect,'Position',[10 z_piston 70 h_piston]);

    % Update 2D piston field
    block = models(1).T_piston(:,i);
    T_field = reshape(block, Nz_p, Nr_p);
    T_vis = T_field(:, idx_r);
    set(h_img,'CData',T_vis);

    h_title.String = sprintf( ...
        '1D temperature profiles (15 min) - t = %.2f h', ...
        t_hours(i));

    drawnow
    writeVideo(v, getframe(gcf));
end

close(v);

disp(['Saved video: ' videoname]);