%% illustration_fk_linear_combination_academic.m
% Visualization of algebraic FK update variants for a limiting case.
% Units: temperature = 1 unit, height = 1 unit

clear; close all; clc;

%%
run(fullfile('..','..','..', 'configs', 'paths_local.m'));

RESULTS_VIS_BASE  = fullfile(RESULTS_BASE, "03_illustrations");
RESULTS_SNAPSHOTS = fullfile(RESULTS_VIS_BASE, "fk_paper_ready");

%% Figure settings
fig = figure('Color','w', 'Units','centimeters', ...
    'Position',[3 3 14 18]);

t = tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Reserve space at bottom for legend
t.OuterPosition = [0 0.1 1 0.9];


colGray = [0.86 0.86 0.86];

zMax = 11;
TMax = 11;

% 2) Default: empty labels
labels = repmat({''}, 1, 10);

% 3) Add labels only at relevant positions
labels{1}  = '1';
labels{2}  = '2';
labels{3}  = '3';
labels{4}  = '4';
labels{5}  = 'T(z_{in})';
labels{6}  = '6';
labels{7}  = '7';
labels{8}  = 'T(z_{mix})';
labels{9}  = '9';
labels{10}  = 'T_{in}';


labelsY = repmat({''}, 1, 10);
labelsY{1}  = '1';
labelsY{2}  = '2';
labelsY{3}  = '3';
labelsY{4}  = '4';
labelsY{5}  = 'z_{in}';
labelsY{6}  = '6';
labelsY{7}  = '7';
labelsY{8}  = 'z_{mix}';
labelsY{9}  = '9';
labelsY{10} = 'z_{mix}^*';

% --- Geometry (from sketch) ---
z_max = 11;
T_max = 11;


%% Geometry and reference quantities
z_in  = 5;
z_mix = 8;

T_old = [5 5 8 8];
z_old = [5 7.9 7.9 8];

Qin2  = 1;
Dz    = z_mix - z_in;
dT_fm = Qin2 / Dz;

%% =========================================================
% UPPER LEFT: Before mixing
% =========================================================
ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');

% Old temperature profile and extrapolation reference
plot(ax1, [0 5], [0 5], 'k-', 'LineWidth',0.7);
plot(ax1, [5 5], [5 7.9], 'k-', 'LineWidth',0.7);
plot(ax1, [5 8], [7.9 7.9], 'k-', 'LineWidth',0.7);
plot(ax1, [0 8], [8 8], 'k-', 'LineWidth',0.7);
plot(ax1, [8 8], [7.9 8], 'k-', 'LineWidth',0.7);

plot(ax1, [5 10], [5 10], '--', ...
    'Color',[0.45 0.45 0.45], 'LineWidth',0.7);
plot(ax1, [0 10], [10 10], '--', ...
    'Color',[0.45 0.45 0.45], 'LineWidth',0.7);

text(ax1, 6.45, 9.45, '\textit{extrapolation}', ...
    'Interpreter','latex', 'Color',[0.35 0.35 0.35], 'FontSize',8);

% Infinitesimal height note
annotation('arrow', ...
    [0.27 0.22], ...
    [0.86 0.89], ...
    'LineWidth',0.7, ...
    'HeadLength',5, ...
    'HeadWidth',5);

text(ax1, 0.7, 9.4, 'infinitesimal height', ...
    'Interpreter','latex', 'FontSize',8);

% Old energy block
fill(ax1, [0 5 5 0], [5 5 8 8], 'w', ...
    'EdgeColor','k', 'LineWidth',0.7);
addHatch(ax1, [0 5 5 0], [5 5 8 8], 'vertical', 0.35);

text(ax1, 1.2, 6.65, '$Q_{\mathrm{old}}=15$', ...
    'Interpreter','latex', 'FontSize',8);

% Bottom energy
fill(ax1, [0 6 6], [0 7 7], colGray, ...
    'EdgeColor','none', 'FaceAlpha',0.45);
text(ax1, 0.2, 3.7, '$Q_{\mathrm{bottom}}\approx 12.5$', ...
    'Interpreter','latex', 'FontSize',8);

% Outflow region: dotted only
plot(ax1, [0 0 0.2 0], [0 0.2 0.2 0], 'k:', 'LineWidth',0.7);
text(ax1, 0.5, 0.35, '$Q_{\mathrm{out}}=1/50\approx 0$', ...
    'Interpreter','latex', 'FontSize',8);

% Inlet energy regions
fill(ax1, [0 5 5 0], [5.2 5.2 5 5], 'w', ...
    'EdgeColor','k', 'LineStyle','--', 'LineWidth',0.5);
addHatch(ax1, [0 5 5 0], [5.2 5.2 5 5], 'horizontal', 0.08);

fill(ax1, [5 10 10 5], [5.2 5.2 5 5], 'w', ...
    'EdgeColor','k', 'LineStyle','--', 'LineWidth',0.5);
addHatch(ax1, [5 10 10 5], [5.2 5.2 5 5], 'horizontal', 0.08);

text(ax1, 1.4, 5.35, '$Q_{\mathrm{in},1}=1$', ...
    'Interpreter','latex', 'FontSize',8);
text(ax1, 7.1, 5.35, '$Q_{\mathrm{in},2}=1$', ...
    'Interpreter','latex', 'FontSize',8);

text(ax1, 3.1, 2.3, '$Q_{\mathrm{total}}\approx 12.5+15+2=29.5$', ...
    'Interpreter','latex', 'FontSize',8);

formatAxesNumeric(ax1, TMax, zMax);

title(ax1, '\textbf{Before mixing}', ...
    'Interpreter','latex', 'FontSize',9.5);

% xticks(1:10); yticks(1:10);
% xticklabels(labels); yticklabels(labelsY);
% xlim([0 T_max]); ylim([0 z_max]);

%% =========================================================
% UPPER RIGHT: Linear approach profile
% =========================================================
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');

% Lower temperature profile
plot(ax2, [0.2 5], [0 4.8], 'k-', 'LineWidth',0.7);
plot(ax2, [5 5], [4.8 5], 'k-', 'LineWidth',0.7);

% New linear approach profile
plot(ax2, [3.333 7.333], [5 8], 'k-', 'LineWidth',0.7);
plot(ax2, [3.333 5], [5 5], 'k-', 'LineWidth',0.7);

% Extrapolation reference
plot(ax2, [7.333 10], [8 10], '--', ...
    'Color',[0.45 0.45 0.45], 'LineWidth',0.7);
plot(ax2, [0 10], [10 10], '--', ...
    'Color',[0.45 0.45 0.45], 'LineWidth',0.7);

text(ax2, 6.45, 9.45, '\textit{extrapolation}', ...
    'Interpreter','latex', 'Color',[0.35 0.35 0.35], 'FontSize',8);

% New energy block
fill(ax2, [0 3.333 7.333 0], [5 5 8 8], 'w', ...
    'EdgeColor','k', 'LineWidth',0.7);
addHatch(ax2, [0 3.333 7.333 0], [5 5 8 8], 'diagonal', 0.35);

text(ax2, 0.15, 7.45, '$Q_{\mathrm{new}}=15+1=16$', ...
    'Interpreter','latex', 'FontSize',8);

% Reference and inflow part
plot(ax2, [5 5], [5 8], 'k--', 'LineWidth',0.7);
plot(ax2, [0 5 5 0 0], [5 5 4.8 4.8 5], 'k--', 'LineWidth',0.7);
addHatch(ax2, [0 5 5 0], [4.8 4.8 5 5], 'horizontal', 0.08);

text(ax2, 1.6, 4.45, '$Q_{\mathrm{in},1}=1$', ...
    'Interpreter','latex', 'FontSize',8);

text(ax2, 0.2, 3.4, '$Q_{\mathrm{bottom}}\approx 12.5$', ...
    'Interpreter','latex', 'FontSize',8);

text(ax2, 3.1, 2.3, '$Q_{\mathrm{total}}\approx 12.5+16+1=29.5$', ...
    'Interpreter','latex', 'FontSize',8);

formatAxesNumeric(ax2, TMax, zMax);

title(ax2, '\textbf{Linear approach profile}', ...
    'Interpreter','latex', 'FontSize',9.5);

%% =========================================================
% LOWER LEFT: Uniform shift profile
% =========================================================
ax3 = nexttile; hold(ax3,'on'); box(ax3,'on');

% Lower temperature profile
plot(ax3, [0.2 5], [0 4.8], 'k-', 'LineWidth',0.7);
plot(ax3, [5 5], [4.8 5], 'k-', 'LineWidth',0.7);

% Uniform shift profile
plot(ax3, [5+dT_fm 5+dT_fm], [5 7.9], 'k-', 'LineWidth',0.7);
plot(ax3, [5+dT_fm 8+dT_fm], [7.9 7.9], 'k-', 'LineWidth',0.7);
plot(ax3, [5 5+dT_fm], [5 5], 'k-', 'LineWidth',0.7);
plot(ax3, [0 8+dT_fm], [8 8], 'k-', 'LineWidth',0.7);
plot(ax3, [8+dT_fm 8+dT_fm], [7.9 8], 'k-', 'LineWidth',0.7);

% Old profile reference
plot(ax3, [5 5], [5 7.9], '--', ...
    'Color',[0.45 0.45 0.45], 'LineWidth',0.7);
plot(ax3, [5 5+dT_fm], [7.9 7.9], '--', ...
    'Color',[0.45 0.45 0.45], 'LineWidth',0.7);

% New energy block
fill(ax3, [0 5+dT_fm 5+dT_fm 8+dT_fm 8+dT_fm 0], ...
          [5 5       7.9     7.9     8       8], 'w', ...
    'EdgeColor','k', 'LineWidth',0.7);
% Lower-left main rectangle
addHatch(ax3, ...
    [0 5+dT_fm 5+dT_fm 0], ...
    [5 5       7.9     7.9], ...
    'diagonal', 0.35);

% Upper-right small strip
addHatch(ax3, ...
    [0 8+dT_fm 8+dT_fm 0], ...
    [7.9 7.9   8       8], ...
    'diagonal', 0.35);

text(ax3, 0.15, 7.45, '$Q_{\mathrm{new}}=15+1=16$', ...
    'Interpreter','latex', 'FontSize',8);

% Inflow part
plot(ax3, [0 5 5 0 0], [5 5 4.8 4.8 5], 'k--', 'LineWidth',0.7);
addHatch(ax3, [0 5 5 0], [4.8 4.8 5 5], 'horizontal', 0.08);

text(ax3, 1.6, 4.45, '$Q_{\mathrm{in},1}=1$', ...
    'Interpreter','latex', 'FontSize',8);

text(ax3, 0.2, 3.4, '$Q_{\mathrm{bottom}}\approx 12.5$', ...
    'Interpreter','latex', 'FontSize',8);

formatAxesNumeric(ax3, TMax, zMax);

title(ax3, '\textbf{Uniform shift profile}', ...
    'Interpreter','latex', 'FontSize',9.5);

%% =========================================================
% LOWER RIGHT: Lambda-weighted profile
% =========================================================
ax4 = nexttile; hold(ax4,'on'); box(ax4,'on');

% FK profile
T_fk = T_old;
mask_mix = (z_old >= z_in);

T_fk(mask_mix) = interp1( ...
    [z_in z_mix], ...
    [3.333 7.333], ...
    z_old(mask_mix));

% Fully mixed profile
T_fm = T_old;
T_fm(mask_mix) = T_old(mask_mix) + dT_fm;

% Lambda from boundary constraints
T_old_in  = 5;
T_old_mix = 8;

Tfk_in  = 3.333;
Tfk_mix = 7.333;

Tfm_in  = 5 + dT_fm;
Tfm_mix = 8 + dT_fm;

lambda_in  = (Tfm_in  - T_old_in ) / (Tfm_in  - Tfk_in);
lambda_mix = (Tfm_mix - T_old_mix) / (Tfm_mix - Tfk_mix);

lambda = min([lambda_in, lambda_mix, 1]);
lambda = max(lambda,0);

% Lambda blend
T_new = T_old;
T_new(mask_mix) = lambda*T_fk(mask_mix) + (1-lambda)*T_fm(mask_mix);

% New energy block
fill(ax4, [0 5 5.6443 8.16664 8.16664 0], ...
          [5 5 7.9    7.9     8       8], 'w', ...
    'EdgeColor','k', 'LineWidth',0.7);

text(ax4, 0.15, 7.45, '$Q_{\mathrm{new}}=15+1=16$', ...
    'Interpreter','latex', 'FontSize',8);

% Profiles
hWeighted = plot(ax4, T_new, z_old, 'k-',  'LineWidth',0.7);
hFk       = plot(ax4, T_fk,  z_old, 'k-.', 'LineWidth',0.7);
hFm       = plot(ax4, T_fm,  z_old, 'k:',  'LineWidth',0.7);

% Additional lambda-result top line
plot(ax4, [0 8.16664], [8 8], 'k-', 'LineWidth',0.7);
plot(ax4, [5 5.333], [5 5], 'k:', 'LineWidth',0.7);

% Inflow part
plot(ax4, [0 5 5 0 0], [5 5 4.8 4.8 5], 'k--', 'LineWidth',0.7);

% Lower-left main block
addHatch(ax4, ...
    [0 5 5.6443 0], ...
    [5 5      7.9    7.9], ...
    'diagonal', 0.35);

% Upper-right small strip
addHatch(ax4, ...
    [0 8.16664 8.16664 0], ...
    [7.9 7.9    8       8], ...
    'diagonal', 0.35);


text(ax4, 1.6, 4.45, '$Q_{\mathrm{in},1}=1$', ...
    'Interpreter','latex', 'FontSize',8);

% Lower temperature profile
plot(ax4, [0.2 5], [0 4.8], 'k-', 'LineWidth',0.7);
plot(ax4, [5 5], [4.8 5], 'k-', 'LineWidth',0.7);

text(ax4, 0.2, 3.4, '$Q_{\mathrm{bottom}}\approx 12.5$', ...
    'Interpreter','latex', 'FontSize',8);

% Boundary conditions
hBc = plot(ax4, [5 5], [0 5], 'r--', 'LineWidth',0.6);
plot(ax4, [8 8], [0 8], 'r--', 'LineWidth',0.6);

% Local subplot legend
lgd = legend(ax4, [hWeighted hFk hFm hBc], ...
    {'Weighted profile', ...
     'Linear approach profile', ...
     'Uniform shift profile', ...
     'Boundary conditions'}, ...
    'Interpreter','latex', ...
    'FontSize',6.5, ...
    'Location','north', ...
    'Box','off');

lgd.Units = 'normalized';
pos = lgd.Position;
pos(2) = pos(2) - 0.02;   % nach oben schieben
lgd.Position = pos;

formatAxesNumeric(ax4, TMax, zMax);

title(ax4, '\textbf{$\lambda$-weighted profile}', ...
    'Interpreter','latex', 'FontSize',9.5);

%% =========================================================
% Global custom legend with hatched patches
% =========================================================

axLeg = axes('Position',[0.15 0.01 0.15 0.12]);
axis(axLeg,'off'); hold(axLeg,'on');

y0 = 0.25;
dy = 0.15;
xPatch = [0.1 0.4];
xText  = 0.43;

drawLegendPatch(axLeg, xPatch, [y0 y0+0.1], 'vertical');
text(axLeg, xText, y0+0.05, ...
    'Thermal energy equivalent of initial temperature profile (mixing zone)', ...
    'Interpreter','latex','FontSize',8, 'VerticalAlignment','middle');

drawLegendPatch(axLeg, xPatch, [y0-dy y0-dy+0.1], 'horizontal');
text(axLeg, xText, y0-dy+0.05, ...
    'Cumulative thermal energy inflow during time step', ...
    'Interpreter','latex','FontSize',8, 'VerticalAlignment','middle');

drawLegendPatch(axLeg, xPatch, [y0-2*dy y0-2*dy+0.1], 'diagonal');
text(axLeg, xText, y0-2*dy+0.05, ...
    'Thermal energy equivalent after FK update', ...
    'Interpreter','latex','FontSize',8, 'VerticalAlignment','middle');

%% Export
version = 'v1';
export_png_dpi = 900;
dateTag = datestr(now, 'yyyymmdd');
baseName = sprintf('%s_illustration_fk_linear_combination_academic_%s', dateTag, version);

out_pdf = fullfile(RESULTS_SNAPSHOTS, [baseName '.pdf']);

exportgraphics(fig, out_pdf, 'ContentType', 'vector', 'Resolution', export_png_dpi);

disp(['Saved PDF: ' out_pdf]);

%% =========================================================
% Local helper functions
% =========================================================
function formatAxesNumeric(ax, TMax, zMax)
    axis(ax, 'equal');
    xlim(ax, [0 TMax]);
    ylim(ax, [0 zMax]);

    xlabel(ax, '$T$', 'Interpreter','latex', 'FontSize',8);
    ylabel(ax, '$z$', 'Interpreter','latex', 'FontSize',8);

    ax.XTick = [0 1 2 3 4 5 6 7 8 9 10 11];
    ax.YTick = [0 1 2 3 4 5 6 7 8 9 10 11];

    ax.XTickLabel = { ...
        '$0$', ...
        '$1$', ...
        '$2$', ...
        '$3$', ...
        '$4$', ...
        '$T(z_{\mathrm{in}})$', ...
        '$6$', ...
        '$7$', ...
        '$T(z_{\mathrm{mix}})$', ...
        '$9$', ...
        '$T_{\mathrm{in}}$',...
        '11'};

    ax.YTickLabel = { ...
        '$0$', ...
        '$1$', ...
        '$2$', ...
        '$3$', ...
        '$4$', ...
        '$z_{\mathrm{in}}$', ...
        '$6$', ...
        '$7$', ...
        '$z_{\mathrm{mix}}$', ...
        '$9$', ...
        '$z_{\mathrm{mix}}^*$',...
        '11'};

    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 8;
    ax.LineWidth = 0.8;

    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';

    ax.XAxis.MinorTickValues = 0:1:10;
    ax.YAxis.MinorTickValues = 0:1:10;

    grid(ax, 'on');
    grid(ax, 'minor');

    ax.GridAlpha = 0.1;
    ax.MinorGridAlpha = 0.04;
    ax.MinorGridLineStyle = '-';
    ax.GridLineStyle = '-';
end

function addHatch(ax, xPoly, yPoly, mode, spacing)
    % Add simple hatch lines clipped approximately to polygon area.

    xMin = min(xPoly); xMax = max(xPoly);
    yMin = min(yPoly); yMax = max(yPoly);

    hatchColor = [0.55 0.55 0.55];
    hatchWidth = 0.25;

    switch mode
        case 'vertical'
            xs = xMin:spacing:xMax;
            for x = xs
                [y1, y2] = polygonLineIntersection(xPoly, yPoly, x, 'vertical');
                if ~isempty(y1)
                    plot(ax, [x x], [y1 y2], '-', ...
                        'LineWidth',hatchWidth, 'Color',hatchColor);
                end
            end

        case 'horizontal'
            ys = yMin:spacing:yMax;
            for y = ys
                [x1, x2] = polygonLineIntersection(xPoly, yPoly, y, 'horizontal');
                if ~isempty(x1)
                    plot(ax, [x1 x2], [y y], '-', ...
                        'LineWidth',hatchWidth, 'Color',hatchColor);
                end
            end

        case 'diagonal'
            cVals = (yMin - xMax):spacing:(yMax - xMin);
            for c = cVals
                [x1, y1, x2, y2] = diagonalPolygonIntersection(xPoly, yPoly, c);
                if ~isempty(x1)
                    plot(ax, [x1 x2], [y1 y2], '-', ...
                        'LineWidth',hatchWidth, 'Color',hatchColor);
                end
            end
    end
end

function [a1, a2] = polygonLineIntersection(xPoly, yPoly, val, direction)
    n = numel(xPoly);
    pts = [];

    for i = 1:n
        j = mod(i,n) + 1;
        x1 = xPoly(i); x2 = xPoly(j);
        y1 = yPoly(i); y2 = yPoly(j);

        if strcmp(direction,'vertical')
            if (val >= min(x1,x2)) && (val <= max(x1,x2)) && (x1 ~= x2)
                t = (val - x1) / (x2 - x1);
                if t >= 0 && t <= 1
                    pts(end+1) = y1 + t*(y2-y1); %#ok<AGROW>
                end
            end
        else
            if (val >= min(y1,y2)) && (val <= max(y1,y2)) && (y1 ~= y2)
                t = (val - y1) / (y2 - y1);
                if t >= 0 && t <= 1
                    pts(end+1) = x1 + t*(x2-x1); %#ok<AGROW>
                end
            end
        end
    end

    pts = unique(round(pts,10));

    if numel(pts) >= 2
        a1 = min(pts);
        a2 = max(pts);
    else
        a1 = [];
        a2 = [];
    end
end

function [xA, yA, xB, yB] = diagonalPolygonIntersection(xPoly, yPoly, c)
    % Intersections with line y = x + c.

    n = numel(xPoly);
    pts = [];

    for i = 1:n
        j = mod(i,n) + 1;
        x1 = xPoly(i); x2 = xPoly(j);
        y1 = yPoly(i); y2 = yPoly(j);

        dx = x2 - x1;
        dy = y2 - y1;
        denom = dy - dx;

        if abs(denom) > 1e-12
            t = (x1 + c - y1) / denom;
            if t >= 0 && t <= 1
                x = x1 + t*dx;
                y = y1 + t*dy;
                pts(end+1,:) = [x y]; %#ok<AGROW>
            end
        end
    end

    if size(pts,1) >= 2
        pts = unique(round(pts,10), 'rows');
        xA = pts(1,1); yA = pts(1,2);
        xB = pts(end,1); yB = pts(end,2);
    else
        xA = []; yA = []; xB = []; yB = [];
    end
end

function drawLegendPatch(ax, xRange, yRange, mode)
    x = [xRange(1) xRange(2) xRange(2) xRange(1)];
    y = [yRange(1) yRange(1) yRange(2) yRange(2)];

    fill(ax, x, y, 'w', ...
        'EdgeColor','k', ...
        'LineWidth',0.8);

    addHatch(ax, x, y, mode, 0.06);
end