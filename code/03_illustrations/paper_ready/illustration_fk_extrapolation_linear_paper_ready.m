%% illustration_fk_extrapolation_academic.m
% Visualization of energy conservation using area arguments
% Extrapolated FK update case
% Units: temperature = 1 unit, height = 1 unit

clear; close all; clc;

%%
run(fullfile('..','..','..', 'configs', 'paths_local.m'));

RESULTS_VIS_BASE  = fullfile(RESULTS_BASE, "03_illustrations");
RESULTS_SNAPSHOTS = fullfile(RESULTS_VIS_BASE, "fk_paper_ready");

%% Figure settings
fig = figure('Color','w', 'Units','centimeters', ...
    'Position',[3 3 14 10]);

tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

colGray = [0.86 0.86 0.86];

zMax = 11;
TMax = 11;

%% =========================================================
% LEFT: Before extrapolated mixing
% =========================================================
ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');

% Temperature profile and extrapolation reference
plot(ax1, [0 8], [0 8], 'k-', 'LineWidth',0.7);
plot(ax1, [0 8], [8 8], 'k-', 'LineWidth',0.7);

plot(ax1, [8 10], [8 10], ...
    'LineWidth',0.7, 'LineStyle','--', 'Color',[0.45 0.45 0.45]);
plot(ax1, [0 10], [10 10], ...
    'LineWidth',0.7, 'LineStyle','--', 'Color',[0.45 0.45 0.45]);

text(ax1, 6.4, 9.45, '\textit{extrapolation}', ...
    'Interpreter','latex', 'Color',[0.35 0.35 0.35], 'FontSize',8);

% Old energy block
fill(ax1, [0 5 8 0], [5 5 8 8], 'w', ...
    'EdgeColor','k', 'FaceAlpha',0.35, 'LineWidth',0.7);
addHatch(ax1, [0 5 8 0], [5 5 8 8], 'vertical', 0.35);

text(ax1, 1.1, 6.65, '$Q_{\mathrm{old}}=19.5$', ...
    'Interpreter','latex', 'FontSize',8);

% Bottom energy
fill(ax1, [0 6 6], [0 7 7], colGray, ...
    'EdgeColor','none', 'FaceAlpha',0.45);
text(ax1, 0.2, 3.7, '$Q_{\mathrm{bottom}}=12$', ...
    'Interpreter','latex', 'FontSize',8);

% Outflow region: dotted only
plot(ax1, [0 0 1 0], [0 1 1 0], 'k:', 'LineWidth',0.7);
text(ax1, 0.5, 0.35, '$Q_{\mathrm{out}}=0.5$', ...
    'Interpreter','latex', 'FontSize',8, 'Color','k');

% Inlet energy regions, dashed boundary
fill(ax1, [0 5 5 0], [6 6 5 5], 'w', ...
    'EdgeColor','k', 'LineStyle','--', 'FaceAlpha',0.35, 'LineWidth',0.5);
addHatch(ax1, [0 5 5 0], [6 6 5 5], 'horizontal', 0.30);

fill(ax1, [5 10 10 5], [6 6 5 5], 'w', ...
    'EdgeColor','k', 'LineStyle','--', 'FaceAlpha',0.35, 'LineWidth',0.5);
addHatch(ax1, [5 10 10 5], [6 6 5 5], 'horizontal', 0.30);

text(ax1, 1.4, 5.35, '$Q_{\mathrm{in},1}=5$', ...
    'Interpreter','latex', 'FontSize',8);
text(ax1, 7.1, 5.35, '$Q_{\mathrm{in},2}=5$', ...
    'Interpreter','latex', 'FontSize',8);

% Annotation
text(ax1, 2.7, 2.0, '$Q_{\mathrm{total}}=19.5+12+0.5+10=42$', ...
    'Interpreter','latex', 'FontSize',8);

formatAxesExtrapolated(ax1, TMax, zMax);

title(ax1, '\textbf{Before mixing}', ...
    'Interpreter','latex', 'FontSize',9.5);

%% =========================================================
% RIGHT: After extrapolated FK update
% =========================================================
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');

% Main updated temperature profile
plot(ax2, [1 5], [0 4], 'k-', 'LineWidth',0.7);
plot(ax2, [5 5], [4 5], 'k-', 'LineWidth',0.7);

% Horizontal inflow part after mixing
addHatch(ax2, [0 5 5 0], [4 4 5 5], 'horizontal', 0.35);
plot(ax2, [0 5 5 0 0], [5 5 4 4 5], 'k--', 'LineWidth',0.7);
text(ax2, 1.6, 4.45, '$Q_{\mathrm{in},1}=5$', ...
    'Interpreter','latex', 'FontSize',8);

% New extrapolated profile reference
plot(ax2, [5 7.38], [5 5], ...
    'LineWidth',0.7, 'LineStyle','--', 'Color',[0.45 0.45 0.45]);
plot(ax2, [7.38 10], [5 10], ...
    'LineWidth',0.7, 'LineStyle','--', 'Color',[0.45 0.45 0.45]);
plot(ax2, [0 10], [10 10], ...
    'LineWidth',0.7, 'LineStyle','--', 'Color',[0.45 0.45 0.45]);

% Old profile inside actual mix zone as dashed reference
plot(ax2, [5 8], [5 8], 'k--', 'LineWidth',0.7);
plot(ax2, [0 8.95], [8 8], 'k--', 'LineWidth',0.7);

% New energy block
fill(ax2, [0 7.38 8.95 0], [5 5 8 8], 'w', ...
    'EdgeColor','k', 'FaceAlpha',0.35, 'LineWidth',0.7);
addHatch(ax2, [0 7.38 8.95 0], [5 5 8 8], 'diagonal', 0.35);

text(ax2, 0.3, 6.55, '$Q_{\mathrm{new}}=19.5+5=24.5$', ...
    'Interpreter','latex', 'FontSize',8);

% Lower contribution
text(ax2, 0.2, 3.4, '$Q_{\mathrm{bottom}}=12$', ...
    'Interpreter','latex', 'FontSize',8);

% Annotation
text(ax2, 3.2, 1.0, '$Q_{\mathrm{total}}=24.5+5+12=41.5$', ...
    'Interpreter','latex', 'FontSize',8);

annotation('arrow', ...
    [0.78 0.78], ...
    [0.55 0.45], ...
    'LineWidth',0.7, ...
    'HeadLength',6, ...
    'HeadWidth',6);
text(ax2, 5.35, 3.2, {'$Q_{\mathrm{in},2}=5$ within', 'actual mix zone'}, ...
    'Interpreter','latex', 'FontSize',8);

formatAxesExtrapolated(ax2, TMax, zMax);

title(ax2, '\textbf{After extrapolated FK update}', ...
    'Interpreter','latex', 'FontSize',9.5);

%% =========================================================
% Custom legend with hatched patches
% =========================================================

axLeg = axes('Position',[0.15 0.02 0.15 0.12]);
axis(axLeg,'off'); hold(axLeg,'on');

y0 = 0.8;
dy = 0.4;
xPatch = [0.1 0.4];
xText  = 0.43;

% --- Vertical hatch ---
drawLegendPatch(axLeg, xPatch, [y0 y0+0.2], 'vertical');
text(axLeg, xText, y0+0.09, ...
    'Thermal energy equivalent of initial temperature profile (mixing zone)', ...
    'Interpreter','latex','FontSize',8, 'VerticalAlignment','middle');

% --- Horizontal hatch ---
drawLegendPatch(axLeg, xPatch, [y0-dy y0-dy+0.2], 'horizontal');
text(axLeg, xText, y0-dy+0.09, ...
    'Cumulative thermal energy inflow during time step', ...
    'Interpreter','latex','FontSize',8, 'VerticalAlignment','middle');

% --- Diagonal hatch ---
drawLegendPatch(axLeg, xPatch, [y0-2*dy y0-2*dy+0.2], 'diagonal');
text(axLeg, xText, y0-2*dy+0.09, ...
    'Thermal energy equivalent after FK update', ...
    'Interpreter','latex','FontSize',8, 'VerticalAlignment','middle');

%% Export
version = 'v1';
export_png_dpi = 900;
dateTag = datestr(now, 'yyyymmdd');
baseName = sprintf('%s_illustration_fk_extrapolation_academic_%s', dateTag, version);

out_pdf = fullfile(RESULTS_SNAPSHOTS, [baseName '.pdf']);

exportgraphics(fig, out_pdf, 'ContentType', 'vector', 'Resolution', export_png_dpi);

disp(['Saved PDF: ' out_pdf]);

%% =========================================================
% Local helper functions
% =========================================================
function formatAxesExtrapolated(ax, TMax, zMax)
    axis(ax, 'equal');
    xlim(ax, [0 TMax]);
    ylim(ax, [0 zMax]);

    xlabel(ax, '$T$', 'Interpreter','latex', 'FontSize',9);
    ylabel(ax, '$z$', 'Interpreter','latex', 'FontSize',9);

    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 8;
    ax.LineWidth = 0.8;

    ax.XTick = [0 5 8 10];
    ax.YTick = [0 5 8 10];

    ax.XTickLabel = { ...
        '$0$', ...
        '$T(z_{\mathrm{in}})$', ...
        '$T(z_{\mathrm{mix}})$', ...
        '$T_{\mathrm{in}}$'};

    ax.YTickLabel = { ...
        '$0$', ...
        '$z_{\mathrm{in}}$', ...
        '$z_{\mathrm{mix}}$', ...
        '$z_{\mathrm{mix}}^*$'};

    % --- Minor ticks ---
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';

    ax.XAxis.MinorTickValues = 0:1:10;
    ax.YAxis.MinorTickValues = 0:1:10;

    % --- Grid ---
    grid(ax, 'on');
    grid(ax, 'minor');

    ax.GridAlpha = 0.1;
    ax.MinorGridAlpha = 0.04;

    ax.MinorGridLineStyle = '-';
    ax.GridLineStyle = '-';
end

function addHatch(ax, xPoly, yPoly, mode, spacing)
    % Add simple hatch lines clipped approximately to polygon area.
    % This helper is sufficient for schematic scientific illustrations.

    xMin = min(xPoly); xMax = max(xPoly);
    yMin = min(yPoly); yMax = max(yPoly);

    switch mode
        case 'vertical'
            xs = xMin:spacing:xMax;
            for x = xs
                [y1, y2] = polygonLineIntersection(xPoly, yPoly, x, 'vertical');
                if ~isempty(y1)
                    plot(ax, [x x], [y1 y2], 'k-', ...
                        'LineWidth',0.25, 'Color',[0 0 0 0.45]);
                end
            end

        case 'horizontal'
            ys = yMin:spacing:yMax;
            for y = ys
                [x1, x2] = polygonLineIntersection(xPoly, yPoly, y, 'horizontal');
                if ~isempty(x1)
                    plot(ax, [x1 x2], [y y], 'k-', ...
                        'LineWidth',0.25, 'Color',[0 0 0 0.45]);
                end
            end

        case 'diagonal'
            cVals = (yMin - xMax):spacing:(yMax - xMin);
            for c = cVals
                [x1, y1, x2, y2] = diagonalPolygonIntersection(xPoly, yPoly, c);
                if ~isempty(x1)
                    plot(ax, [x1 x2], [y1 y2], 'k-', ...
                        'LineWidth',0.25, 'Color',[0 0 0 0.45]);
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
    % Intersections with line y = x + c
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