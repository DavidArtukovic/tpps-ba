%% illustration_fk_replacement_academic.m
% Visualization of energy conservation using area arguments
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

% Colors
colGray   = [0.86 0.86 0.86];

zMax = 11;
TMax = 11;

%% =========================================================
% LEFT: Before mixing
% =========================================================
ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');

% Main temperature profile
plot(ax1, [0 10], [0 10], 'k-', 'LineWidth', 0.7);

% Old energy block
fill(ax1, [0 6 10 0], [6 6 10 10], 'w', ...
    'EdgeColor','k', 'FaceAlpha',0.35, 'LineWidth',0.7);
addHatch(ax1, [0 6 10 0], [6 6 10 10], 'vertical', 0.35);

text(ax1, 2.0, 8, '$Q_{\mathrm{old}}=32$', ...
    'Interpreter','latex', 'FontSize',8);

% Bottom energy
fill(ax1, [0 6 6], [0 7 7], colGray, ...
    'EdgeColor','none', 'FaceAlpha',0.45);
text(ax1, 0.2, 4.2, '$Q_{\mathrm{bottom}}=17.5$', ...
    'Interpreter','latex', 'FontSize',8);

% Outflow region: dotted only, no orange fill
plot(ax1, [0 0 1 0], [0 1 1 0], 'k:', 'LineWidth',0.7);
text(ax1, 0.5, 0.35, '$Q_{\mathrm{out}}=0.5$', ...
    'Interpreter','latex', 'FontSize',8, 'Color','k');

% Inlet energy regions, dashed boundary
fill(ax1, [0 6 6 0], [7 7 6 6], 'w', ...
    'EdgeColor','k', 'LineStyle','--', 'FaceAlpha',0.35, 'LineWidth',0.5);
addHatch(ax1, [0 6 6 0], [7 7 6 6], 'horizontal', 0.30);

fill(ax1, [6 10 10 6], [7 7 6 6], 'w', ...
    'EdgeColor','k', 'LineStyle','--', 'FaceAlpha',0.35, 'LineWidth',0.5);
addHatch(ax1, [6 10 10 6], [7 7 6 6], 'horizontal', 0.30);

text(ax1, 2.4, 6.35, '$Q_{\mathrm{in},1}=6$', ...
    'Interpreter','latex', 'FontSize',8);
text(ax1, 7.0, 6.35, '$Q_{\mathrm{in},2}=4$', ...
    'Interpreter','latex', 'FontSize',8);

% Annotation
text(ax1, 2.7, 2.3, '$Q_{\mathrm{total}} = 32+17.5+0.5+10=60$', ...
    'Interpreter','latex', 'FontSize',8);

formatAxes(ax1, TMax, zMax);

title(ax1, '\textbf{Before mixing}', ...
    'Interpreter','latex', 'FontSize',9.5);

%% =========================================================
% RIGHT: After mixing
% =========================================================
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');

% Main mixed temperature profile
plot(ax2, [1 6], [0 5], 'k-', 'LineWidth',0.7);
plot(ax2, [6 6], [5 6], 'k-', 'LineWidth',0.7);
addHatch(ax2, [0 6 6 0], [5 5 6 6], 'horizontal', 0.35);

% New energy block
fill(ax2, [0 8 10 0], [6 6 10 10], 'w', ...
    'EdgeColor','k', 'FaceAlpha',0.35, 'LineWidth',0.7);
addHatch(ax2, [0 8 10 0], [6 6 10 10], 'diagonal', 0.35);

text(ax2, 1.7, 8, '$Q_{\mathrm{new}}=32+4=36$', ...
    'Interpreter','latex', 'FontSize',8);

% Inlet part after mixing, dashed
plot(ax2, [0 6 6 0 0], [6 6 5 5 6], 'k--', 'LineWidth',0.7);
text(ax2, 2.4, 5.45, '$Q_{\mathrm{in},1}=6$', ...
    'Interpreter','latex', 'FontSize',8);

% Lower contribution
text(ax2, 0.2, 4.2, '$Q_{\mathrm{bottom}}=17.5$', ...
    'Interpreter','latex', 'FontSize',8);

% Annotation
text(ax2, 3.9, 2.3, '$Q_{\mathrm{total}}=36+6+17.5=59.5$', ...
    'Interpreter','latex', 'FontSize',8);

formatAxes(ax2, TMax, zMax);

title(ax2, '\textbf{After free convection mixing}', ...
    'Interpreter','latex', 'FontSize',9.5);

% =========================================================
% Custom legend with hatched patches
% =========================================================

% Create invisible axes for legend patches (clean separation)
axLeg = axes('Position',[0.15 0.02 0.15 0.12]); % bottom side
axis(axLeg,'off'); hold(axLeg,'on');

y0 = 0.8;
dy = 0.4;
xPatch = [0.1 0.4];   % statt [0.02 0.07]
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

version = 'v2';
export_png_dpi = 900;
dateTag = datestr(now, 'yyyymmdd');
baseName = sprintf('%s_illustration_fk_replacement_academic_%s', dateTag, version);

out_pdf = fullfile(RESULTS_SNAPSHOTS, [baseName '.pdf']);

exportgraphics(fig, out_pdf, 'ContentType', 'vector', 'Resolution', export_png_dpi);

disp(['Saved PDF: ' out_pdf]);


%% =========================================================
% Local helper functions
% =========================================================
function formatAxes(ax, TMax, zMax)
    axis(ax, 'equal');
    xlim(ax, [0 TMax]);
    ylim(ax, [0 zMax]);

    xlabel(ax, '$T$', 'Interpreter','latex', 'FontSize',9);
    ylabel(ax, '$z$', 'Interpreter','latex', 'FontSize',9);

    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 8;
    ax.LineWidth = 0.8;
    ax.XTick = [0 6 10];
    ax.YTick = [0 6 10];

    ax.XTickLabel = {'$0$', '$T(z_{\mathrm{in}})$', '$T(z_{\mathrm{mix}})$'};
    ax.YTickLabel = {'$0$', '$z_{\mathrm{in}}$', '$z_{\mathrm{mix}}$'};

    % --- Minor ticks (for dense grid) ---
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    
    % define spacing manually (key part!)
    ax.XAxis.MinorTickValues = 0:1:10;
    ax.YAxis.MinorTickValues = 0:1:10;
    
    % --- Grid ---
    grid(ax, 'on');        % major grid
    grid(ax, 'minor');     % minor grid
    
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

    % Border
    fill(ax, x, y, 'w', ...
        'EdgeColor','k', ...
        'LineWidth',0.8);

    % Hatch inside
    addHatch(ax, x, y, mode, 0.06); % smaller spacing for legend
end