%% illustration_fk_extrapolation.m
% Visualization of energy conservation using area arguments
% Units: temperature = 1 unit, height = 1 unit

clear; close all; clc;

figure('Color','w','Position',[100 100 1100 450]);

%% =========================================================
% LEFT: Decomposition into old energy + inlet contribution
% =========================================================
subplot(1,2,1); hold on; axis equal;
title('Left: Temperature Profile and New Water Before Mixing');
xlabel('T'); ylabel('z');

% 1) Set numeric ticks from 1 to 10
xticks(1:10);

% 2) Default: empty labels
labels = repmat({''}, 1, 10);

% 3) Add labels only at relevant positions
labels{5}  = 'T(z_{in})';
labels{8}  = 'T(z_{mix})';
labels{10}  = 'T_{in}';

% 4) Apply labels
xticklabels(labels);

yticks(1:10);

labelsY = repmat({''}, 1, 10);
labelsY{5} = 'z_{in}';
labelsY{8} = 'z_{mix}';
labelsY{10} = 'z_{mix}^*';

yticklabels(labelsY);

% --- Geometry (from sketch) ---
z_max = 11;
T_max = 11;

% Black and Gray Lines
plot([0 8],[0 8],'k','LineWidth',1.5);
plot([0 8],[8 8],'k','LineWidth',1.5);
plot([8 10],[8 10],'LineWidth',1.5, 'LineStyle','--', 'Color',[0.6 0.6 0.6]);
plot([0 10],[10 10],'LineWidth',1.5, 'LineStyle','--', 'Color',[0.6 0.6 0.6]);
text(6.8, 9.5,'Extrapolation','Color',[0.6 0.6 0.6],'FontSize',10);

% Old energy block (purple) Q_old = 19.5
fill([0 5 8 0], [5 5 8 8],[0.8 0.6 0.9],'EdgeColor','k', 'FaceAlpha', 0.4);
text(1.2, 6.7,'Q_{old} = 19.5','Color',[0.4 0 0.6],'FontSize',10);

% Middle trapezoid (gray) = 12
fill([0 6 6],[0 7 7],[0.85 0.85 0.85],'EdgeColor','none');
text(0.7,3.6,'Q_{bottom} = 12','FontSize',10);

% Bottom triangle (orange) Q_out = 1/2
fill([0 0 1],[0 1 1],[1 0.7 0.4],'EdgeColor','k', 'FaceAlpha', 0.4);
text(0.6,0.3,'Q_{out}=1/2','Color',[0.8 0.4 0],'FontSize',9);

% Inlet energy (blue hatched area) Qin_2 =5
fill([5 10 10 5],[6 6 5 5],[0.6 0.8 1],'EdgeColor','k', 'FaceAlpha', 0.4);
text(7.1,5.4,'Q_{in,2} = 5','Color',[0 0.3 0.8],'FontSize',10);

% Inlet energy (blue hatched area) Qin_1=5
fill([0 5 5 0],[6 6 5 5],[0.6 0.8 1],'EdgeColor','k', 'FaceAlpha', 0.4);
text(1.4,5.4,'Q_{in,1} = 5','Color',[0 0.3 0.8],'FontSize',10);

% Axes & annotation
xlim([0 T_max]); ylim([0 z_max]);
text(4, 2.5,'12 + 19.5 + 0.5 + 10 = 42','FontSize',10);

grid on;

%% =========================================================
% RIGHT: Equivalent energy after extrapolated mixing
% =========================================================
subplot(1,2,2); hold on; axis equal;
title('Right: Temperature Profile and New Water After Mixing');
xlabel('T'); ylabel('z');


% 1) Set numeric ticks from 1 to 10
xticks(1:10);

% 2) Default: empty labels
labels = repmat({''}, 1, 10);

% 3) Add labels only at relevant positions
labels{5}  = 'T(z_{in})';
labels{8}  = 'T(z_{mix})';
labels{10}  = 'T_{in}';

% 4) Apply labels
xticklabels(labels);

yticks(1:10);

labelsY = repmat({''}, 1, 10);
labelsY{5} = 'z_{in}';
labelsY{8} = 'z_{mix}';
labelsY{10} = 'z_{mix}^*';

yticklabels(labelsY);

% Main diagonal
plot([1 5],[0 4],'k','LineWidth',1.5);

% vertical part
plot([5 5],[4 5],'k','LineWidth',1.5);

% New Extrapolated Part
plot([7.38, 10],[5 10],'LineWidth',1.5, 'LineStyle','-', 'Color',[0.6 0.6 0.6]);
plot([5, 7.38],[5 5],'LineWidth',1.5, 'LineStyle','-', 'Color',[0.6 0.6 0.6]);
plot([0, 10],[10 10],'LineWidth',1.5, 'LineStyle','-', 'Color',[0.6 0.6 0.6]);

% Old Part
plot([5, 8],[5 8], 'k', 'LineWidth',1.5, 'LineStyle','--');
plot([0, 8.95],[8 8],'k', 'LineWidth',1.5, 'LineStyle','--');

% Extrapolated energy block (green) Q_new = 36
fill([0 7.38 8.95 0], [5 5 8 8],[0.6 0.9 0.6],'EdgeColor','k', 'FaceAlpha', 0.4);
text(0.2,6.5,'Q_{new} = 19.5 + 5 = 24.5','Color',[0.1 0.4 0.1]  ,'FontSize',10);

% Middle rectangle (dashed) = 6
plot([0 5 5 0],[5 5 4 4],'k--');
text(1.6,4.5,'Q_{in,1}=5','FontSize',10);

% Lower triangle = 12.5
text(0.5,2.4,'Q_{bottom} = 12','FontSize',10);

% Axes & annotation
xlim([0 T_max]); ylim([0 z_max]);
text(4.3, 2.5,'24.5 + 5 + 12 = 41.5','FontSize',10);


annotation('arrow', ...
    [0.78 0.78], ...
    [0.55 0.45], ...
    'LineWidth',1.5);
text(5.5, 4,'Q_{in,2} = 5 within actual mix-zone','FontSize',10);

grid on;
