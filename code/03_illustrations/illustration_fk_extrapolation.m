%% plot_energy_conservation_boxes.m
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

% Desired tick positions
xticks([5 8 10])
yticks([5 8 10])
% Custom tick labels
xticklabels({'T(z_{in})','T(z_{mix})','T_{in}'})
yticklabels({'z_{in}','z_{mix}','z_{mix}^*'})

% --- Geometry (from sketch) ---
z_max = 11;
T_max = 11;

% Main diagonal (black)
plot([0 8],[0 8],'k','LineWidth',1.5);
plot([0 8],[8 8],'k','LineWidth',1.5);

% Old energy block (purple) Q_old = 19.5
fill([0 5 8 0], [5 5 8 8],[0.8 0.6 0.9],'EdgeColor','k', 'FaceAlpha', 0.4);
text(2.2,8.5,'Q_{old} = 32','Color',[0.4 0 0.6],'FontSize',10);


% Middle trapezoid (gray) = 17.5
fill([0 6 6],[0 7 7],[0.85 0.85 0.85],'EdgeColor','none');
text(0.7,3.6,'Q_{bottom} = 17.5','FontSize',10);

% Bottom triangle (orange) Q_out = 1/2
fill([0 0 1],[0 1 1],[1 0.7 0.4],'EdgeColor','k', 'FaceAlpha', 0.4);
text(0.4,0.3,'Q_{out}=1/2','Color',[0.8 0.4 0],'FontSize',9);

% Inlet energy (blue hatched area) Qin_2 ≈ 4
fill([5 10 10 5],[6 6 5 5],[0.6 0.8 1],'EdgeColor','k', 'FaceAlpha', 0.4);
text(7.1,6.4,'Q_{in,2} \approx 4','Color',[0 0.3 0.8],'FontSize',10);

% Inlet energy (blue hatched area) Qin_1 ≈ 6
fill([0 6 6 0],[7 7 6 6],[0.6 0.8 1],'EdgeColor','k', 'FaceAlpha', 0.4);
text(3.1,6.4,'Q_{in,1} \approx 6','Color',[0 0.3 0.8],'FontSize',10);

% Axes & annotation
xlim([0 T_max]); ylim([0 z_max]);
text(4, 2.5,'12 + 17.5 + 0.5 + 10 = 60','FontSize',10);

grid on;

%% =========================================================
% RIGHT: Equivalent energy after extrapolated mixing
% =========================================================
subplot(1,2,2); hold on; axis equal;
title('Left: Temperature Profile and New Water After Mixing');
xlabel('T'); ylabel('z');

% Main diagonal
plot([1 6],[0 5],'k','LineWidth',1.5);

% vertical part
plot([6 6],[6 5],'k','LineWidth',1.5);

% Extrapolated energy block (green) Q_new = 36
fill([0 8 10 0], [6 6 10 10],[0.6 0.9 0.6],'EdgeColor','k', 'FaceAlpha', 0.4);
text(2.2,8.5,'Q_{new} = 32+4=36','Color',[0.1 0.4 0.1]  ,'FontSize',10);

% Middle rectangle (dashed) = 6
plot([0 6 6 0],[6 6 5 5],'k--');
text(2.6,5.5,'Q_{in,1}=6','FontSize',10);

% Lower triangle = 17.5
text(1.7,3.6,'Q_{bottom}=17.5','FontSize',10);

% Axes & annotation
xlim([0 T_max]); ylim([0 z_max]);
text(4, 2.5,'36 + 6 + 17.5 = 59.5','FontSize',10);

grid on;
