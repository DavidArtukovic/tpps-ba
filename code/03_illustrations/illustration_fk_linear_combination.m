%% illustration_fk_linear_combination.m
% Visualization of algebraic problem, when energy is to small.
% Units: temperature = 1 unit, height = 1 unit

clear; close all; clc;

figure('Color','w','Position',[100 10 1100 1000]);

%% =========================================================
% UPPER LEFT: Decomposition into old energy + inlet contribution
% =========================================================
subplot(2,2,1); hold on; axis equal;
title('Extreme Concave Temperature Profile and New Water Before Mixing');
xlabel('T'); ylabel('z');

% 1) Set numeric ticks from 1 to 10
xticks(1:10);

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

% 4) Apply labels
xticklabels(labels);

yticks(1:10);

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

yticklabels(labelsY);

% --- Geometry (from sketch) ---
z_max = 11;
T_max = 11;

% Black and Gray Lines
plot([0 5],[0 5],'k','LineWidth',1.5);
plot([5 5],[5 7.9],'k','LineWidth',1.5);
plot([5 8],[7.9 7.9],'k','LineWidth',1.5);
plot([0 8],[8 8],'k','LineWidth',1.5);
plot([8 8],[7.9 8],'k','LineWidth',1.5);
plot([5 10],[5 10],'LineWidth',1.5, 'LineStyle','--', 'Color',[0.6 0.6 0.6]);
plot([0 10],[10 10],'LineWidth',1.5, 'LineStyle','--', 'Color',[0.6 0.6 0.6]);
text(6.8, 9.5,'Extrapolation','Color',[0.6 0.6 0.6],'FontSize',10);

annotation('arrow', ...
    [0.34 0.28], ...
    [0.83 0.87], ...
    'LineWidth',1);
text(1.9, 9.4,'Infinitesimal height','FontSize',8);


% Old energy block (purple) Q_old = 19.5
fill([0 5 5 0], [5 5 8 8],[0.8 0.6 0.9],'EdgeColor','k', 'FaceAlpha', 0.4);
text(1.2, 6.7,'Q_{old} = 15','Color',[0.4 0 0.6],'FontSize',10);

% Middle trapezoid (gray) = 12
fill([0 6 6],[0 7 7],[0.85 0.85 0.85],'EdgeColor','none');
text(0.7,3.6,'Q_{bottom} \approx 12.5','FontSize',10);

% Bottom triangle (orange) Q_out = 1/2
fill([0 0 0.2],[0 0.2 0.2],[1 0.7 0.4],'EdgeColor','k', 'FaceAlpha', 0.4);
text(0.6,0.3,'Q_{out}=1/50 \approx 0','Color',[0.8 0.4 0],'FontSize',9);

% Inlet energy (blue hatched area) Qin_2 =5
fill([5 10 10 5],[5.2 5.2 5 5],[0.6 0.8 1],'EdgeColor','k', 'FaceAlpha', 0.4);
text(7.1,5.4,'Q_{in,2} = 1','Color',[0 0.3 0.8],'FontSize',10);

% Inlet energy (blue hatched area) Qin_1=5
fill([0 5 5 0],[5.2 5.2 5 5],[0.6 0.8 1],'EdgeColor','k', 'FaceAlpha', 0.4);
text(1.4,5.4,'Q_{in,1} = 1','Color',[0 0.3 0.8],'FontSize',10);

% Axes & annotation
xlim([0 T_max]); ylim([0 z_max]);
text(4, 2.5,'12.5 + 15 + 2 = 29.5','FontSize',10);

grid on;

%% =========================================================
% UPPER RIGHT: Concave extrapolation case (3.7.6)
% =========================================================
subplot(2,2,2); hold on; axis equal;
title('Extrapolated Linear Approach Profile After Mixing');
xlabel('T'); ylabel('z');

% Axes
xticks(1:10); yticks(1:10);
xticklabels(labels); yticklabels(labelsY);
xlim([0 T_max]); ylim([0 z_max]);

% Old profile (concave: constant temperature)
% Main diagonal
plot([0.2 5],[0 4.8],'k','LineWidth',1.5);
% vertical part
plot([5 5],[4.8 5],'k','LineWidth',1.5);

% New profile in actual mixing zone
plot([3.333 7.333],[5 8],'k','LineWidth',1.5);  % T_w^*(z)
plot([3.333 5],[5 5],'k','LineWidth',1.5);

% Geometric extrapolation (infinitesimal slope)
plot([7.333 10],[8 10],'--','Color',[0.6 0.6 0.6],'LineWidth',1.5);
plot([0 10],[10 10],'LineWidth',1.5, 'LineStyle','--', 'Color',[0.6 0.6 0.6]);
text(6.8, 9.5,'Extrapolation','Color',[0.6 0.6 0.6],'FontSize',10);

% New energy block (Q_new = 16)
fill([0 3.333 7.333 0],[5 5 8 8],[0.6 0.9 0.6],'EdgeColor','k','FaceAlpha',0.45);
text(0.0,7.5,'Q_{new}=Q_{old}+Q_{in,2}=15+1=16','Color',[0.1 0.4 0.1],'FontSize',10);

% Vertical dashed line 
plot([5 5 ],[5 8],'k--');

% Middle rectangle (dashed) = 6
plot([0 5 5 0],[5 5 4.8 4.8],'k--');
text(1.6,4.5,'Q_{in,1}=1','FontSize',10);

% Lower triangle = 12.5
text(0.2,3.4,'Q_{bottom} \approx 12.5','FontSize',10);

% Axes
xlim([0 11]); ylim([0 11]);
grid on;

text(4, 2.5,'12.5 + 15 + 2 = 29.5','FontSize',10);

%% =========================================================
% LOWER LEFT: Uniform shift
% =========================================================
% Old concave temperature profile (polyline)
T_old = [5 5 8 8];
z_old = [5 7.9 7.9 8];

z_in  = 5;
z_mix = 8;
Qin2  = 1;

Dz = z_mix - z_in;
dT_fm = Qin2 / Dz;   % uniform temperature shift


subplot(2,2,3); hold on; axis equal;

title('Uniform Shift Approach Profile After Mixing');
xlabel('T'); ylabel('z');

% Axes
xticks(1:10); yticks(1:10);
xticklabels(labels); yticklabels(labelsY);
xlim([0 T_max]); ylim([0 z_max]);

% Lower Part
% Main diagonal
plot([0.2 5],[0 4.8],'k','LineWidth',1.5);
% vertical part
plot([5 5],[4.8 5],'k','LineWidth',1.5);

% New temperature profile
plot([5+dT_fm 5+dT_fm],[5 7.9],'k','LineWidth',1.5);
plot([5+dT_fm 8+dT_fm],[7.9 7.9],'k','LineWidth',1.5);
plot([5 5+dT_fm],[5 5],'k','LineWidth',1.5);

plot([0 8+dT_fm],[8 8],'k','LineWidth',1.5);
plot([8+dT_fm 8+dT_fm],[7.9 8],'k','LineWidth',1.5);

% Old Profile dashed line 
plot([5 5 ],[5 7.9],'LineStyle','--', 'Color', [0.6 0.6 0.6], 'LineWidth',0.3);
plot([5 5+dT_fm ],[7.9 7.9], 'LineStyle','--', 'Color', [0.6 0.6 0.6], 'LineWidth',0.3);

% New energy block (Q_new = 16)
fill([0 5+dT_fm 5+dT_fm 8+dT_fm 8+dT_fm 0],[5 5 7.9 7.9 8 8],[0.6 0.9 0.6],'EdgeColor','k','FaceAlpha',0.45);
text(0.0,7.5,'Q_{new}=Q_{old}+Q_{in,2}=15+1=16','Color',[0.1 0.4 0.1],'FontSize',10);

% Middle rectangle (dashed) = 5
plot([0 5 5 0],[5 5 4.8 4.8],'k--');
text(1.6,4.5,'Q_{in,1}=1','FontSize',10);

% Lower triangle = 12.5
text(0.2,3.4,'Q_{bottom} \approx 12.5','FontSize',10);

% Plot axes limits
z_max = 11;
T_max = 11;
xlim([0 T_max]); ylim([0 z_max]);

grid on;

%% =========================================================
% LOWER RIGHT: lambda-weighted profile
% =========================================================
subplot(2,2,4); hold on; axis equal;
title('\lambda-Weighted Profile');
xlabel('T'); ylabel('z');

% same axes and ticks
xticks(1:10); yticks(1:10);
xticklabels(labels); yticklabels(labelsY);
xlim([0 T_max]); ylim([0 z_max]);


% ---------- FK profile (only mixing zone) ----------
% FK line goes from (T,z) = (3.333,5) to (7.333,7.9)
T_fk = T_old;
mask_mix = (z_old >= 5);

T_fk(mask_mix) = interp1( ...
    [5 8], ...
    [3.333 7.333], ...
    z_old(mask_mix));

% ---------- fully mixed profile (reuse logic) ----------
T_fm = T_old;
T_fm(mask_mix) = T_old(mask_mix) + dT_fm;

% ---------- compute lambda from boundary constraints ----------
T_old_in  = 5;
T_old_mix = 8;

% FK boundary values (known analytically)
Tfk_in  = 3.333;
Tfk_mix = 7.333;

% Fully mixed boundary values
Tfm_in  = 5 + dT_fm;
Tfm_mix = 8 + dT_fm;

lambda_in  = (Tfm_in  - T_old_in ) / (Tfm_in  - Tfk_in);
lambda_mix = (Tfm_mix - T_old_mix) / (Tfm_mix - Tfk_mix);

lambda = min([lambda_in, lambda_mix, 1]);
lambda = max(lambda,0);

% ---------- lambda blend ----------
T_new = T_old;
T_new(mask_mix) = lambda*T_fk(mask_mix) + (1-lambda)*T_fm(mask_mix);

% ---------- plot new profiles ----------
plot(T_new, z_old, 'k','LineWidth',1.5);      % lambda result
plot(T_fk,  z_old, 'k-.','LineWidth',1.0);    % FK reference
plot(T_fm,  z_old, 'k:','LineWidth',1.0);     % fully mixed reference

% additional lines
plot([5 5.333], [5 5], 'k:','LineWidth',1.0)
plot([0 8.16664], [8 8], 'k','LineWidth',1.5);      % lambda result

% Middle rectangle (dashed) = 5
plot([0 5 5 0],[5 5 4.8 4.8],'k--');
text(1.6,4.5,'Q_{in,1}=1','FontSize',10);

% Lower Part
% Main diagonal
plot([0.2 5],[0 4.8],'k','LineWidth',1.5);
% vertical part
plot([5 5],[4.8 5],'k','LineWidth',1.5);

% Lower triangle = 12.5
text(0.2,3.4,'Q_{bottom} \approx 12.5','FontSize',10);


% New energy block (Q_new = 16)
fill([0 5 5.6443 8.16664 8.16664 0],[5 5 7.9 7.9 8 8],[0.6 0.9 0.6],'EdgeColor','k','FaceAlpha',0.45);
text(0.0,7.5,'Q_{new}=Q_{old}+Q_{in,2}=15+1=16','Color',[0.1 0.4 0.1],'FontSize',10);

% indicators for boundary temperatures
plot([5 5],[0 5],'r--', 'LineWidth',0.5);
plot([8 8],[0 8],'r--', 'LineWidth',0.5);


% Add legend to the plot for clarity
legend('Weighted Profile', 'Linear Approach Profile', 'Uniform Shift Profile', 'Location', 'Best');
grid on;
