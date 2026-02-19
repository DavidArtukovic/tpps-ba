%% plot_fk_update_illustration.m
% -------------------------------------------------------------------------
% Illustration of discrete free-convection (FK) update in TPPS model
%
% Left : Temperature profile T(z) at multiple time levels
% Right: Geometrically scaled TPPS cross-section with temperature coloring
%
% The blue patch indicates the discrete FK update region
% (z_lower_bypass <= z <= z_membrane).
%
% Author: David Artukovic
% -------------------------------------------------------------------------

clear; clc; close all;

%% ------------------------------------------------------------------------
% 1) Geometric parameters (all lengths in meters)
% -------------------------------------------------------------------------
H_tot = 10.0;              % total water height

z_inlet        = 10.0;
z_upper_bypass = 8.5;
z_upper_rep    = 7.0;
z_membrane     = 5.5;
z_lower_bypass = 4.0;
z_lower_rep    = 2.0;
z_outlet       = 0.0;

z_marks = [ ...
    z_inlet
    z_upper_bypass
    z_upper_rep
    z_membrane
    z_lower_bypass
    z_lower_rep
    z_outlet ];

z_labels = { ...
    'z_{inlet}'
    'z_{upper bypass}'
    'z_{upper rep}'
    'z_{membrane}'
    'z_{lower bypass}'
    'z_{lower rep}'
    'z_{outlet}' };


Nz = 500;
z  = linspace(z_outlet, z_inlet, Nz).';

% Relevant indices
idx_inlet = z(end);
[~, idx_upper_bypass] = min(abs(z - 8.5));



T1 = zeros(size(z));

% --- Segment D: lower volume (almost isothermal) ---
T1(z <= z_lower_rep) = 45;

% --- Segment C: weak gradient (lower bypass region) ---
idx_C = z > z_lower_rep & z <= z_inlet;
T1(idx_C) = linspace(45, 65, sum(idx_C));

T_inlet = 80;
T1(end)= T_inlet;
T_upper_bypass = T1(idx_upper_bypass);

%% ------------------------------------------------------------------------
% 1) Temperature profile T2
% -------------------------------------------------------------------------
% z  -> column vector of heights (same grid as T profile)
% T1 -> original temperature profile [°C]

T2 = zeros(size(T1));     % new profile
dT = zeros(size(T1));     % temperature shift ΔT(z)


idx_low = z <= z_upper_bypass;
dT(idx_low) = 0.5;


idx_top = z > z_upper_bypass;


% Boundary shifts
dT1 = 4.0;                                   % shift at upper bypass

% Distance measured from the top downward
eta = (z_inlet - z(idx_top));  % 1 at bypass, 0 at inlet

k = 2;   % curvature control (larger = flatter near bypass, steeper near inlet)

% Monotonic convex exponential rise toward dT2
dT(idx_top) = (80 - 64) .* exp(-k * eta);
dT(end) = 0;
T2 = T1 + dT;
T2(idx_upper_bypass) = T2(idx_upper_bypass)+dT1;





%% ------------------------------------------------------------------------
% 1) Region below upper bypass: constant +0.5 °C
% -------------------------------------------------------------------------
idx_low = z <= z_upper_bypass;
dT(idx_low) = 0.5;

%% ------------------------------------------------------------------------
% 2) Region between upper bypass and inlet: exponential rise
% -------------------------------------------------------------------------
idx_top = z > z_upper_bypass;

z1 = z_upper_bypass;
z2 = z_inlet;

% Boundary conditions for the shift:
dT1 = 4.0;                         % shift at upper bypass
T_target_top = 80;                 % final temperature at inlet
dT2 = T_target_top - T1(z == z2);  % required shift at inlet

% If z grid does not hit z_inlet exactly, use last value:
if isempty(dT2)
    dT2 = T_target_top - T1(end);
end

% Normalized height coordinate ξ in [0,1]
xi = (z(idx_top) - z1) / (z2 - z1);

% Exponential (concave) interpolation
k = 3;  % curvature parameter (>0 → concave increasing)

dT(idx_top) = dT1 + (dT2 - dT1) * (1 - exp(-k*xi)) / (1 - exp(-k));

%% ------------------------------------------------------------------------
% 3) Build new temperature profile
% -------------------------------------------------------------------------
T2 = T1 + dT;


%% ------------------------------------------------------------------------
% 5) FK update region as trapezoid (aligned to T^n)
% -------------------------------------------------------------------------
idx_fk = z >= z_lower_bypass & z <= z_membrane;

T_fk_top = Tn(idx_fk);
T_fk_bot = T_fk_top - 4.0;   % width of trapezoid (visual, schematic)

z_fk = z(idx_fk);

%% ------------------------------------------------------------------------
% 3) Figure & layout (shared z-axis)
% -------------------------------------------------------------------------
figure('Color','w','Position',[100 100 1200 700]);

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

%% ------------------------------------------------------------------------
% 4) LEFT: z–T diagram
% -------------------------------------------------------------------------
ax1 = nexttile;
hold on; box on; grid on;

plot(Tn,    z, 'k',  'LineWidth',2);
plot(Tn_dt, z, 'Color',[0.5 0.5 0.5],'LineWidth',2);
plot(Tn_2dt,z, 'Color',[0.8 0.8 0.8],'LineWidth',2);

% FK update region (blue patch)

% FK patch (drawn first)
patch([T_fk_top; flipud(T_fk_bot)], ...
      [z_fk;     flipud(z_fk)], ...
      [0.8 0.9 1.0], ...
      'EdgeColor','none','FaceAlpha',0.45);

% horizontal reference lines
for i = 1:numel(z_marks)
    yline(z_marks(i),'k:');
end

xlabel('Temperature T [^\circC]');
ylabel('Height z [m]');
title('Temperature evolution and FK update region');

legend({'T^n','T^{n+\Delta t}','T^{n+2\Delta t}','Discrete FK update'}, ...
       'Location','southoutside');

ylim([z_outlet z_inlet]);

%% ------------------------------------------------------------------------
% 5) RIGHT: TPPS cross-section (scaled in z)
% -------------------------------------------------------------------------
ax2 = nexttile;
hold on; box on;

% geometric widths (arbitrary but fixed)
w_core   = 1.0;
w_gap    = 0.3;
x0       = 0.0;

% temperature colormap
cmap = jet(256);
colormap(ax2,cmap);

% draw water column as colored rectangles
for i = 1:Nz-1
    dz_loc = z(i+1) - z(i);
    Tloc   = T_n(i);
    cidx   = round(1 + (Tloc-45)/(80-45)*(size(cmap,1)-1));
    cidx   = max(min(cidx,size(cmap,1)),1);
    
    rectangle('Position',[x0, z(i), w_core, dz_loc], ...
              'FaceColor',cmap(cidx,:), ...
              'EdgeColor','none');
end

% membrane
plot([x0 x0+w_core],[z_membrane z_membrane],'k','LineWidth',2);
text(x0+w_core+0.05,z_membrane,'membrane','VerticalAlignment','bottom');

% bypass openings
plot([x0+w_core x0+w_core+w_gap],[z_upper_bypass z_upper_bypass],'k','LineWidth',2);
plot([x0+w_core x0+w_core+w_gap],[z_lower_bypass z_lower_bypass],'k','LineWidth',2);

text(x0+w_core+w_gap+0.05,z_upper_bypass,'upper bypass','VerticalAlignment','middle');
text(x0+w_core+w_gap+0.05,z_lower_bypass,'lower bypass','VerticalAlignment','middle');

% residual energy box
rectangle('Position',[x0+w_core+w_gap+0.6, z_lower_bypass+0.3, 1.6, 1.0], ...
          'EdgeColor','k','LineWidth',1.5);
text(x0+w_core+w_gap+0.65,z_lower_bypass+0.9, ...
    {'Residual thermal energy','stored for discrete FK update'}, ...
    'FontSize',9);

annotation('arrow', ...
    [0.77 0.85],[0.55 0.55]);

axis equal;
xlim([-0.2 4.0]);
ylim([z_outlet z_inlet]);

xlabel('Radial direction (schematic)');
title('TPPS cross-section and algorithmic interpretation');

cb = colorbar;
cb.Label.String = 'Temperature [^\circC]';

linkaxes([ax1 ax2],'y');

%% ------------------------------------------------------------------------
% End of script
% -------------------------------------------------------------------------
