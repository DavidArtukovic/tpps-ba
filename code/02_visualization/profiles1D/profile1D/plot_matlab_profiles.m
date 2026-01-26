clc
clear
close all
%%%------------------------------------------%%%
% 01. Load paths and data
%%%------------------------------------------%%%

run(fullfile('..','..','..', 'configs', 'paths_local.m'));

RESULTS_VIS_BASE   = fullfile(RESULTS_BASE, "02_visualizations");
RESULTS_PROFILES_1D = fullfile(RESULTS_VIS_BASE, "profiles1D");

DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_OWN  = fullfile(DATA_BASE, 'scenario1_freeConv');

load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));
load(fullfile(DATA_BASE, 'scenario1/SzenarioComsol.mat')) % for plotting the piston

data_no_fk = load(fullfile(DATA_SCEN1_BASE, ...
    'd18_h18_Res_Matlab_d18_18.mat'));

data_fk_v6 = load(fullfile(DATA_SCEN1_OWN, ...
    '260116_d18_h18_Res_Matlab_FK_v11.mat'));

data_fk_v7 = load(fullfile(DATA_SCEN1_OWN, ...
    '260117_d18_h18_Res_Matlab_FK_v13.mat'));

%%%------------------------------------------%%%
% 02. Spatial grid (MATLAB only)
%%%------------------------------------------%%%

dz_M = 0.005;                           % MATLAB resolution
Nz   = size(data_no_fk.Res_System_d18_18,1) - 400;
z_M_star  = flip((0:Nz-1)' * dz_M);
z_M = -1*(37.33 - z_M_star)+0.665;
inlet_idx = 4288;
%%
 
%%%------------------------------------------%%%
% 03. Extract full 15-min resolution (no reduction)
%%%------------------------------------------%%%

T_noFK = data_no_fk.Res_System_d18_18(401:end,:);
T_fk6  = data_fk_v6.Res_System_d18_18_FK(401:end,:);
T_fk7  = data_fk_v7.Res_System_d18_18_FK(401:end,:);

%%%------------------------------------------%%%
% 04. Time axis (14 min per step)
%%%------------------------------------------%%%

dt_min = 15;
t_hours = (0:size(T_noFK,2)-1) * dt_min / 60;

%%%------------------------------------------%%%
% 05. Select operating points (indices!)
%%%------------------------------------------%%%
% Example: every 6 hours
p = round([31 31.25 31.5 38.5 54] * 60 / dt_min);

%%%------------------------------------------%%%
% 06. Model container
%%%------------------------------------------%%%

models(1).name  = 'MATLAB (no FK)';
models(1).T     = T_noFK;
models(1).style = '-';
models(1).color = [0.8500 0.3250 0.0980];

models(2).name  = 'MATLAB (FK v11)';
models(2).T     = T_fk6;
models(2).style = '--';
models(2).color = [0.4660 0.6740 0.1880];

models(3).name  = 'MATLAB (FK v13)';
models(3).T     = T_fk7;
models(3).style = '--';
models(3).color = [0.4940 0.1840 0.5560];

%%%------------------------------------------%%%
% 07. Plot
%%%------------------------------------------%%%

figure('Color','w','Name','1D temperature profiles (15 min resolution)');

for m = 1:length(p)
    subplot(1,length(p),m)
    hold on
    
    % Physical piston position (lower edge)
    z_piston = t_900(2,p(m))*(-19)+(1-t_900(2,p(m)))*(-36);
    disp(z_piston);
    rectangle('Position', [40, z_piston, 40, 18], ...
              'FaceColor', [1 1 1 0.1], ...
              'EdgeColor', [0 0 0]);


    for mm = 1:numel(models)
        plot(models(mm).T(:,p(m)), z_M, ...
            models(mm).style, ...
            'Color', models(mm).color,...
            'LineWidth',1);
    end

    grid on
    axis([40 80 min(z_M) max(z_M)])
    xticks([40 50 60 70 80])
    yline(z_M(inlet_idx),'k--','LineWidth',1);

    title(sprintf('t = %.2f h', t_hours(p(m))), ...
        'Interpreter','Latex','FontSize',10)

    if m == 1
        ylabel('System depth [m]', ...
            'Interpreter','Latex','FontSize',11)
        xlabel('Temperature [$^{\circ}$C]', ...
            'Interpreter','Latex','FontSize',11)
    end

    if m == length(p)
        legend({models.name}, ...
            'Interpreter','Latex','FontSize',9, ...
            'Location','best')
    end
end

sgtitle('1D temperature profiles (15 min resolution)', ...
    'Interpreter','Latex','FontSize',12)
