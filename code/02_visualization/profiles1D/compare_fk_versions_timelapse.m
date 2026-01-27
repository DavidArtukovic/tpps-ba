clc
clear
close all
%%%------------------------------------------%%%
% 01. Load Scenario and Initialization Simulation
%%%------------------------------------------%%%

% Load local paths (per-machine config)
run(fullfile('..','..','..', 'configs', 'paths_local.m'));

RESULTS_VIS_BASE   = fullfile(RESULTS_BASE, "02_visualizations");
RESULTS_TIMELAPSE  = fullfile(RESULTS_VIS_BASE, "timelapse");

if ~exist(RESULTS_TIMELAPSE, "dir")
    mkdir(RESULTS_TIMELAPSE);
end

% Build data subfolder for this configuration
DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_FK  = fullfile(DATA_BASE, 'scenario1_freeConv');
DATA_SCEN1  = fullfile(DATA_BASE, 'scenario1v');
%%%------------------------------------------%%%
% 02. Load relevant data sets
%%%------------------------------------------%%%

load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));
data_no_fk   = load(fullfile(DATA_SCEN1_BASE, 'd18_h18_Res_Matlab_d18_18.mat'));
data_comsol  = load(fullfile(DATA_SCEN1_BASE, '1D_05_TPPS_18_18.mat'));
data_fk_v1   = load(fullfile(DATA_SCEN1_OWN,  '260108_d18_h18_Res_Matlab_FK_v4.mat'));
data_fk_v2   = load(fullfile(DATA_SCEN1_OWN,  '260109_d18_h18_Res_Matlab_FK_v6.mat'));

%%%------------------------------------------%%%
% 03. Reduce Matlab solutions to hourly resolution
%%%------------------------------------------%%%

k = 0;
for i = 1:1920
    if mod(i,4) == 0
        k = k + 1;
        T_Sys_M(:,k) = data_no_fk.Res_System_d18_18(401:end,i);
    end
end

k = 0;
for i = 1:1920
    if mod(i,4) == 0
        k = k + 1;
        T_Sys_M_fk_v1(:,k) = data_fk_v1.Res_System_d18_18_FK(401:end,i);
    end
end

k = 0;
for i = 1:1920
    if mod(i,4) == 0
        k = k + 1;
        T_Sys_M_fk_v2(:,k) = data_fk_v2.Res_System_d18_18_FK(401:end,i);
    end
end

%%%------------------------------------------%%%
% 04. Spatial alignment MATLAB vs COMSOL
%%%------------------------------------------%%%

dz   = 0.05;
dz_M = 0.005;

h_M     = (length(T_Sys_M)-1)*dz_M;
h_Sys_C = (length(data_comsol.z)-1)*dz;
h_diff  = abs(h_Sys_C - h_M);

z_Mo = data_comsol.z(1)   - h_diff/2;
z_Mu = data_comsol.z(end) + h_diff/2;
z_M  = flip(z_Mu:dz_M:z_Mo)';

%%%------------------------------------------%%%
% 05. Map Matlab time steps to COMSOL output times
%%%------------------------------------------%%%

o = 0;
l = 0;

for i = 1:480
    if t_hour(3,i) > 0
        o = o + 1;
        l = l + t_hour(3,i);
        t_Sys_M(1,o) = l;
    end
end

for i = 1:length(t_Sys_M)
    T_M(:,i)        = T_Sys_M(:,t_Sys_M(1,i));
    T_M_fk_v1(:,i)  = T_Sys_M_fk_v1(:,t_Sys_M(1,i));
    T_M_fk_v2(:,i)  = T_Sys_M_fk_v2(:,t_Sys_M(1,i));
end

%%%------------------------------------------%%%
% 06. Assemble COMSOL system temperature
%%%------------------------------------------%%%

for i = 1:length(t_Sys_M)
    zp(1,i) = 1 + (data_comsol.z(1) - data_comsol.z_p(1,i)) / dz;
end

T_Sys_C = zeros(length(data_comsol.T1_W), length(t_Sys_M));

for i = 1:length(t_Sys_M)
    T_Sys_C(:,i) = data_comsol.T1_W(:,i);
    T_Sys_C(zp(1,i)+1 : zp(1,i)+length(data_comsol.z_p)-2 , i) = ...
        data_comsol.T1_P(2:end-1,i);
end

%%%------------------------------------------%%%
% 07. Collect models for plotting
%%%------------------------------------------%%%

models(1).name  = 'COMSOL';
models(1).T     = T_Sys_C;
models(1).z     = data_comsol.z;
models(1).style = '-';
models(1).color = [0.3010 0.7450 0.9330];

models(2).name  = 'MATLAB (no FK)';
models(2).T     = T_M;
models(2).z     = z_M;
models(2).style = '-';
models(2).color = [0.8500 0.3250 0.0980];

models(3).name  = 'MATLAB (FK v4)';
models(3).T     = T_M_fk_v1;
models(3).z     = z_M;
models(3).style = '--';
models(3).color = [0.4660 0.6740 0.1880];

models(4).name  = 'MATLAB (FK v6)';
models(4).T     = T_M_fk_v2;
models(4).z     = z_M;
models(4).style = '--';
models(4).color = [0.4940 0.1840 0.5560];


%%%------------------------------------------%%%
% Video writer initialization
%%%------------------------------------------%%%

videoname = fullfile(RESULTS_TIMELAPSE, ...
    'timelapse_1D_temperature_profiles.mp4');

v = VideoWriter(videoname,'MPEG-4');
v.FrameRate = 10;
open(v);

%%%------------------------------------------%%%
% 08. Dynamic time-lapse visualization
%%%------------------------------------------%%%

figure('Color','w','Name','1D temperature profiles – time-lapse');
hold on
grid on

for mm = 1:numel(models)
    h(mm) = plot(models(mm).T(:,1), models(mm).z, ...
        models(mm).style, ...
        'Color', models(mm).color, ...
        'LineWidth',1.5);
end

xlabel('Temperature in TPPS [$^{\circ}$C]','Interpreter','Latex')
ylabel('System depth [m]','Interpreter','Latex')
legend({models.name},'Interpreter','Latex','Location','best')

axis([40 80 z_Mu z_Mo])

% Piston rectangle (initialized)
z_piston = data_comsol.z_p(end,1);

h_rect = rectangle('Position',[40 z_piston 40 18], ...
                   'FaceColor',[1 1 1 0.1],'EdgeColor',[0 0 0]);
uistack(h_rect,'bottom');   % <-- IMPORTANT
title_handle = title('', 'Interpreter','Latex');

%%%------------------------------------------%%%
% 09. Animation loop
%%%------------------------------------------%%%

for i = 1:length(t_Sys_M)

    for mm = 1:numel(models)
        set(h(mm), 'XData', models(mm).T(:,i), ...
                   'YData', models(mm).z);
    end

    z_piston = data_comsol.z_p(end,i);

    set(h_rect, 'Position', [40 z_piston 40 18]);

    title_handle.String = ...
        sprintf('1D temperature profiles - t = %d h', i)

    drawnow
    writeVideo(v, getframe(gcf));
    pause(0.1)
end
close(v);