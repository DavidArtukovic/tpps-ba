clc
clear
close all
%%%------------------------------------------%%%
% 01. Load Scenario and Initialization Simulation
%%%------------------------------------------%%%

% Load local paths (per-machine config)
run(fullfile('..','..','..', 'configs', 'paths_local.m'));

RESULTS_VIS_BASE = fullfile(RESULTS_BASE, "02_visualizations");
RESULTS_PROFILES_1D = fullfile(RESULTS_VIS_BASE, "profiles1D");

%%
% Build data subfolder for this configuration
DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_OWN  = fullfile(DATA_BASE, 'scenario1_freeConv');

%%%------------------------------------------%%%
% 02. Load relevant data sets
%%%------------------------------------------%%%

load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));
data_no_fk = load(fullfile(DATA_SCEN1_BASE, 'd18_h18_Res_Matlab_d18_18.mat'));    % Matlab reference (no FK)
data_comsol = load(fullfile(DATA_SCEN1_BASE, '1D_05_TPPS_18_18.mat'));             % COMSOL reference
data_fk_v1 = load(fullfile(DATA_SCEN1_OWN,  '260124_d18_h18_Res_Matlab_FK_v16.mat')); % Matlab FK version
data_fk_v2 = load(fullfile(DATA_SCEN1_OWN,  '260206_d18_h18_Res_Matlab_FK_v17.mat')); % Matlab FK version

%%
%%%------------------------------------------%%%
% 03. Reduce Matlab solutions to hourly resolution
%%%------------------------------------------%%%

%%% Matlab No FK Version
indices_lower_volume = data_no_fk.Res_900_d18_18(7,:); % number of cells in lower volume

% take every 4th timestep (columns)
cols = 4:4:size(data_no_fk.Res_System_d18_18, 2);

% skip first 400 rows (insulation)
T_Sys_M = data_no_fk.Res_System_d18_18(401:end, cols);

      
% Part for Ringgap
blockLen = 389;

A        = data_no_fk.Res_Wasser_d18_18;     % [Nz x Nt]
startIdx = data_no_fk.Res_900_d18_18(7,:);   % [1 x Nt]

A_red        = A(:, cols);            % reduced time matrix
startIdx_red = startIdx(cols);        % matching start indices

nCols = length(cols);

% --- build row index matrix ---
rowOffsets = (0:blockLen-1)';                 % [389 x 1]
rows = startIdx_red(:)' + rowOffsets;         % [389 x Nt_reduced]

colsMat = repmat(1:nCols, blockLen, 1);       % [389 x Nt_reduced]
linInd = sub2ind(size(A_red), rows, colsMat);
T_ring_no_fk = A_red(linInd);                % [389 x Nt_reduced]


%%% Matlab FK Version 1

T_Sys_M_fk_v1 = data_fk_v1.Res_System_d18_18_FK(401:end,cols); % FK version 1

% Part for Ringgap
blockLen = 389;

A        = data_no_fk.Res_Wasser_d18_18;     % [Nz x Nt]
startIdx = data_no_fk.Res_900_d18_18(7,:);   % [1 x Nt]

A_red        = A(:, cols);            % reduced time matrix
startIdx_red = startIdx(cols);        % matching start indices

nCols = length(cols);

% --- build row index matrix ---
rowOffsets = (0:blockLen-1)';                 % [389 x 1]
rows = startIdx_red(:)' + rowOffsets;         % [389 x Nt_reduced]

colsMat = repmat(1:nCols, blockLen, 1);       % [389 x Nt_reduced]
linInd = sub2ind(size(A_red), rows, colsMat);
T_ring_fk_v1 = A_red(linInd);                % [389 x Nt_reduced]

%%% matlab FK Version 2
T_Sys_M_fk_v2 = data_fk_v2.Res_System_d18_18_FK(401:end,cols); % FK version 2

A        = data_no_fk.Res_Wasser_d18_18;     % [Nz x Nt]
startIdx = data_no_fk.Res_900_d18_18(7,:);   % [1 x Nt]

A_red        = A(:, cols);            % reduced time matrix
startIdx_red = startIdx(cols);        % matching start indices

nCols = length(cols);

% --- build row index matrix ---
rowOffsets = (0:blockLen-1)';                 % [389 x 1]
rows = startIdx_red(:)' + rowOffsets;         % [389 x Nt_reduced]

colsMat = repmat(1:nCols, blockLen, 1);       % [389 x Nt_reduced]
linInd = sub2ind(size(A_red), rows, colsMat);
T_ring_fk_v2 = A_red(linInd);                % [389 x Nt_reduced]
%%
%%%------------------------------------------%%%
% 04. Spatial alignment MATLAB vs COMSOL
%%%------------------------------------------%%%

dz   = 0.05;     % Spatial resolution COMSOL
dz_M = 0.005;    % Spatial resolution MATLAB

h_M     = (length(T_Sys_M)-1)*dz_M;
h_Sys_C = (length(data_comsol.z)-1)*dz;
h_diff  = abs(h_Sys_C - h_M);

z_Mo = data_comsol.z(1)   - h_diff/2;
z_Mu = data_comsol.z(end) + h_diff/2;
z_M  = flip(z_Mu:dz_M:z_Mo)';


%%%------------------------------------------%%%
% 05. Map Matlab time steps to COMSOL output times
%%%------------------------------------------%%%

% --- Ring gap vertical mapping (replacement volume -> real piston) ---

dz_M = 0.005;          % already used in script
Nz_ring = size(T_ring_fk_v1,1);   % = Nz(5)

h_ring_eq   = Nz_ring * dz_M;     % ~1.945 m (replacement volume height)
h_piston    = 18.0;               % real piston length [m]

z_ring_eq   = (0:Nz_ring-1)' * dz_M;        % 0 ... h_ring_eq
z_ring_real = z_ring_eq * (h_piston / h_ring_eq);  % 0 ... 18 m


%%%------------------------------------------%%%
% 06. Map Matlab time steps to COMSOL output times
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


Nt_plot = length(t_Sys_M);

% --- preallocate ---
T_M        = zeros(size(T_Sys_M,1),        Nt_plot);
T_M_fk_v1  = zeros(size(T_Sys_M_fk_v1,1),  Nt_plot);
T_M_fk_v2  = zeros(size(T_Sys_M_fk_v2,1),  Nt_plot);

T_ring_no_fk_plot  = zeros(size(T_ring_no_fk,1), Nt_plot);
T_ring_fk_v1_plot = zeros(size(T_ring_fk_v1,1), Nt_plot);
T_ring_fk_v2_plot = zeros(size(T_ring_fk_v2,1), Nt_plot);

% --- extract matching time steps (same logic, now extended) ---
for i = 1:Nt_plot
    idx = t_Sys_M(1,i);

    T_M(:,i)        = T_Sys_M(:,idx);
    T_M_fk_v1(:,i)  = T_Sys_M_fk_v1(:,idx);
    T_M_fk_v2(:,i)  = T_Sys_M_fk_v2(:,idx);

    T_ring_no_fk_plot(:,i) = T_ring_no_fk(:,idx);
    T_ring_fk_v1_plot(:,i) = T_ring_fk_v1(:,idx);
    T_ring_fk_v2_plot(:,i) = T_ring_fk_v2(:,idx);
end


%%%------------------------------------------%%%
% 07. Assemble COMSOL system temperature
%%%------------------------------------------%%%

for i = 1:length(t_Sys_M)
    zp(1,i) = 1 + (data_comsol.z(1) - data_comsol.z_p(1,i)) / dz;
end

T_Sys_C = zeros(length(data_comsol.T1_W), length(t_Sys_M));

for i = 1:length(t_Sys_M)
    T_Sys_C(:,i) = data_comsol.T1_W(:,i); % water volume
    T_Sys_C(zp(1,i)+1 : zp(1,i)+length(data_comsol.z_p)-2 , i) = ...
        data_comsol.T1_P(2:end-1,i);     % piston
end
%%
%%%------------------------------------------%%%
% 08. Collect models for plotting (version-safe)
%%%------------------------------------------%%%

models(1).name  = 'COMSOL';
models(1).T     = T_Sys_C;
models(1).z     = data_comsol.z;
models(1).style = '-';
models(1).color = [0.3010 0.7450 0.9330];
models(1).T_gap  = [];          % no ring gap
models(1).z_gap  = [];

models(2).name  = 'MATLAB (no FK)';
models(2).T     = T_M;
models(2).z     = z_M;
models(2).style = '-';
models(2).color = [0.8500 0.3250 0.0980];
models(2).T_gap  = T_ring_no_fk_plot;
models(2).z_gap  = z_ring_real;


models(3).name  = 'MATLAB (FK 16)';
models(3).T     = T_M_fk_v1;
models(3).z     = z_M;
models(3).style = '-';
models(3).color = [0.4660 0.6740 0.1880];
models(3).T_gap  = T_ring_fk_v1_plot;
models(3).z_gap  = z_ring_real;

models(4).name  = 'MATLAB (FK 17)';
models(4).T     = T_M_fk_v2;
models(4).z     = z_M;
models(4).style = '-';
models(4).color = [0.4940 0.1840 0.5560];
models(4).T_gap  = T_ring_fk_v2_plot;
models(4).z_gap  = z_ring_real;

%%%------------------------------------------%%%
% 09. Selected operating points for plotting
%%%------------------------------------------%%%

p = [1, 2, 3, 4, 5];


%%%------------------------------------------%%%
% 10. Static 1D temperature profile plots
%%%------------------------------------------%%%

figure('Name','Temperature profiles','Color','w')

for m = 1:length(p)
    subplot(1,length(p),m)
    hold on

    % Physical piston position (lower edge)
    z_piston = data_comsol.z_p(end,p(m));
    disp(z_piston);
    rectangle('Position', [40, z_piston, 40, 18], ...
              'FaceColor', [1 1 1 0.1], ...
              'EdgeColor', [0 0 0]);

    % Plot all model variants
    for mm = 1:numel(models)
        plot(models(mm).T(:,p(m)), models(mm).z, ...
            models(mm).style, ...
            'Color', models(mm).color, ...
            'LineWidth',1);
        % --- plot ring gap water if available ---
        if ~isempty(models(mm).T_gap)
            plot( ...
                models(mm).T_gap(:,p(m)), ...
                z_piston + models(mm).z_gap, ...
                '--', ...
                'Color', models(mm).color, ...
                'LineWidth',1);
        end
    end

    grid on
    xticks([40 50 60 70 80])
    axis([40 80 z_Mu z_Mo])

    if m == 1
        ylabel('System depth [m]','Interpreter','Latex','Fontsize',12)
        xlabel('Temperature in TPPS [$^{\circ}$C]','Interpreter','Latex','Fontsize',12)
    end

    if m == length(p)
        legend({models.name},'Interpreter','Latex','Fontsize',9)
    end
end

sgtitle(['1D temperature profiles at operating times: ' num2str(t_Sys_M(p)) ], ...
    'Interpreter','Latex','Fontsize',12)

%%
% datestr_run = datestr(now, "yymmdd");
% version = 'v2';
% filename = sprintf('%s_d18_h18_profiles1D_FK_compare_%s.png', datestr_run, version);
% 
% exportgraphics(gcf, fullfile(RESULTS_PROFILES_1D, filename), ...
%     'Resolution', 600);