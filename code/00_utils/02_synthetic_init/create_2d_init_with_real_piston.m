 % Inspect old

clc
clear

%%%------------------------------------------%%%
% 01. Load Scenario and Initialization Simulation
%%%------------------------------------------%%%

% Load local paths (per-machine config)
run(fullfile('..','..','..', 'configs', 'paths_local.m'));

% Build data subfolder for this configuration
DATA_INIT = fullfile(DATA_BASE, 'init');
%%
% Load init and scenario files
InitOut_snythetic = load(fullfile(DATA_INIT, '20260301_d18_hp16.0_gap0.5_2D_chunk008_v2_synth.mat'));   % Geometry, material values and initial values
InitOut_piston = load(fullfile(DATA_INIT, '20260330_Init_piston_d18_hp18.0_gap0.5_Init2D_chunk009_v1.mat'));       % Scenario control flags
%%

% Extract data
T_RPf = InitOut_piston.InitOut.state.T_RPf_end;   % [Nz(3) x Nz(14)]
z_RP  = InitOut_piston.InitOut.grid.z_RP;         % radial grid
dz    = InitOut_piston.InitOut.grid.dz;

% Build axes
z = (0:size(T_RPf,1)-1) * dz;      % axial coordinate
r = z_RP(1,:);                     % radial coordinate

% Plot
figure
imagesc(r, z, T_RPf)
set(gca, 'YDir','normal')
colorbar
xlabel('r [m]')
ylabel('z [m]')
title('2D Piston Temperature Field')

%%

InitOut = InitOut_snythetic.InitOut;

%% Overwrite with real piston 2D data
T_RPf_real = InitOut_piston.InitOut.state.T_RPf_end;   % [Nz(3) x Nz(14)]

% Store real 2D piston field explicitly
InitOut.state.T_RPf_end = T_RPf_real;

n_piston_2d = numel(T_RPf_real);
InitOut.state.IC_Sys_end(end-n_piston_2d+1:end) = T_RPf_real(:);

save(fullfile(DATA_INIT, '20260401_d18_hp18.0_gap0.5_2D_chunk009.mat'), 'InitOut', '-v7.3')
