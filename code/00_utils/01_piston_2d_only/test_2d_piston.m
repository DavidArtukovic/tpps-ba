% Build full file path
filename = '20260319_Init_piston_d18_hp18.0_gap0.5_Init2D_chunk001.mat';
fullpath = fullfile(DATA_INIT, filename);

% Load struct
S = load(fullpath);
InitOut = S.InitOut;

% Extract data
T_RPf = InitOut.state.T_RPf_end;   % [Nz(3) x Nz(14)]
z_RP  = InitOut.grid.z_RP;         % radial grid
dz    = InitOut.grid.dz;

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