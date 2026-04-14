 % Inspect old

clc
clear

%%%------------------------------------------%%%
% 01. Load Scenario and Initialization Simulation
%%%------------------------------------------%%%

% Load local paths (per-machine config)
run(fullfile('..','..', 'configs', 'paths_local.m'));

% Build data subfolder for this configuration
DATA_SCEN1 = fullfile(DATA_BASE, 'scenario1');

% Load init and scenario files
load(fullfile(DATA_SCEN1, 'Init_d18_h18_time8.mat'));   % Geometry, material values and initial values
load(fullfile(DATA_SCEN1, 'SzenarioComsol.mat'));       % Scenario control flags

% Build data subfolder for this configuration
DATA_INIT = fullfile(DATA_BASE, 'init');

%% ------------------------------------------------------------
%  02. Minimal 2D visualization of radial soil temperature field
% -------------------------------------------------------------

dz = d(2);                                 % axial grid spacing [m]
Nr = Nz(12);                                % number of radial nodes (soil)
Nz_ax = Nz(13);                             % number of axial nodes (soil)

sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4;
% Extract soil temperatures from IC_Sys
idx_start = sys_top_idx + 1;
idx_end   = idx_start + Nr*Nz_ax - 1;

T_RE_vec = IC_Sys(idx_start:idx_end);

% Reshape into 2D field (axial x radial)
T_RE = reshape(T_RE_vec, Nz_ax, Nr);

% Axial coordinate
z_ax = (0:Nz_ax-1) * dz;

% Radial coordinates (non-equidistant)
r = z_RE(1,:);

% Create 2D plot
figure
pcolor(r, z_ax, T_RE)
shading flat
colorbar

xlabel('r [m]')
ylabel('z [m]')
title('Radial soil temperature field (IC)')

%% ------------------------------------------------------------
%  03. 1D piston temperature profile (IC)
% -------------------------------------------------------------

dz = 0.005;                     % axial grid spacing [m]
Nz_pist = Nz(3);                % number of piston axial nodes

% Extract piston temperatures
T_pist = IC_Sys(1:Nz_pist);

% Axial coordinate
z_pist = (0:Nz_pist-1) * dz;

% Plot
figure
plot(T_pist, z_pist, 'LineWidth', 2)
set(gca, 'YDir', 'reverse')     % optional: if z=0 is top
grid on

xlabel('T [°C]')
ylabel('z [m]')
title('1D piston temperature profile (IC)')

%% ------------------------------------------------------------
%  04. Create 2D piston geometry
% -------------------------------------------------------------

n_radial_piston_cells = 10;                     % 10 cells
n_radial_piston_nodes = n_radial_piston_cells + 1;

Nz(14) = n_radial_piston_nodes+1;                 % store node count pölus interface cell

% geometric progression (doubling inward)
dr0_piston = 8.5 / (2^n_radial_piston_nodes - 1);

dr_vec = dr0_piston * 2.^(0:n_radial_piston_nodes-1);   % cell widths

% radial node positions (0 ... R)
r_nodes = [0, cumsum(flip(dr_vec))];                  % length = 11, first = 0, last = r_pist

% --- store like soil grid
z_RP = zeros(4, Nz(14));
z_RP(1,:) = r_nodes;                            % nodes radial positions

% Area of first circular ring segment [m²]
z_RP(2,1) = round(pi*(z_RP(1,2)-z_RP(1,1))*((z_RP(1,2)-z_RP(1,1))/4 + z_RP(1,1)), 10); 

% Area of last circular ring segment [m²]                                            
z_RP(2,end) = round(-pi*(z_RP(1,end-1)-z_RP(1,end)) * ((z_RP(1,end-1)-z_RP(1,end))/4 + z_RP(1,end)), 10);
                                             
% Area of intermediate ring segments [m²]
for i = 2:Nz(14)-1
    z_RP(2,i) = round(pi*(z_RP(1,i+1)-z_RP(1,i-1)) * ...
                     (z_RP(1,i-1) + 0.5*((z_RP(1,i+1)-z_RP(1,i-1))/2 ...
                      + z_RP(1,i) - z_RP(1,i-1))) ...
                      ,10);                                     
end

for i = 2:Nz(14)-1
    z_RP(3,i) = round(z_RP(1,i+1) - z_RP(1,i-1), 10);  % symmetric spacing to next/previous node
    z_RP(4,i) = round(z_RP(1,i+1) - z_RP(1,i),   10); % forward spacing to next node
    z_RP(5,i) = round(z_RP(1,i)   - z_RP(1,i-1), 10); % backward spacing to previous node
end

%% ------------------------------------------------------------
%  05. Synthetic 2D piston temperature field
% -------------------------------------------------------------

R      = 8.5;                 % piston radius
RE = 9;
T_wall = 45;                  % boundary temperature [°C]

Nz_ax  = Nz(3);               % axial piston nodes
Nr     = Nz(14);              % radial nodes (including r=0 and R)

dz     = 0.005;

% --- 1D reference profile
T_1D = IC_Sys(1:Nz_ax);

% --- radial coordinates (non-equidistant)
r = z_RP(1,:);                % size Nr

% --- radial shape parameters
alpha = 0.4;                  % camel intensity (0...1)
p     = 2;                    % edge decay exponent

% --- build radial shape function φ(r)
phi = (1-alpha)*(1 - (r/R).^p) ...
      + alpha*(4*(r/R).*(1 - r/R));

% enforce exact boundary conditions numerically
phi(1)   = 1;                 % r = 0
phi(end) = 0;                 % r = R

% --- build 2D temperature field
T_2D = zeros(Nz_ax, Nr);

for i = 1:Nz_ax
    A = T_1D(i) - T_wall;
    T_2D(i,:) = T_wall + A * phi;
end

%% ------------------------------------------------------------
%  06. Energy matching
% -------------------------------------------------------------

% --- radial cell widths (from nodes)
dr = diff(r);                             % length Nr-1
r_mid = 0.5*(r(1:end-1)+r(2:end));        % cell centers

% --- 1D energy (per unit density*cp)
E_1D = sum(T_1D) * dz * (pi*R^2);

% --- 2D energy
E_2D = 0;

for i = 1:Nz_ax
    for j = 1:Nr-1
        dV = 2*pi*r_mid(j) * dr(j) * dz;
        E_2D = E_2D + T_2D(i,j) * dV;
    end
end

% --- scale amplitude relative to wall temperature
scale = E_1D / E_2D;

T_2D = T_wall + scale*(T_2D - T_wall);

%% ------------------------------------------------------------
%  07. Quick figure for Check
% -------------------------------------------------------------
figure
pcolor(r, (0:Nz_ax-1)*dz, T_2D)
shading flat
colorbar
xlabel('r [m]')
ylabel('z [m]')
title('Synthetic 2D piston IC')

%% ------------------------------------------------------------
%  08. Add material parameters
% -------------------------------------------------------------


SW(4,4) = (2*SW(2,1)) / ((RE^2-R^2) * SW(1,2) * SW(1,3) * log((R+0.005)/R));
                                             % Radial loss factor: ring-gap water <-> soil, referenced to ring-gap water domain [1/s]

SW(4,5) = (2*SW(2,1)) / ((RE^2-R^2) * SW(1,2) * SW(1,3) * log(R/(R - dr0_piston)));
                                             % Radial loss factor: ring-gap water <-> piston, referenced to ring-gap water domain [1/s]



%% ------------------------------------------------------------
%  09. Save synthetic data based on structure from init.m
% -------------------------------------------------------------

T_RPf_vec = T_2D(:).';
% ------------------------------------------------------------
%  Append 2D piston to IC_Sys
% -------------------------------------------------------------

IC_Sys_old = IC_Sys(1:idx_end);                      % backup

IC_Sys = [IC_Sys_old , T_RPf_vec];


disp(['Length IC_Sys = ', num2str(length(IC_Sys))])
disp(['Expected length = ', num2str( ...
    Nz(10) + Nz(13)*Nz(12) + Nz(3)*Nz(14))]);


% ------------------------------------------------------------
%  Build InitOut structure
% -------------------------------------------------------------
version = 'v2_synth';
modeTag = '2D';
dateTag = datestr(now,'yyyymmdd');
o = 8;

InitOut = struct();

% --- Meta ---
InitOut.meta.dateTag  = dateTag;
InitOut.meta.modeTag  = modeTag;
InitOut.meta.version  = version;
InitOut.meta.chunk_o  = o;
InitOut.meta.diameter = d(1);
InitOut.meta.h_pist   = 16;
InitOut.meta.r_gap    = 0.5;

% --- Grid ---
InitOut.grid.Nz   = Nz;
InitOut.grid.dz   = dz;
InitOut.grid.z_RE = z_RE;
InitOut.grid.z_RP = z_RP;

% --- Geometry ---
InitOut.geom.H      = H;
InitOut.geom.d      = d;
InitOut.geom.z_OD   = z_OD;
InitOut.geom.z_UD   = z_UD;
InitOut.geom.z_W    = z_W;
InitOut.geom.z_Pist = z_Pist;
InitOut.geom.z_SSys = z_SSys;
InitOut.geom.z_Sys  = z_Sys;

% --- Material ---
InitOut.param.SW     = SW;
InitOut.param.T0init = [];

% --- State ---
InitOut.state.IC_Sys_end = IC_Sys;
InitOut.state.T_V_end    = T_V;
InitOut.state.T_W_end    = [];
InitOut.state.T_RPf_end  = T_2D;
InitOut.state.T_REf_end  = T_RE;

% --- Save ---
filenameSIM = sprintf('%s_d%d_hp%.1f_gap%.1f_%s_chunk%03d_%s.mat', ...
                      dateTag, d(1), 16, 0.5, modeTag, o, version);

fullpathSIM = fullfile(DATA_INIT, filenameSIM);
save(fullpathSIM, "InitOut");