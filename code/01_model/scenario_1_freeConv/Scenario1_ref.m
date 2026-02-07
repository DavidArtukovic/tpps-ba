%% SZENARIO1.m
% ---------------------------------------------------------------
% PURPOSE:
%   Run a full thermo-hydraulic TPPS simulation for a predefined operating
%   schedule (charging, discharging, idle periods).
%
% DESCRIPTION:
%   - Loads geometry, material data, and initial temperature fields from
%     the initialization .mat file.
%   - Loads the scenario time plan (15-minute control flags) defining the
%     flow direction and operating mode for each timestep.
%   - Preprocesses the control schedule: piston position, discretization
%     in the water column, and charge/discharge flow rates.
%   - For each timestep, calls the coupled heat transfer
%     model (HeattransferSzen) and updates the system state.
%   - Uses a separate helper function to compute energy and exergy
%     balances consistently across the script.
%   - Stores vertical temperature profiles, water temperatures, and
%     energy/exergy-related result vectors for later evaluation and plotting.
%
% ---------------------------------------------------------------

clc
clear

%%%------------------------------------------%%%
% 01. Load Scenario and Initialization Simulation
%%%------------------------------------------%%%

% Load local paths (per-machine config)
run(fullfile('..','..','..', 'configs', 'paths_local.m'));

% relative paths, needed for fk modules
run(fullfile('..','..','..', 'configs','paths_relative.m'));

% Build data subfolder for this configuration
DATA_SCEN1 = fullfile(DATA_BASE, 'scenario1');
DATA_SCEN1_FK = fullfile(DATA_BASE, 'scenario1_freeConv');

% Load init and scenario files
load(fullfile(DATA_SCEN1, 'Init_d18_h18_time8.mat'));   % Geometry, material values and initial values
load(fullfile(DATA_SCEN1, 'SzenarioComsol.mat'));       % Scenario control flags
%%
%%%------------------------------------------%%%
% 02. Initialize Arrays and Geometry
%%%------------------------------------------%%%

% Result arrays
Res_900  = zeros(10, length(t_900)  + 1);   % Values for quarter hours
Res_hour = zeros(10, length(t_hour) + 1);   % Values for full hours

% Temperature vector for the whole system
T_SSys = zeros(Nz(13), 1);

T_V_900 = zeros(Nz(13), length(t_900));     % Temperature array: system grid x time
T_W_900 = zeros(Nz(6),  length(t_900));     % Temperature array: water grid x time

% Spatial discretization
dz = d(4);        % Vertical step size
dr = d(2);        % Radial step size

% Time discretization
dt  = 900;        % Time in seconds for every procedure step
Nt2 = 2;          % Number of timesteps for every procedure step
t2  = (0 : dt/Nt2 : dt);

% Reference temperature for exergy calculation
K = 273.15;

% Cross-section of the water volume [m^2]
Spz = round(SW(1,2) * SW(1,3) * d(1)^2 * pi/4, 10);

% Geometric radii and areas
d_ST   = d(1);           % Diameter in m
r_ST   = d_ST / 2;       % Radius in m
r_gap  = 0.5;            % Annular gap thickness in m
r_pist = r_ST - r_gap;   % Piston radius in m

A(1) = d(1)^2 * pi/4;    % Area of the cylindrical storage
A(2) = r_pist^2 * pi;    % Area of the piston
A(3) = A(1) - A(2);      % Area of the annulus (ring gap)

% Initialize positions indices of bypass inlets and mebrane
[idx_bypass_lower_vec, idx_bypass_upper_vec, idx_membrane_vec] = compute_bypass_indices(t_900(2,:), H(3), H(11), Nz, dz);
bypass_indices = [idx_bypass_lower_vec; idx_bypass_upper_vec; idx_membrane_vec];

%%
%%%------------------------------------------%%%
% 03. Preprocess Scenario: piston position and flow rates
%%%------------------------------------------%%%

[Res_900, LC, DC, Lflow, Dflow] = prepare_scenario_flow(Res_900, t_900, Nz, H, d);
%%
%%%------------------------------------------%%%
% 04. Set Initial Temperature Fields
%%%------------------------------------------%%%

% Set temperature vector to initial state
T_Sys = IC_Sys;

DTRE = zeros(1, Nz(12));   % Vector of radial temperature differences in radial soil

% Temperature in the whole vertical system
T_V(1, :) = T_Sys(1, 1:Nz(10));

% Temperature in water system
T_W(1, :) = T_Sys(1, Nz(3)+Nz(8):Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3);

% Exergy density of the water (initialized, will be overwritten)
wEX = zeros(1, length(T_W));

%%
%%%------------------------------------------%%%
% 05. Set Initial Energies
%%%------------------------------------------%%%

% Compute initial energy and exergy balances using helper function
[Heat_insu, Heat_Wasser, Heat_piston, Heat_Vsoil, ...
    Heat_Rsoil, Heat_Rinsu, WEX, wEX, DTRE] = ...
    compute_energy_balances(T_V, T_W, T_Sys, T_REf, ...
                            Nz, dz, A, SW, z_RE, z_W, K, Spz);

% Store initial values in hourly result array
Res_hour(1,1)  = Heat_insu;
Res_hour(2,1)  = Heat_Wasser;
Res_hour(3,1)  = Heat_piston;
Res_hour(4,1)  = Heat_Vsoil;
Res_hour(5,1)  = Heat_Rsoil;
Res_hour(6,1)  = Heat_Rinsu;
Res_hour(11,1) = WEX;

% Store initial values in 15-minute result array
Res_900(1,1)  = Heat_insu;
Res_900(2,1)  = Heat_Wasser;
Res_900(3,1)  = Heat_piston;
Res_900(4,1)  = Heat_Vsoil;
Res_900(5,1)  = Heat_Rsoil;
Res_900(6,1)  = Heat_Rinsu;
Res_900(11,1) = WEX;

% Initial temperatures in °C (for boundary conditions)
T0init(1) = 45;                 % Initial storage temperature
T0init(2) = 80;                 % Supply temperature
T0init(3) = 45;                 % Return temperature
T0init(4) = 20;                 % Air temperature
T0init(5) = 11;                 % Ground temperature
T0init(6) = T(3,1) - 273.15;    % Boundary temperature for initial radial soil
T0init(7) = 11;                 % Initial insulation temperature

%%
%%%------------------------------------------%%%
% 06. Time Loop over Scenario (t_900/2 for 10 days)
%%%------------------------------------------%%%


% Initialize logging
version = 'v18';
dateTag = datestr(now, "yymmdd");

FK_LOG_BASE = fullfile(DATA_SCEN1_FK, '01_logs');
NOFK_LOG_BASE = fullfile(DATA_SCEN1, '01_logs');

fk_logfile = fullfile(FK_LOG_BASE, ...
    [dateTag '_scenario1_FK_' version '.log']);

fid_fk = fopen(fk_logfile,'w');

fklog = @(s) fprintf(fid_fk,'%s\n',string(s));

fklog(sprintf([ ...
    '=== Scenario 1 Free Convection Simulation ===\n' ...
    'Description:\n' ...
    'FK with shifting membrane adapted convection diffusion  \n'...
    'Added upwind scheme in water \n'...

]));
fklog(sprintf('Date and Time: %s', dateTag));
fklog(sprintf('Version: %s', version));

% create fk code history vector
fk_code_hist = zeros(1, length(t_900));

for i = 1:(length(t_900))
    tic
    fklog('|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||');
    fklog(sprintf('t = %.2f h', i * 15 / 60));
    fklog(sprintf('Temperature in middle of ring gap:  %.2f °C', T_Sys(end, Nz(3)+Nz(8)+Nz(2)+round(Nz(5)/2-3))));


    % Update vertical water discretization according to piston position
    Nz(2) = Res_900(7, i+1);   % Number of cells in lower pressure zone
    Nz(4) = Res_900(8, i+1);   % Number of cells in upper pressure zone

    % Flow velocity in the storage system [m/s] (charging or discharging)
    flow = Res_900(9, i+1) + Res_900(10, i+1);
    disp(['flow        = ' num2str(flow)]);

    % Compute heat transfer and mass transport for the current 15-min step
    [T_Sys, T_REf, fk_code] = HeattransferSzen(t2, IC_Sys, Nz, dz, bypass_indices(:,i),...
                                               flow, T0init, SW, A, z_RE, T_REf, Nt2,fklog);

    % Update initial condition for the next procedure step
    IC_Sys = T_Sys(end, :);

    % Updated temperatures in the vertical system
    T_V(1, :) = T_Sys(end, 1:Nz(10));

    % Updated water temperature profile
    T_W(1, :) = T_Sys(end, Nz(3)+Nz(8):Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3);
    disp(['temperature below membrane: ' num2str(T_Sys(end, bypass_indices(3,i)-3)) ' °C']);
    
    % Recompute energy and exergy balances using helper function
    [Heat_insu, Heat_Wasser, Heat_piston, Heat_Vsoil, ...
        Heat_Rsoil, Heat_Rinsu, WEX, wEX, DTRE] = ...
        compute_energy_balances(T_V, T_W, T_Sys, T_REf, ...
                                Nz, dz, A, SW, z_RE, z_W, K, Spz);

    % Store quarter-hourly results
    Res_900(1, 1+i)  = Heat_insu;
    Res_900(2, 1+i)  = Heat_Wasser;
    Res_900(3, 1+i)  = Heat_piston;
    Res_900(4, 1+i)  = Heat_Vsoil;
    Res_900(5, 1+i)  = Heat_Rsoil;
    Res_900(6, 1+i)  = Heat_Rinsu;
    Res_900(11, 1+i) = WEX;          % Exergy in water volume

    % Store instantaneous water temperature (flipped for plotting)
    T_W_900(:, i) = flip(T_W);

    % Assemble system temperature vector for plotting (piston, zones, insulation)
    T_SSys(1:Nz(2)) = T_V(Nz(3)+Nz(8):Nz(3)+Nz(8)+Nz(2)-1);   % Lower pressure zone

    T_SSys(Nz(2)+Nz(3)-1 : Nz(2)+Nz(3)+Nz(4)-2) = ...
        T_V(Nz(3)+Nz(8)+Nz(2)+Nz(5)-2 : Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3); % Upper pressure zone

    T_SSys(Nz(2) : Nz(2)+Nz(3)-1) = T_V(1:Nz(3));                               % Piston
    T_SSys(Nz(2)+Nz(3)+Nz(4)-2:end) = ...
        T_V(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3:end);                               % Vertical insulation

    T_V_900(:, i) = flip(T_SSys);

    % Store hourly aggregated values every 4th 15-min step
    if mod(i, 4) == 0
        idx_hr = 1 + i/4;
        Res_hour(1, idx_hr)  = Heat_insu;
        Res_hour(2, idx_hr)  = Heat_Wasser;
        Res_hour(3, idx_hr)  = Heat_piston;
        Res_hour(4, idx_hr)  = Heat_Vsoil;
        Res_hour(5, idx_hr)  = Heat_Rsoil;
        Res_hour(6, idx_hr)  = Heat_Rinsu;
        Res_hour(7, idx_hr)  = Res_900(7, 1+i);   % Number of cells lower water volume
        Res_hour(8, idx_hr)  = Res_900(8, 1+i);   % Number of cells upper water volume
        Res_hour(9, idx_hr)  = Res_900(9, 1+i);   % Charging velocity [m/s]
        Res_hour(10, idx_hr) = Res_900(10,1+i);   % Discharging velocity [m/s]
        Res_hour(11,idx_hr)  = Res_900(11,1+i);   % Exergy in the water volume
    end
    fk_code_hist(i) = fk_code;
    disp(i * 15 / 60);   % Elapsed simulation time in hours
    fk_code_hist(i) = fk_code;
    toc
end

fclose(fid_fk);
%%
%%%------------------------------------------%%%c
% 07. Store Scenario Results
%%%------------------------------------------%%%
fk = true;

ResOut = struct();

ResOut.meta.date        = dateTag;
ResOut.meta.diameter    = d_ST;
ResOut.meta.height      = H(3);
ResOut.meta.version     = version;
ResOut.meta.use_fk      = fk;

ResOut.res.series_900  = Res_900;
ResOut.res.series_hour = Res_hour;

ResOut.temperature.water  = T_W_900;
ResOut.temperature.system = T_V_900;

ResOut.geometry.bypass_indices = bypass_indices;


if fk
    ResOut.fk.logging_code = fk_code_hist;
    modeTag   = 'FK';
    dataPath  = DATA_SCEN1_FK;
else
    ResOut.fk.logging_code = [];
    modeTag   = 'noFK';
    dataPath  = DATA_SCEN1;
end

filenameSIM = sprintf('%s_d%d_h%d_Res_Matlab_%s_%s.mat', dateTag, d_ST, H(3), modeTag, version);
fullpathSIM = fullfile(dataPath, filenameSIM);
save(fullpathSIM, "ResOut");


%%
%%% ============================================================ %%%
%                           LOCAL FUNCTIONS                       
%%% ============================================================ %%%

function [Res_900, LC, DC, Lflow, Dflow] = prepare_scenario_flow(Res_900, t_900, Nz, H, d)
    % PREPARE_SCENARIO_FLOW
    %   Preprocess the scenario schedule:
    %   - Compute piston position and number of cells in upper/lower pressure zones.
    %   - Count total charge and discharge hours.
    %   - Compute average charge/discharge flow rates.
    %   - Compute per-step charge/discharge velocities in the storage.
    %
    % INPUT:
    %   Res_900 - result array for 15-minute values (modified in-place)
    %   t_900   - scenario matrix with flow control flags (2 x Nt)
    %   Nz      - vector of grid sizes in the system
    %   H, d    - parameter vectors (geometry and operating parameters)
    %
    % OUTPUT:
    %   Res_900 - updated with piston index, zone sizes and flow velocities
    %   LC      - total charging hours
    %   DC      - total discharging hours
    %   Lflow   - average charge flow [m/s]
    %   Dflow   - average discharge flow [m/s]

    LC = 0;
    DC = 0;

    % Compute piston position and zone heights for each 15-min step
    for k = 1:length(t_900)
        % Number of vertical cells in the lower pressure zone
        % Lower dead-zone + piston position (t_900=1 implies at the top)
        Res_900(7, k+1) = round(1 + (H(7) + t_900(2,k) * H(11)) / d(4));
        % Number of vertical cells in upper pressure zone
        Res_900(8, k+1) = Nz(1) - Res_900(7, k+1) - Nz(3) + 2;

        % Count charging/discharging quarter-hours
        if t_900(1,k) > 0
            LC = LC + 1;
        elseif t_900(1,k) < 0
            DC = DC + 1;
        end
    end

    % Convert to full charge/discharge hours
    LC = LC / 4;
    DC = DC / 4;

    % Average flow for a full charge or discharge cycle [m/s]
    Lflow = H(6) / (LC * 3600);   % Average charge flow
    Dflow = H(6) / (DC * 3600);   % Average discharge flow

    % Per-step flow velocities for each 15-minute control step
    for k = 1:length(t_900)
        if t_900(1,k) > 0
            Res_900(9, k+1) = -Lflow * t_900(1,k);   % Charging velocity, negative value
        elseif t_900(1,k) < 0
            Res_900(10, k+1) = -Dflow * t_900(1,k);  % Discharging velocity, positive value
        end
    end
end

function [Heat_insu, Heat_Wasser, Heat_piston, Heat_Vsoil, ...
          Heat_Rsoil, Heat_Rinsu, WEX, wEX, DTRE] = ...
          compute_energy_balances(T_V, T_W, T_Sys, T_REf, ...
                                  Nz, dz, A, SW, z_RE, z_W, K, Spz)
    % COMPUTE_ENERGY_BALANCES
    %   Compute all thermal energy contents and the exergy of the water volume
    %   based on the current temperature fields in the TPPS model.
    %
    % INPUT:
    %   T_V    - vertical temperature vector of the whole system
    %   T_W    - water temperature vector
    %   T_Sys  - full 1D system state vector (including flatten radial soil)
    %   T_REf  - radial soil/insulation temperature matrix (radial x vertical)
    %   Nz     - vector with numbers of grid points in each subdomain
    %   dz     - vertical grid spacing
    %   A      - areas (storage, piston, annulus)
    %   SW     - material properties (density * heat capacity)
    %   z_RE   - radial grid information for soil/insulation
    %   z_W    - vertical grid of the water volume
    %   K      - reference temperature for exergy computation
    %   Spz    - water cross-section area
    %
    % OUTPUT:
    %   Heat_insu   - thermal energy in vertical insulation [J]
    %   Heat_Wasser - thermal energy in water column [J]
    %   Heat_piston - thermal energy in the piston [J]
    %   Heat_Vsoil  - thermal energy in vertical soil [J]
    %   Heat_Rsoil  - thermal energy in radial soil (before insulation) [J]
    %   Heat_Rinsu  - thermal energy in radial insulation [J]
    %   WEX         - total exergy in the water volume [J]
    %   wEX         - exergy density along the water column
    %   DTRE        - auxiliary vector with radial soil contributions

    % --- Vertical insulation energy [J] ---
    Heat_insu = SW(3,2) * SW(3,3) * dz * A(1) * ...
        (sum(T_V(1, Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2:end-1)) + ...
         0.5 * (T_V(1, Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3) + T_V(1,end)));

    % --- Water column energy [J] ---
    Heat_Wasser = SW(1,2) * SW(1,3) * dz * A(1) * ...
        (sum(T_W(1,2:end-1)) + 0.5 * (T_W(1,1) + T_W(1,end)));

    % --- Piston energy [J] ---
    % Note: The sum uses the first row of T_Sys (consistent with the
    % original loop implementation), while the boundary uses the last row.
    Heat_piston = SW(2,2) * SW(2,3) * dz * A(2) * ...
        (sum(T_Sys(1, 2:Nz(3)-1)) + ...
         0.5 * (T_Sys(end,1) + T_Sys(end, Nz(3))));

    % --- Vertical soil energy below the system [J] ---
    Heat_Vsoil = SW(2,2) * SW(2,3) * dz * A(1) * ...
        (sum(T_V(1, Nz(3)+2:Nz(3)+Nz(8)-1)) + ...
         0.5 * (T_V(1, Nz(3)+1) + T_V(1, Nz(3)+Nz(8))));

    % --- Radial soil energy (before insulation) [J] ---
    DTRE = zeros(1, Nz(12));   % Will be overwritten each call
    for m = 1:Nz(12)
        DTRE(m) = dz * z_RE(2,m) * ...
            (sum(T_REf(2:end-Nz(9), m)) + ...
             0.5 * (T_REf(1,m) + T_REf(end-Nz(9)+1, m)));   % 'K*m^3'
    end
    Heat_Rsoil = SW(2,2) * SW(2,3) * sum(DTRE);

    % --- Radial insulation energy [J] ---
    for m = 1:Nz(12)
        DTRE(m) = dz * z_RE(2,m) * ...
            (sum(T_REf(end-Nz(9)+2:end-1, m)) + ...
             0.5 * (T_REf(end-Nz(9)+1, m) + T_REf(end, m))); % 'K*m^3'
    end
    Heat_Rinsu = SW(3,2) * SW(3,3) * sum(DTRE);

    % --- Exergy in the water volume [J] ---
    wEX = zeros(1, length(z_W));
    for p = 1:length(z_W)
        wEX(1,p) = T_W(p) - K * log((K + T_W(p)) / K);
    end
    WEX = Spz * dz * ...
        (sum(wEX(1,2:end-1)) + 0.5 * (wEX(1,1) + wEX(1,end)));
end

function [idx_bypass_lower_vec, idx_bypass_upper_vec, idx_membrane_vec] = compute_bypass_indices(soc_vec, h_pist, h_lift, Nz, dz)
    % COMPUTE_BYPASS_INDICES
    % --------------------------------------------------------------
    % Computes lower and upper bypass inlet indices and mebrane index
    % based on piston position, following the piecewise definition from the sketch.
    %
    % INPUT:
    %   soc_vec          - normalized vector of piston position in [0,1]
    %                     (soc=1: piston top, soc=0: piston bottom)
    %   h_pist          - total piston stroke height [m]
    %   h_lift          - height of maximum lift [m]
    %   Nz              - grid size vector
    %   dz              - vertical grid spacing [m]
    %
    % OUTPUT:
    %   idx_bypass_lower - global index of lower bypass inlet
    %   idx_bypass_upper - global index of upper bypass inlet
    %   idx_mebrane_vec  - global index of mebrane node
    %
    % NOTES:
    %   - Lower bypass always remains inside replacement volume.
    %   - Upper bypass is inside replacement volume for large s and
    %     transitions into upper water volume once a threshold is crossed.
    %   - No dead-zone correction applied (by design).


    % Number of cells in replacement volume of ringgap
    Nrep = Nz(5);
    idx_rep_begin = Nz(3) + Nz(8) + Nz(2) - 1; % Start of replacement volume
    idx_upper_begin = Nz(3) + Nz(8) + Nz(2) + Nz(5) -2 ; % Start of upper water volume

    % ----------------------------------------------------------
    % Lower bypass: always inside replacement volume
    % ----------------------------------------------------------
    idx_bypass_lower_vec = idx_rep_begin + round((1 - soc_vec) .* Nrep .* h_lift ./ h_pist);

    % ----------------------------------------------------------
    % Upper bypass: piecewise definition
    % ----------------------------------------------------------

    % Threshold for leaving replacement volume
    % upper bypass ist 8m below upper dead zone. Thus as soon as the piston drives
    % down more than 8m, the bypass enters the upper water volume
    % This is given by s_crit = 1 - 8/h_lift
    s_crit = 1.0 - 8/h_lift;

    mask_rep = soc_vec >= s_crit;      % logical 1xN
    mask_up  = ~mask_rep;

    % Case 1: upper bypass still in replacement volume
    % upper bypass is always 10m above lower bypass, thus +10/18 of Nrep
    idx_bypass_upper_vec(mask_rep) = idx_bypass_lower_vec(mask_rep) + round(Nrep *10/18); 

    % Case 2: upper bypass in upper water volume
    delta_s = s_crit - soc_vec(mask_up);
    idx_bypass_upper_vec(mask_up) = idx_upper_begin + round(delta_s .* h_lift ./ dz);

    % ----------------------------------------------------------
    % Membrane index
    % ----------------------------------------------------------

    idx_membrane_vec = idx_rep_begin ...
    + round( ...
        Nrep/2 ...                          % mid at soc=1
      + (1 - soc_vec) .* (17/h_pist) .* (Nrep/2) ...
    );
end
