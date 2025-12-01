clc 
clear


%%%------------------------------------------%%%
% 01 Initialize Time and Temperature
%%%------------------------------------------%%%

t0 = 0;
t1 = 2890*3600;                         % 120 days in s
t_end = 35064*3600+3600*1700;           % 1532 days in s or 4.2 years

time_short = t0:900:t1;                 % short time array in 15 minutes steps
time_long = t0:900:t_end;               % long time array in 15 minutes steps

T = zeros(3,length(time_long));         % 3-row temperature array (only third row relevant)

% Custom linear temperature profile for 120 days
for k = 1:length(time_short)
    T(3,k) = 11+273.15 + time_long(k)*(286.25-(11+273.15))/(t1);
end

% Custom temperature profile for 4.2 years, as function of sin and exponential terms.
for k = length(time_short):length(time_long)
    T(1,k) = 11+79*(1-exp(-time_long(k)/(35064*3600/3))); % an 1
    T(2,k) = 11+44*(1-exp(-time_long(k)/(35064*3600/3))); % an 1
    T(3,k) = 264.75+((T(1,k)-T(2,k))/2)*sin(time_long(k)*(2*pi)/(3600*8766)+3600*8764/4)+((T(1,k)-T(2,k))/2)+T(2,k);
end


%%%------------------------------------------%%%
% 02 Storage Basic Geometry
%%%------------------------------------------%%%

d_ST = 18;                              % diameter in m
r_ST = d_ST/2;                          % radius in m
h_tot = 1;                              % dead zone in m
h_lift = 16;                            % lifting height in m
h_Water = h_lift + 2*h_tot;             % cylindric water column in m
V_WZ = h_Water*pi/4*d_ST^2;             % volume cylindirc water columns in m³

%Prallplatte Einlass/Auslass
d_in = 3.6; % Durchmesser Prallplatte in m
r_in = d_in/2; % Radius Prallplatte in m
h_in = 0.05; % Höhe Prallplatte in Meter
V_PP = h_in*d_in^2/4*pi; % Volumen Prallplatte in m³

% Geometry of dome
h_kupp = 1;                             % height dome in m
V_Kupp = (2/3)*r_ST^2*h_kupp*pi;        % volume half dome in m³

% Piston geometry
h_pist = 18;                            % piston heigth in m
r_gap = 0.5;                            % ring gap (Ringspalt) in m
r_pist = r_ST-r_gap;                    % piston radius in m
V_pist = h_pist*pi*r_pist^2;            % volume piston in m³
V_gap = h_pist*pi*(r_ST^2-r_pist^2);    % volume ring gap (Ringspalt) in m³;

% Surrounding area
h_insu = 2;                             % isolation height in m
r_soil = 9;                             % radial soil in m 
h_soil = 1+r_soil;                      % height vertical soil in m

% Membrane
h_mem = 0.05;                           % height of membrane in m
V_mem = h_mem*pi*(r_ST^2-r_pist^2);     % volume of membrane in m

% Pressure tunnel (Druckstollen)
r_out = 0.05;                                   % radius circular ring (Kreisring) in m
r_out_gap = 0.01;                               % radius outlet negativ in m;
h_out = 0.045;                                  % height outlet
V_To = pi^2*r_ST*r_out^2;                       % volume circular ring (Kreisring) in m³ (Approx)
V_Out = h_out*pi*((r_ST+r_out_gap)^2-r_ST^2);   % volume outlet negativ in m³

% Volume water
V_W_gap = 2*(V_To-V_Out)-V_mem+V_gap;           % water volume in ring gap (Ringspalt) in m³
V_W_ST = 2*(V_Kupp-V_PP)+V_WZ;                  % water volume in cylinder and dome in m³
V_W = V_W_gap+V_W_ST;                           % Total water volume in m³


%%%------------------------------------------%%%
% 03 Energy Stored
%%%------------------------------------------%%%

% Potential energy 
E_pot = V_pist*(2400-1000)*9.81*h_lift/(3600*1e3);        % m*g*h [kWh] with m_eff​=Vpist​(ρ_piston​−ρ_water​)

% Thermal energy
E_HWS = V_W*1000*4200*(80-11)/(3600*1e3);                 % V_w*ρ*cp*Dt [kWh]


%%%------------------------------------------%%%
% 04 Equivalents
%%%------------------------------------------%%%

h_1D_zyl = V_W_ST*4/(d_ST^2*pi);                % 1D height equivalent of the cylindrical water section [m]
h_1D_ps  = V_W_gap*4/(d_ST^2*pi);               % 1D height equivalent of the pressure shaft section [m]
h_1D     = h_1D_zyl + h_1D_ps;                  % Total 1D height equivalent of the system for identical cross-section [m]

d_ps_1D  = (V_W_gap*4/(pi*h_pist))^0.5;         % 1D equivalent diameter of the pressure shaft [m]

h_1D_tot = h_tot + (V_Kupp - V_PP)/(d_ST^2*pi/4); % Total 1D height incl. cupola & piston chamber volume [m]


%%%------------------------------------------%%%
% 05 Substance Values
%%%------------------------------------------%%%

% Water
SW(1,1) = 0.6;                                  % Thermal conductivity λ of water [W/mK]
SW(1,2) = 1000;                                 % Density ρ of water [kg/m³]
SW(1,3) = 4200;                                 % Specific heat capacity cp of water [J/kgK]
SW(1,4) = round(sqrt(SW(1,1)*SW(1,2)*SW(1,3)),10); % Thermal diffusivity coefficient √(λρcp) [J/m²K·s^0.5]
SW(1,5) = SW(1,1)/(SW(1,2)*SW(1,3));            % Thermal diffusivity α = λ/(ρcp) [m²/s]

% Soil mass (approximated as concrete-like)
SW(2,1) = 2.1;                                  % Thermal conductivity λ of soil mass [W/mK]
SW(2,2) = 2400;                                 % Density ρ of soil mass [kg/m³]
SW(2,3) = 1000;                                 % Specific heat capacity cp of soil mass [J/kgK]
SW(2,4) = round(sqrt(SW(2,1)*SW(2,2)*SW(2,3)),10); % Thermal diffusivity coefficient √(λρcp) [J/m²K·s^0.5]
SW(2,5) = SW(2,1)/(SW(2,2)*SW(2,3));            % Thermal diffusivity α = λ/(ρcp) [m²/s]

% Insulation material (glass foam granulate)
SW(3,1) = 0.08;                                 % Thermal conductivity λ of insulation [W/mK]
SW(3,2) = 170;                                  % Density ρ of insulation [kg/m³]
SW(3,3) = 850;                                  % Specific heat capacity cp of insulation [J/kgK]
SW(3,4) = round(sqrt(SW(3,1)*SW(3,2)*SW(3,3)),10); % Thermal diffusivity coefficient √(λρcp) [J/m²K·s^0.5]
SW(3,5) = SW(3,1)/(SW(3,2)*SW(3,3));            % Thermal diffusivity α = λ/(ρcp) [m²/s]


%%%------------------------------------------%%%
% 06 Discretization
%%%------------------------------------------%%%

dz = 0.005;                              % Vertical grid spacing [m]

d(1) = d_ST;                             % Cylinder diameter [m]
d(2) = round(0.005,3);                   % Radial grid spacing [m]
d(3) = d_ps_1D;                          % 1D equivalent penstock diameter [m]
d(4) = round(0.005,3);                   % Vertical grid spacing for 1D representation [m]

H(3)  = h_pist;                          % Piston height [m]
H(11) = h_lift;                          % Maximum lifting height [m]
H(7)  = dz * round(h_1D_tot/dz);         % Upper dead-zone incl. dome volume [m]

A(1) = d(1)^2*pi/4;                      % Cylinder cross-section area [m²]
A(2) = r_pist^2*pi;                      % Piston area [m²]
A(3) = A(1) - A(2);                      % Annular gap (ring gap) area [m²]

%--------------------------------------------------------------
% 6.1. Horizontal Geometry
%--------------------------------------------------------------

SoC_pot(1) = 1;                           % State of charge during initialization (time step 0)

H(1) = 2*H(7) + H(3) + H(11);             % Total cylinder height: 2 dead zones + piston + lifting height [m]
Nz(1) = round(1 + H(1)/dz);               % Number of grid points in full vertical domain
z_SSys = (0:dz:H(1));                     % Vertical coordinate vector (full system)

H(2) = H(7) + SoC_pot(1)*H(11);           % Height of lower pressure zone [m]
Nz(2) = round(1 + H(2)/dz);               % Grid points in lower pressure zone
z_UD = (0:dz:H(2));                       % Vertical vector lower pressure zone

Nz(3)  = round(1 + H(3)/dz);              % Grid points within piston
Nz(11) = 1 + ((Nz(3)-1)/2);               % Grid points in half-piston
z_Pist = (H(2):dz:H(2)+H(3));             % Vertical vector piston

H(4) = H(1) - H(2) - H(3);                % Upper pressure zone height [m]
Nz(4) = Nz(1) - Nz(2) - Nz(3) + 2;        % Grid points upper pressure zone
z_OD = (H(2)+H(3):dz:H(1));               % Vertical vector upper pressure zone

H(5) = round((A(3)*H(3))/A(1),2);         % Equivalent height of suction/pressure pipe (replacement volume) [m]
Nz(5) = round(1 + H(5)/dz);               % Grid points in replacement volume

H(6) = H(4) + H(5) + H(2);                % Total water inventory height [m]
Nz(6) = round(1 + H(6)/dz);               % Grid points water domain
z_W = (0:dz:H(6));                        % Vertical vector water inventory

H(8) = h_soil;                            % Soil height [m]
Nz(8) = round(1 + H(8)/dz);               % Grid points soil domain

H(9) = h_insu;                            % Insulation layer height [m]
Nz(9) = round(1 + H(9)/dz);               % Grid points insulation domain

H(10) = H(3) + H(8) + H(6) + H(9);        % Total height piston + free volume + soil + water + insulation [m]
Nz(10) = round(2 + H(10)/dz);             % Grid points for full vertical system
z_Sys = 0:dz:H(10)+dz;                    % Vertical coordinate vector full system

H(13) = H(2) + H(3) + H(4) + H(9);        % Height soil + piston + pressure zones + insulation [m]
Nz(13) = round(1 + H(13)/dz);             % Grid points vertical soil zone
z_ZE = 0:dz:H(13);                        % Vertical coordinate soil-only domain


%--------------------------------------------------------------
% 6.2. Radial Geometry
%--------------------------------------------------------------

Nz(12) = 13;                              % Number of radial nodes
z_RE = zeros(5, Nz(12));                  % Matrix for radial geometry/storage

z_RE(1,1) = d(1)/2;                       % Inner radius (outer cylinder radius)

% Hard-coded radial node distribution for soil domain [m]
z_RE(1,:) = [d(1)/2 , ...
             d(1)/2+0.005 , d(1)/2+0.015 , d(1)/2+0.035 , d(1)/2+0.075 , ...
             d(1)/2+0.15 , d(1)/2+0.3 , d(1)/2+0.6 , d(1)/2+1.2 , ...
             d(1)/2+2.4 , d(1)/2+4.6 , d(1)/2+6.8 , d(1)/2+9];


%--------------------------------------------------------------
% 6.3. Area of Single Circular Ring Segments 
%--------------------------------------------------------------

z_RE(2,1) = round(pi*(z_RE(1,2)-z_RE(1,1)) * ((z_RE(1,2)-z_RE(1,1))/4 + z_RE(1,1)), 10); 
                                             % Area of first circular ring segment [m²]

z_RE(2,end) = round(-pi*(z_RE(1,end-1)-z_RE(1,end)) * ...
                         ((z_RE(1,end-1)-z_RE(1,end))/4 + z_RE(1,end)), 10);
                                             % Area of last circular ring segment [m²]

for i = 2:Nz(12)-1
    z_RE(2,i) = round( ...
        pi*(z_RE(1,i+1)-z_RE(1,i-1)) * ...
        (z_RE(1,i-1) + 0.5*((z_RE(1,i+1)-z_RE(1,i-1))/2 + z_RE(1,i) - z_RE(1,i-1))), ...
    10);                                      % Area of intermediate ring segments [m²]
end

for i = 2:length(z_RE)-1
    z_RE(3,i) = round(z_RE(1,i+1)-z_RE(1,i-1), 10); % Spacing to next/previous node (symmetric) [m]
    z_RE(4,i) = round(z_RE(1,i+1)-z_RE(1,i),   10); % Forward node spacing [m]
    z_RE(5,i) = round(z_RE(1,i)-z_RE(1,i-1),   10); % Backward node spacing [m]
end    

dr = d(2);                                % Radial grid spacing [m]


%--------------------------------------------------------------
% 6.4. Temporal Discreization
%--------------------------------------------------------------

dt  = 900;                                 % Time step size [s]
Nt2 = 2;                                    % Number of sub-steps per procedure
t2  = (0:dt/Nt2:dt);                        % Time discretization for a single procedure



%%%------------------------------------------%%%
% 07 Radial heat-loss factors (new analytical approximation)
%%%------------------------------------------%%%

RE   = d(1)/2;      % Outer cylinder radius [m]
R_PS = d(3)/2;      % Penstock (1D equivalent) radius [m]


SW(4,1) = (2*SW(3,1)) / (RE^2 * SW(3,2) * SW(3,3) * log((RE+dr)/RE));
                                             % Radial loss factor for insulation [1/s]

SW(4,2) = (2*SW(2,1)) / ((RE^2 - R_PS^2) * SW(2,2) * SW(2,3) * log((RE+dr)/RE));
                                             % Radial loss factor for soil (referenced to piston radius) [1/s]

SW(4,3) = (2*SW(2,1)) / (RE^2 * SW(1,2) * SW(1,3) * log((RE+dr)/RE));
                                             % Radial loss factor for soil (referenced to water domain) [1/s]


%%%------------------------------------------%%%
% 08 Initial Conditions and Time Marching
%%%------------------------------------------%%%

% Initial temperatures in the system [°C]
T0init(1) = 45;                      % Initial bulk storage temperature
T0init(2) = 80;                      % Supply (feed) temperature
T0init(3) = 45;                      % Return temperature
T0init(4) = 20;                      % Ambient air temperature
T0init(5) = 11;                      % Undisturbed ground temperature
T0init(6) = T(3,1) - 273.15;         % Initial boundary temperature for ground (from external data)
T0init(7) = 11;                      % Initial insulation temperature


%--------------------------------------------------------------
% Initialization of 1D temperature field
%--------------------------------------------------------------

IC_Sys = T0init(1) * ones(1, Nz(10) + Nz(12)*Nz(13));      % Initial temperature for full 1D system

IC_Sys(1:Nz(3)) = T0init(5);                               % Initial piston temperature
IC_Sys(Nz(3)+1 : Nz(3)+Nz(8)) = T0init(5);                 % Initial soil temperature below the storage

IC_Sys(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4) = T0init(4); 
                                                           % Initial air temperature (single air node)

IC_Sys(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2 : ...
       Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-5) = T0init(7); % Initial insulation temperature

IC_Sys(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-3 : end) = T0init(5);
                                                           % Initial temperature in radial soil domain


% Helper matrices for post-processing
T_V   = zeros(1, Nz(10));                                   % Vertical system temperature profile
T_W   = zeros(1, Nz(6));                                    % Water volume temperature profile
T_REf = T0init(5) * ones(Nz(13), Nz(12));                   % Radial soil temperature field (vertical x radial nodes)


%--------------------------------------------------------------
% Time loop: conductive pre-initialization
%--------------------------------------------------------------

filenameSIM = ['Init_d' num2str(d_ST) '_h' num2str(h_pist) '_g' num2str(r_gap) '.mat'];
o = 1;

for i = 1:length(time_long)
    tic

    % Example for varying boundary conditions (kept as reference):
    % T0init(4) = ZR_900(k,6);   % Ambient temperature
    % Nz(2)     = ZR_900(k,11);  % Grid points lower pressure zone
    % Nz(4)     = ZR_900(k,12);  % Grid points upper pressure zone
    % flow      = h_1D/(120*3600); % Flow velocity in storage [m/s] for charging/discharging

    T0init(6) = T(3,i) - 273.15;                        % Updated ground boundary temperature [°C]

    % Compute transient heat conduction and heat transfer
    [T_Sys, T_REf] = HeattransferInit( ...
        t2, IC_Sys, Nz, dz, flow, T0init, SW, A, z_RE, T_REf, Nt2);

    % Update initial condition for next time step (use last time level)
    IC_Sys = T_Sys(end,:);

    % Extract vertical system temperature profile
    T_V(1,:) = T_Sys(end,1:Nz(10));

    % Extract water inventory temperature profile
    T_W(1,:) = T_Sys(end, Nz(3)+Nz(8) : Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3);

    toc


    %--------------------------------------------------------------
    % Periodic saving of intermediate simulation results
    %--------------------------------------------------------------
    if i == 28800 || i == 41800 || i == 61400 || i == 77500 || ...
       i == 95800 || i == 112800 || i == 130800 || i == 147057

        filenameSIM = ['Init_d' num2str(d_ST) ...
                       '_h' num2str(h_pist) ...
                       '_time' num2str(o) '.mat'];

        % Save current simulation state
        save(filenameSIM, "H", "z_OD", "z_UD", "z_W", "z_Pist", ...
                         "z_SSys", "z_Sys", "IC_Sys", "T_REf", ...
                         "d", "SW", "Nz", "T", "z_RE", "T_V")

        o = o + 1;
    end
end
