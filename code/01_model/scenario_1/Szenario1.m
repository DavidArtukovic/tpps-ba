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
%   - For each timestep, calls the coupled heat transfer
%     model (HeattransferSzen) and updates the system state.
%   - Stores vertical temperature profiles, water temperatures, and
%     energy/exergy-related result vectors for later evaluation and plotting.
%
% ---------------------------------------------------------------

clc 
clear
%%%------------------------------------------%%%
% 01. Load Scenario and Initilization Simulation
%%%------------------------------------------%%%

% Load local paths (per-machine config)
run(fullfile('..','..', 'configs', 'paths_local.m'));

% Now build data subfolder for this configuration:
DATA_SCEN1 = fullfile(DATA_BASE, 'Matlab_d18_h18_Szenario1');
%%
% Load init and scenario files
load(fullfile(DATA_SCEN1, 'Init_d18_h18_time8.mat')); % Load geometry, material values and initial values from initial simulation
load(fullfile(DATA_SCEN1, 'SzenarioComsol.mat')); % Load the scenarion

%%%------------------------------------------%%%
% 02. Initialize Arrays
%%%------------------------------------------%%%
%%
% Result arrays
Res_900 = zeros(10,length(t_900)+1); % Values for quarter hours
Res_hour = zeros(10,length(t_hour)+1); % Values for full hours

% 1D Temperature vector for whole system (soil zone)
T_SSys = zeros(Nz(13),1);

T_V_900 = zeros (Nz(13),length(t_900)); % Temperature array, consisting of grid(system) x time
T_W_900 = zeros (Nz(6),length(t_900)); % Temperature array, consisting of grid(water) x time

% Step width spacial discretization
dz = d(4); % Spacial discretization in vertical direction
dr = d(2); % Spacial discretization in radial direction

% Time discretization
dt = 900; % Time in seconds for every procedure step
Nt2 = 2; % Number of timesteps for every procedure step
t2 = (0:dt/Nt2:dt);
K = 273.15; % reference temperature for exergy calculation
Spz = round(SW(1,2)*SW(1,3)*d(1)^2*pi/4,10); 

d_ST = d(1); % diameter in m
r_ST = d_ST/2; % radius in m
r_gap = 0.5; % Ringspalt in m
r_pist = r_ST-r_gap; % piston radius in m
A(1) = d(1)^2*pi/4;% area cylinder
A(2) = r_pist^2*pi; % area piston
A(3) = A(1)-A(2); % area ringspalt

% Charge and Discharge counters
LC = 0;
DC = 0;

for k = 1:length(t_900)
    Res_900(7,k+1) = round(1+(H(7)+t_900(2,k)*H(11))/d(4)); % position of the piston end in 
    Res_900(8,k+1) = Nz(1)-Res_900(7,k+1)-Nz(3)+2; % Koordinate der Oberen Druckzone
    if t_900(1,k) > 0
        LC = LC+1;
    elseif t_900(1,k) < 0
        DC = DC+1; 
    end
end
LC = LC/4; % Volle Ladestunden
DC = DC/4; % Volle Entladestunden

% Average flow for a full charge or discharge cycle
Lflow = H(6)/(LC*3600); % average charge flow in m/s
Dflow = H(6)/(DC*3600); % average discharge flow in m/s

% Average charge and discharge flow volume in meter
for k = 1:length(t_900)
    if t_900(1,k) > 0
        Res_900(9,k+1) = -Lflow*t_900(1,k); % Ladegeschwindigkeit in m/s per Zeitschritt
    elseif t_900(1,k) < 0
        Res_900(10,k+1) = -Dflow*t_900(1,k); % Entladegeschindigkeit in m/s per Zeitschritt
    end 
end

% Set tempererature vector to initial state
T_Sys = IC_Sys;

DTRE = zeros(1,Nz(12)); % vector of radial temperature differences in radial soil
% temperature in whole vertical system
T_V(1,:) = T_Sys(1,1:Nz(10));
% temperature in water system
T_W(1,:) = T_Sys(1,Nz(3)+Nz(8):Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3);

% Exergy density of the water
wEX = zeros(1,length(T_W));

% Thermal energy in insulation (vertical) [J] relative to initial time stamp.
Heat_insu = SW(3,2)*SW(3,3)*dz*A(1)...
*(sum(T_V(1,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2:end-1))...
+0.5*(T_V(1,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3)+T_V(1,end))); 

% Thermal energy in the water column [J]
Heat_Wasser = SW(1,2)*SW(1,3)*dz*A(1)...
*(sum(T_W(1,2:end-1))+0.5*(T_W(1,1)+T_W(1,end)));

% Thermal energy in the piston [J]
Heat_piston = SW(2,2)*SW(2,3)*dz*A(2)*(sum(T_Sys(end,2:Nz(3)-1))...
+0.5*(T_Sys(end,1)+T_Sys(end,Nz(3))));

% Thermal energy in vertical soil below the system [J]
Heat_Vsoil = SW(2,2)*SW(2,3)*dz*A(1)...
*(sum(T_V(1,Nz(3)+2:Nz(3)+Nz(8)-1))...
+0.5*(T_V(1,Nz(3)+1)+T_V(1,Nz(3)+Nz(8))));


% Thermal energy in radial soil (before insulation layer)
% T_REf is the radial x vertical mesh in the soil initialized with 11°C
% is counts from bottom to top, from in to out.
for m = 1 : Nz(12)
    DTRE(m) = dz*z_RE(2,m)*(sum(T_REf(2:end-Nz(9),m))+0.5*(T_REf(1,m)+T_REf(end-Nz(9)+1,m))); %'K*m³
end

% Thermal Energy in radial soil [J] relative to initial temperature.
Heat_Rsoil = SW(2,2)*SW(2,3)*sum(DTRE); 

% Thermal energy in radial insulation [J] relative to initial temperatuer.
for m = 1 : Nz(12)
    DTRE(m) = dz*z_RE(2,m)*(sum(T_REf(end-Nz(9)+2:end-1,m))+0.5*(T_REf(end-Nz(9)+1,m)+T_REf(end,m))); %'K*m³
end

Heat_Rinsu = SW(3,2)*SW(3,3)*sum(DTRE); 

% With z_W as grid of the water volume.
for p = 1:length(z_W) % Voranalyse des Wasservolumens
    wEX(1,p) = T_W(p)-K*log((K+T_W(p))/(K)); %  Exergy density in total water volume referenced to temperature K (273,15 K) 
end 
WEX = Spz*dz*(sum(wEX(1,2:end-1))+0.5*(wEX(1,1)+wEX(1,end))); % Totla exergy in total water volume referenced to temperature K (273,15 K) 

Res_hour(1,1) = Heat_insu;
Res_hour(2,1) = Heat_Wasser;
Res_hour(3,1) = Heat_piston;
Res_hour(4,1) = Heat_Vsoil;
Res_hour(5,1) = Heat_Rsoil;
Res_hour(6,1) = Heat_Rinsu;
Res_hour(11,1) = WEX;

Res_900(1,1) = Heat_insu;
Res_900(2,1) = Heat_Wasser;
Res_900(3,1) = Heat_piston;
Res_900(4,1) = Heat_Vsoil;
Res_900(5,1) = Heat_Rsoil;
Res_900(6,1) = Heat_Rinsu;
Res_900(11,1) = WEX;
                                                      
% % Anfangsbedingungen im Systems, in '°C'
% Initialtemperaturen
T0init(1) = 45; % Starttemperatur                                                            
T0init(2) = 80; % Vorlauftemperatur
T0init(3) = 45; % Ruecklauftemperatur
T0init(4) = 20; % Lufttemperatur
T0init(5) = 11; % Erdtemperatur
T0init(6) = T(3,1)-273.15; %Randtemperatur für Initalbetrachtung
T0init(7) = 11; % Initialtemperatur Dämmung

%%
% t_900/2 für 10 Tage.
for i = 1:length(t_900/2)
tic
    %T0init(4) = ZR_900(k,6); % Umgebungstemperatur festlegen
    Nz(2) = Res_900(7,i+1); % Anzahl der Ortsschritte untere Druckzone
    Nz(4) = Res_900(8,i+1); % Anzahl der Ortsschritte obere Druckzone
    flow = Res_900(9,i+1)+Res_900(10,i+1); % Strömungsgeschwindigkeit im Speichersytem, in 'm/s', für thermisches Be- und Entladen berechnen
    
    % Temperaturübertragung und Stofftransport berechnen
            [T_Sys,T_REf] = HeattransferSzen(t2,IC_Sys,Nz,dz,flow,T0init,SW,A,z_RE,T_REf,Nt2);

            % Neusetzen der Anfangsbedingung für nächsten Prozedurschritt
            IC_Sys = T_Sys(end,:);            
            % Setzen der Temperaturen im vertikalen System
            T_V(1,:) = T_Sys(end,1:Nz(10));
            % Berechnungsmatrix für Wasservolumen
            T_W(1,:) = T_Sys(end,Nz(3)+Nz(8):Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3);

            % Berechnen der Energieinhalte
            
    	    Heat_insu = SW(3,2)*SW(3,3)*dz*A(1)*(sum(T_V(1,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2:end-1))+0.5*(T_V(1,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3)+T_V(1,end))); %Thermische Energie im Wasservolumen in 'J' relativ zum Anfangszeitpunkt;
            Heat_Wasser = SW(1,2)*SW(1,3)*dz*A(1)*(sum(T_W(1,2:end-1))+0.5*(T_W(1,1)+T_W(1,end))); %Thermische Energie im Wasservolumen in 'J' relativ zum Anfangszeitpunkt;
            Heat_piston = SW(2,2)*SW(2,3)*dz*A(2)*(sum(T_Sys(1,2:Nz(3)-1))+0.5*(T_Sys(end,1)+T_Sys(end,Nz(3)))); %Thermische Energie im Kolben in 'J' relativ zum Anfangszeitpunkt;
            Heat_Vsoil = SW(2,2)*SW(2,3)*dz*A(1)*(sum(T_V(1,Nz(3)+2:Nz(3)+Nz(8)-1))+0.5*(T_V(1,Nz(3)+1)+T_V(1,Nz(3)+Nz(8)))); %Thermische Energie im Erdreich unterhalb das Systems in 'J' relativ zum Anfangszeitpunkt;
            % Thermische Energie im radialen Soil bis Dämmung
            for m = 1 : Nz(12)
                DTRE(m) = dz*z_RE(2,m)*(sum(T_REf(2:end-Nz(9),m))+0.5*(T_REf(1,m)+T_REf(end-Nz(9)+1,m))); %'K*m³
            end
            Heat_Rsoil = SW(2,2)*SW(2,3)*sum(DTRE); %Thermische Energie im radialen Erdreich in 'J' relativ zum Anfangszeitpunkt;
            % Thermische Energie in der radialen Dämmung
            for m = 1 : Nz(12)
                DTRE(m) = dz*z_RE(2,m)*(sum(T_REf(end-Nz(9)+2:end-1,m))+0.5*(T_REf(end-Nz(9)+1,m)+T_REf(end,m))); %'K*m³
            end
            Heat_Rinsu = SW(3,2)*SW(3,3)*sum(DTRE); %Thermische Energie im radialen Erdreich in 'J' relativ zum Anfangszeitpunkt;

            for p = 1:length(z_W) % Voranalyse des Wasservolumens
                wEX(1,p) = T_W(p)-K*log((K+T_W(p))/(K)); % Exergie im gesamten Wasservolumen auf RÃ¼cklauftemperatur bezogen ohne aktive Nullung
            end
            WEX = Spz*dz*(sum(wEX(1,2:end-1))+0.5*(wEX(1,1)+wEX(1,end)));

            

            % Berechnungswerte der Viertelstundenaufnahmen
            Res_900(1,1+i) = Heat_insu;
            Res_900(2,1+i) = Heat_Wasser;
            Res_900(3,1+i) = Heat_piston;
            Res_900(4,1+i) = Heat_Vsoil;
            Res_900(5,1+i) = Heat_Rsoil;
            Res_900(6,1+i) = Heat_Rinsu;
            Res_900(11,1+i) = WEX; %Exergie im Wasservolumen

            T_W_900(:,i) = flip(T_W);

            T_SSys(1:Nz(2)) = T_V(Nz(3)+Nz(8):Nz(3)+Nz(8)+Nz(2)-1); % Untere Druckzone besetzen

            T_SSys(Nz(2)+Nz(3)-1:Nz(2)+Nz(3)+Nz(4)-2) = T_V(Nz(3)+Nz(8)+Nz(2)+Nz(5)-2:Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3); % Obere Druckzone besetzen

            T_SSys(Nz(2):Nz(2)+Nz(3)-1) = T_V(1:Nz(3)); % Kolben besetzen
            T_SSys(Nz(2)+Nz(3)+Nz(4)-2:end) = T_V(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3:end); % Dämmung besetzen
            T_V_900(:,i) = flip(T_SSys);


            if mod((i),4) == 0 % Werte für die Stundenaufnahmen
                Res_hour(1,1+i/4) = Heat_insu; % Thermische Energie in der vertikalen Dämmung
                Res_hour(2,1+i/4) = Heat_Wasser; % Thermische Energie im Wasservolumen
                Res_hour(3,1+i/4) = Heat_piston; % Thermische Energie im Kolben
                Res_hour(4,1+i/4) = Heat_Vsoil; % Thermische Energie im vertiakalen Boden
                Res_hour(5,1+i/4) = Heat_Rsoil; % Thermische Energie im radialen Boden
                Res_hour(6,1+i/4) = Heat_Rinsu; % Thermische Energie in der radialen Dämmung
                Res_hour(7,1+i/4) = Res_900(7,1+i); % Anzahl der örtlichen Diskreten des unteren Wasservolumen
                Res_hour(8,1+i/4) = Res_900(8,1+i); % Anzahl der örtlichen Diskreten des oberen Wasservolumen
                Res_hour(9,1+i/4) = Res_900(9,1+i); % Ladegeschindigkeit in m/s
                Res_hour(10,1+i/4) = Res_900(10,1+i); % Entladegeschindigkeit in m/s
                Res_hour(11,1+i/4) = Res_900(11,1+i); % %Exergie im Wasservolumen

            end
    disp(i*15/60);
toc
end

Res_900_d18_18 = Res_900;
Res_hour_d18_18 = Res_hour;
Res_Wasser_d18_18 = T_W_900;
Res_System_d18_18 = T_V_900;


filenameSIM = ['d' num2str(d_ST) '_h' num2str(H(3)) '_Res_Matlab_d18_18.mat'];    
% Speichern des Simulationsergebnisses
% save(filenameSIM,"Res_900_d18_18","Res_hour_d18_18","Res_Wasser_d18_18","Res_System_d18_18")

