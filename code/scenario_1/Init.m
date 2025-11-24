clc 
clear


%%%------------------------------------------%%%
% Setup script for initialization of the thermal storage model
%%%------------------------------------------%%%
% TODISCUSS - Timestep length
t0 = 0;
t1 = 2890*3600; % 120 days in s
t_end = 35064*3600+3600*1700; % 1532 days in s or 4.2 years

time2 = t0:900:t1;
time = t0:900:t_end;

T = zeros(3,length(time));

for k = 1:length(time2)
T(3,k) = 11+273.15 + time(k)*(286.25-(11+273.15))/(t1);
end
% TODISCUSS - Welche 3 Starttemperaturen - Randtemperaturen für alle Zeitschritte?
% Initial temperature + 
for k = length(time2):length(time)
T(1,k) = 11+79*(1-exp(-time(k)/(35064*3600/3))); % an 1
T(2,k) = 11+44*(1-exp(-time(k)/(35064*3600/3))); % an 1
T(3,k) = 264.75+((T(1,k)-T(2,k))/2)*sin(time(k)*(2*pi)/(3600*8766)+3600*8764/4)+((T(1,k)-T(2,k))/2)+T(2,k);
end

%%%------------------------------------------%%%
% Speicher Basis geometrie
%%%------------------------------------------%%%
d_ST = 18; %Durchmesser in m
r_ST = d_ST/2; %Radius in m
h_tot = 1; % Totzone in m
h_lift = 16; % Hubhöhe in m
h_Water = h_lift + 2*h_tot; % zylindrische Wassersauele in m
V_WZ = h_Water*pi/4*d_ST^2; % Volumen zylindrische Wassersaeule in m³

%Prallplatte Einlass/Auslass
d_in = 3.6; % Durchmesser Prallplatte in m
r_in = d_in/2; % Radius Prallplatte in m
h_in = 0.05; % Höhe Prallplatte in Meter
V_PP = h_in*d_in^2/4*pi; % Volumen Prallplatte in m³

% Kuppelgeometrie
h_kupp = 1; % Höhe Kuppel in m
V_Kupp = (2/3)*r_ST^2*h_kupp*pi; % Volumen einer halben Kuppel in m³

% Kolbengeometrie
h_pist = 18; % Kolbenhöhe in m
r_gap = 0.5; % Ringspalt in m
r_pist = r_ST-r_gap; % Kolbenradius in m
V_pist = h_pist*pi*r_pist^2; % Volumen Kolben in m³
V_gap = h_pist*pi*(r_ST^2-r_pist^2); % Volumen Ringspalt in m³;

% Umgebung
h_insu = 2; % Isolationhöhe in m
r_soil = 9; % Radiales erdreich in m 
h_soil = 1+r_soil; % Höhe vertikales Erdreich

% Membran
h_mem = 0.05; % Höhe der Membran in m
V_mem = h_mem*pi*(r_ST^2-r_pist^2); % Volumen der Membran in m

% Druckstollen
r_out = 0.05; % Radius Kreisring in m
r_out_gap = 0.01; % Radius Auslass negativ in m;
h_out = 0.045; % Höhe Auslass
V_To = pi^2*r_ST*r_out^2; % Volumen Kreisring in m³ (Approx)
V_Out = h_out*pi*((r_ST+r_out_gap)^2-r_ST^2); % Volumen Auslass negativ in m³

% Volumen Wasser
V_W_gap = 2*(V_To-V_Out)-V_mem+V_gap; % Wasservolumen im Ringspalt in m³
V_W_ST = 2*(V_Kupp-V_PP)+V_WZ; % Wasservolumen im Zylinder und Kuppel in m³
V_W = V_W_gap+V_W_ST; % Gesamtes Wasservolumen in m³

E_pot = V_pist*(2400-1000)*9.81*h_lift/3600000; % m*g*h [J]
E_HWS = V_W*1000*4200*(80-11)/3600000; %V_w*rho*cp*Dt [J]

% Äquvivalente
h_1D_zyl = V_W_ST*4/(d_ST^2*pi); % Höhen Äquvivalent des zylindrischen Teils in m
h_1D_ps = V_W_gap*4/(d_ST^2*pi); % Höhen Äquvivalent des Druckstollens Teils in m
h_1D = h_1D_zyl+h_1D_ps; % Höhen Äquvivalent des gesamten Systems bei gleichen Querschnitt in m
d_ps_1D = (V_W_gap*4/(pi*h_pist))^0.5; % Durchmesser Äquvivalent für Druckstollen in m
h_1D_tot = h_tot+(V_Kupp-V_PP)/(d_ST^2*pi/4);

%%%------------------------------------------%%%
% Stoffwerte
%%%------------------------------------------%%%

% Wasser
SW(1,1) = 0.6; % Lambda Wasser in 'W/mK'
SW(1,2) = 1000; % Rho Wasser in 'kg/m³'
SW(1,3) = 4200; % Spez. Wärmekapazität Wasser in 'J/KgK'
SW(1,4) = round(sqrt(SW(1,1)*SW(1,2)*SW(1,3)),10); % Eindringkoeffizient Wasser 'J/m²Ks^0,5'
SW(1,5) = SW(1,1)/(SW(1,2)*SW(1,3)); % Lambda/(CP*Rho) Wasser in 'm²/s' 
% Erdmasse (an Beton angelehnt)
SW(2,1) = 2.1; % Lambda Erdmasse in 'W/mK'
SW(2,2) = 2400; % Rho Erdmasse in 'kg/m³'
SW(2,3) = 1000; % Spez. Wärmekapazität Erdwasser in 'J/KgK'
SW(2,4) = round(sqrt(SW(2,1)*SW(2,2)*SW(2,3)),10); % Eindringkoeffizient Erdmasse 'J/m²Ks^0,5'
SW(2,5) = SW(2,1)/(SW(2,2)*SW(2,3)); % Lambda/(CP*Rho) Erdmasse in 'm²/s'
% Dämmmaterial (an Glasschaum Granulat angelehnt)
SW(3,1) = 0.08; % Lambda Dämmmaterial in 'W/mK'
SW(3,2) = 170; % Rho Dämmmaterial in 'kg/m³'
SW(3,3) = 850; % Spez. Wärmekapazität Dämmmaterial in 'J/KgK'
SW(3,4) = round(sqrt(SW(3,1)*SW(3,2)*SW(3,3)),10); % Eindringkoeffizient Dämmmaterial 'J/m²Ks^0,5'
SW(3,5) = SW(3,1)/(SW(3,2)*SW(3,3)); % Lambda/(CP*Rho) Dämmmaterial in 'm²/s'


%%%------------------------------------------%%%
% Discretization
%%%------------------------------------------%%%

dz = 0.005; % verical grid spacing
d(1) = d_ST; % cylinder diameter
d(2) = round(0.005,3); % Ortsdiskretisierung in radialer Richtung
d(3) = d_ps_1D; % Durchmesser Penstock
d(4) = round(0.005,3); % Ortsdiskretisierung in veritkaler Richtung
H(3) = h_pist; % piston height
H(11) = h_lift; % maximal lifting height
H(7) = dz*round(h_1D_tot/dz); % upper dead-zone incl. dome

A(1) = d(1)^2*pi/4;% cylinder surface
A(2) = r_pist^2*pi; % surface piston
A(3) = A(1)-A(2); % surface ringspalt

% Horizontale Geometrie
                                                                            % Ortsdiskretisierung in 'm'
SoC_pot(1) = 1;                                                             % Speicherladezustand in der Initalisierung (Nullter Zeitschritt)

H(1) = 2*H(7)+H(3)+H(11);                                                   % Zylinderhöe (zwei mal Totzone+Kolben+Hub) in 'm'
Nz(1) = round(1+H(1)/dz);                                                   % Anzahl der Ortsschirtte im Gesamten System
z_SSys = (0:dz:H(1));                                                       % Ortsvektor Speichersystem

H(2) = H(7)+SoC_pot(1)*H(11);                                               % Höhe Untere Druckzone relativ zur Hubhöhe
Nz(2) = round(1+H(2)/dz);                                                   % Anzahl der Ortsschirtte in der untere Druckzone
z_UD = (0:dz:H(2));                                                         % Ortsvektor untere Druckzone

Nz(3) = round(1+H(3)/dz);                                                   % Anzahl der Ortsschritte des Kolbens
Nz(11) = 1+((Nz(3)-1)/2);                                                   % Anzahl der Ortsschritte des halben Kolbens
z_Pist = (H(2):dz:H(2)+H(3));                                               % Ortsvektor Kolben

H(4) = H(1)-H(2)-H(3);                                                      % Höhe Obere Druckzone relativ zur Hubhöhe/2
Nz(4) = Nz(1)-Nz(2)-Nz(3)+2;                                                % Anzahl der Ortsschirtte in der oberen Druckzone
z_OD = (H(2)+H(3):dz:H(1));                                                 % Ortsvektor obere Druckzone

H(5) = round((A(3)*H(3))/A(1),2);                                           % Höhe Druck- bzw. Saugrohr (Ersatzvolumen)
Nz(5) = round(1+H(5)/dz);                                                   % Anzahl der Ortsschirtte im Druck- bzw. Saugrohr (Ersatzvolumen)

H(6) = H(4)+H(5)+H(2);                                                      % Höhe des Wasservolumen untere Druckzone, Untere Druckzone, Ersatzvolumen  
Nz(6) = round(1+H(6)/dz);                                                   % Anzahl der Ortsschritte Wasservolumen (Obere Druckzone, Untere Druckzone, Ersatzvolumen)
z_W = (0:dz:H(6));                                                          % Ortsvektor Wasserinventar
                                                           
H(8) = h_soil;                                                              % Höhe des Erdreichs unterhalb des Speichersystems
Nz(8) = round(1+H(8)/dz);                                                   % Ortschritte des Erdreichs

H(9) = h_insu;                                                              % Höhe Dämmung in 'm'
Nz(9) = round(1+H(9)/dz);                                                   % Ortsschritte der Dämmung

H(10) = H(3)+H(8)+H(6)+H(9);                                                % Höhe des vertikalen Systems: Kolben, 1 Freivolumen ,Erdreich, Wasservolumen (untere Druckzone, Ersatzvolumen, obere Druckzone), Dämmung
Nz(10) = round(2+H(10)/dz);                                                 % Anzahl der Orsschritte vertikalen Systems Kolben, 1 Freivolumen ,Erdreich, Wasservolumen (untere Druckzone, Ersatzvolumen, obere Druckzone), Dämmung
z_Sys = 0:dz:H(10)+dz;                                                      % Ortsvektor des Systems

H(13) = H(2)+H(3)+H(4)+H(9);                                                % Höhe System Erdreich in 'm'
Nz(13) = round(1+H(13)/dz);                                                 % Anzahl der vertikalen Ortsschritte Erdreich
z_ZE = 0:dz:H(13);        


% Diskretisierung der radialen Geometrie
Nz(12) = 13;                                                                % Anzahl Knoten in raidaler Richtung
z_RE = zeros(5,Nz(12));                                                     % Radialer Vektor für das Erdreich
z_RE(1,1) = d(1)/2;                                                         % Radius für Außenzylinder
%dr = d(2);                                                                  % Ortsdiskretisierung in radialer Richtung
% for i = 0:Nz(12)-2
%     z_RE(1,i+2) = z_RE(1,i+1)+dr*2^i;                                       % radialer Vektor des Erdreichs in 'm'
% end
z_RE(1,:)= [d(1)/2 , d(1)/2+0.005 , d(1)/2+0.015 , d(1)/2+0.035 , d(1)/2+0.075 , d(1)/2+0.15 , d(1)/2+0.3 , d(1)/2+0.6 , d(1)/2+1.2 , d(1)/2+2.4 , d(1)/2+4.6 , d(1)/2+6.8 , d(1)/2+9];

% Flächeninhalte der einzelnen Kreisringsegmente
z_RE(2,1)=round(pi*(z_RE(1,2)-z_RE(1,1))*(((z_RE(1,2)-z_RE(1,1))/4+z_RE(1,1))),10); % Flächeninhalt des ersten Kreisringsegmentes in 'm²'
z_RE(2,end) = round(-pi*(z_RE(1,end-1)-z_RE(1,end))*((z_RE(1,end-1)-z_RE(1,end))/4+z_RE(1,end)),10); % Flächeninhalt des letzten Kreisringsegmentes in 'm²'
for i = 2:Nz(12)-1
    z_RE(2,i) = round(pi*(z_RE(1,i+1)-z_RE(1,i-1))*(z_RE(1,i-1)+0.5*((z_RE(1,i+1)-z_RE(1,i-1))/2+z_RE(1,i)-z_RE(1,i-1))),10); % Flächeninhalte der Kreisringsegmente in 'm²'
end
for i = 2:length(z_RE)-1
    z_RE(3,i)= round(z_RE(1,i+1)-z_RE(1,i-1),10); % Abstand zwischen den vorherigen und den nachfolgenden Knotenpunkt 'm'
    z_RE(4,i)= round((z_RE(1,i+1)-z_RE(1,i)),10); % Abstand zwischen den vorherigen Knotenpunkt in 'm'
    z_RE(5,i)= round((z_RE(1,i)-z_RE(1,i-1)),10); % Abstand zwischen den nachfolgenden Knotenpunkt in 'm'
end    

dr = d(2); % Ortsdiskretisierung in radialer Richtung

% Zeitdiskretisierung für Berechnungsschritt
dt = 900;                                                                   % Zeitschritt relativ zur vollen Stunde, in 's'
Nt2 = 2;                                                                    % Anzahl der Zeitschritte einer Prozedur
t2 = (0:dt/Nt2:dt);                                                        % Zeitdiskretisierung einer Prozedur

RE = d(1)/2;
R_PS = d(3)/2;

SW(4,1) = (2*SW(3,1))/(RE^2*SW(3,2)*SW(3,3)*log((RE+dr)/RE)); % ganz neuer Ansatz Verlustfaktor radiale Dämmung '1/s'
SW(4,2) = (2*SW(2,1))/((RE^2-R_PS^2)*SW(2,2)*SW(2,3)*log((RE+dr)/RE)); % ganz neuer Ansatz Verlustfaktor radiales Erdreich im verhältnis zu kolben '1/s'
SW(4,3) = (2*SW(2,1))/(RE^2*SW(1,2)*SW(1,3)*log((RE+dr)/RE)); % ganz neuer Ansatz Verlustfaktor radiales Erdreich '1/s' im verhältnis zu wasser


% % Anfangsbedingungen im Systems, in '°C'
% Initialtemperaturen
T0init(1) = 45; % Starttemperatur                                                            
T0init(2) = 80; % Vorlauftemperatur
T0init(3) = 45; % Ruecklauftemperatur
T0init(4) = 20; % Lufttemperatur
T0init(5) = 11; % Erdtemperatur
T0init(6) = T(3,1)-273.15; %Randtemperatur für Initalbetrachtung
T0init(7) = 11; % Initialtemperatur Dämmung

% INITIALISIERUNG DER ANFANGSWERTE
IC_Sys = T0init(1)*ones(1,Nz(10)+Nz(12)*Nz(13));                            % Starttemperatur für Gesamte 1D Simulation
IC_Sys(1:Nz(3)) = T0init(5);                                                % Starttemperatur des Kolbens der Gesamten 1D Simulation
IC_Sys(Nz(3)+1:Nz(3)+Nz(8)) = T0init(5);                                    % Starttemperatur im Boden (Unterhalb vom Speicher) der Gesamten 1D Simulation
IC_Sys(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4) = T0init(4);                  % Starttemperatur der Luft der Gesamten 1D Simulation
IC_Sys(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2:Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-5) = T0init(7); % Starttemperatur Dämmung
IC_Sys(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-3:end) = T0init(5);              % Starttemperatur im radialen Erdreich

T_V = zeros(1,Nz(10));                                                      % Berechnungsmatrix für vertikales System
T_W = zeros(1,Nz(6));                                                       % Berechnungsmatrix für Wasservolumen
T_REf = T0init(5)*ones(Nz(13),Nz(12));                                      % Matrix für Tensor raidales Erdreich radiale Gitterpunkte x vertikale Gitterpunkte.

filenameSIM = ['Init_d' num2str(d_ST) '_h' num2str(h_pist) '_g' num2str(r_gap) '.mat'];
o=1;
for i = 1:length(time)
tic
    %T0init(4) = ZR_900(k,6); % Umgebungstemperatur festlegen
    %Nz(2) = ZR_900(k,11); % Anzahl der Ortsschritte untere Druckzone
    %Nz(4) = ZR_900(k,12); % Anzahl der Ortsschritte obere Druckzone
    %flow = h_1D/(120*3600); % Strömungsgeschwindigkeit im Speichersytem, in 'm/s', für thermisches Be- und Entladen berechnen
    T0init(6) = T(3,i)-273.15; %Randtemperatur für Initalbetrachtung 

    % Temperaturübertragung und Stofftransport berechnen
    [T_Sys,T_REf] = HeattransferInit(t2,IC_Sys,Nz,dz,flow,T0init,SW,A,z_RE,T_REf,Nt2);

    % Neusetzen der Anfangsbedingung für nächsten Prozedurschritt
    IC_Sys = T_Sys(end,:);            
    % Setzen der Temperaturen im vertikalen System
    T_V(1,:) = T_Sys(end,1:Nz(10));
    % Berechnungsmatrix für Wasservolumen
    T_W(1,:) = T_Sys(end,Nz(3)+Nz(8):Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3);

toc

if i == 28800 || i == 41800 || i == 61400 || i == 77500 || i == 95800 || i == 112800 || i== 130800 || i == 147057
filenameSIM = ['Init_d' num2str(d_ST) '_h' num2str(h_pist) '_time' num2str(o) '.mat'];    
% Speichern des Simulationsergebnisses
save(filenameSIM,"H","z_OD","z_UD","z_W","z_Pist","z_SSys","z_Sys","IC_Sys","T_REf","d","SW","Nz","T","z_RE","T_V")
o = o+1;
end

end

