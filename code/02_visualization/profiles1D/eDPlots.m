
clc
clear
close all
%%%------------------------------------------%%%
% 01. Load Scenario and Initialization Simulation
%%%------------------------------------------%%%

% Load local paths (per-machine config)
run(fullfile('..','..','..', 'configs', 'paths_local.m'));

% Build data subfolder for this configuration
DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');
DATA_SCEN1_OWN = fullfile(DATA_BASE, 'scenario1_freeConv');

% Laden der Relevanten Daten
load(fullfile(DATA_SCEN1_BASE, 'd18_h18.mat'));
load(fullfile(DATA_SCEN1_BASE, 'd18_h18_Res_Matlab_d18_18.mat')); % Dominic Matlab Ergebnisse: "Temperatur über Systemhöhe alle 900 sek"
load(fullfile(DATA_SCEN1_BASE, '1D_05_TPPS_18_18.mat')); % Comsol Ergebnisse "Temperatur über Höhe an unterschiedlichen Punkten 
load(fullfile(DATA_SCEN1_OWN, '260108_d18_h18_Res_Matlab_FK_v2.mat'));  % Matlab Ergebnisse mit freier Konvektion v2
load(fullfile(DATA_SCEN1_OWN, '260105_d18_h18_Res_Matlab_FK_v1.mat')); % Matlab Ergebnisse mit freier Konvektion v1

%%
dz = 0.05; %Ortsdiskretisierung Comsol
o=0; %laufvariable für einzelne (Comsol) Simulation im Betriebsszenario (bis anzahl an Comsol ergebnissen)
k=0; %laufvariable für Viertelstunden -> Stunden (Matlab) 
l=0; %laufvariable für jeweilige Stunde im Betriebsszenario für ComsolSimulation


%%
% Schleife Matlab Ergebnisse auf Stundenschritte reduzieren
for i = 1:1920
    if mod(i,4) == 0
       k=k+1;
       T_Sys_M(:,k)=Res_System_d18_18(401:end,i); %Extrahieren der Temperaturwerte ohne Dämmung
    end
end

% Schleife Matlab Ergebnisse (DA free convection) auf Stundenschritte reduzieren
k_fk = 0;
for i = 1:1920
    if mod(i,4) == 0
       k_fk = k_fk + 1;
       T_Sys_M_fk(:,k_fk) = Res_System_d18_18_FK(401:end,i); % without insulation
    end
end

% Schleife Matlab Ergebnisse DA (no fk) auf Stundenschritte reduzieren
k_fk = 0;
for i = 1:1920
    if mod(i,4) == 0
       k_fk = k_fk + 1;
       T_Sys_M_nofk(:,k_fk) = Res_System_d18_18_noFK(401:end,i); % without insulation
    end
end

dz_M = 0.005; %Ortsdisketisierung Matlab
h_M = (length(T_Sys_M)-1)*dz_M; % Systemhöhe Matlab in m
h_Sys_C = (length(z)-1)*dz; % Systemhöhe Comsol in m
h_diff = abs(h_Sys_C-h_M); % Höhendifferenz der Modelle in m 
z_Mo = z(1)-h_diff/2; % oberer Punkt der z-Achse für Matlab
z_Mu = z(end)+h_diff/2; % untere Punkt der z-Achse für Matlab
z_M = flip(z_Mu:dz_M:z_Mo)'; % Z-Achse für Matlab Lösung

% Matlab auf Anzahl der 1D Profile Reduzieren
for i = 1:480
    if t_hour(3,i) >0
        o=o+1;
        l = l+t_hour(3,i);
        t_Sys_M(1,o) = l; %Array zeigt an, zu welche Stunde im Betriebsszenario das jeweilige Comsolergebnis passend ist
    end
end

% Matlabergebnisse auf Zeitpunkte reduzieren, an welchen Comsol ergebnisse vorhanden sint
T_M = zeros(length(T_Sys_M),length(t_Sys_M));
for i = 1:length(t_Sys_M)
    T_M(:,i) = T_Sys_M(:,t_Sys_M(1,i));
end

% Matlab (free convection) auf Zeitpunkte reduzieren, an welchen Comsol-Ergebnisse vorhanden sind
T_M_fk = zeros(length(T_Sys_M_fk),length(t_Sys_M));
for i = 1:length(t_Sys_M)
    T_M_fk(:,i) = T_Sys_M_fk(:,t_Sys_M(1,i));
end

% Matlab (free convection) auf Zeitpunkte reduzieren, an welchen Comsol-Ergebnisse vorhanden sind
T_M_nofk = zeros(length(T_Sys_M_fk),length(t_Sys_M));
for i = 1:length(t_Sys_M)
    T_M_nofk(:,i) = T_Sys_M_nofk(:,t_Sys_M(1,i));
end


% Kolbenposition zu den Zeitpunkten bestimmen
zp = zeros(1,length(t_Sys_M));
for i = 1:length(t_Sys_M)
    zp(1,i) = 1+(z(1)-z_p(1,i))/dz; %Kolbenposition im array
end

% Temperaturen aus Comsolsystem zusammensetzen
T_Sys_C = zeros(length(T1_W),length(t_Sys_M));
for i = 1:length(t_Sys_M)
    T_Sys_C(:,i) = T1_W(:,i); %Temperatur im Wasservolumen
    T_Sys_C(zp(1,i)+1:zp(1,i)+length(z_p)-2,i) = T1_P(2:end-1,i); %Temperatur im Kolben
end

% Array für Plots festlegen "Diese Ergebnisse werden geplottet"
% von 1 bis 21 Thermal Charge (Ungerade = Stilstand, Gerade = Thermisch Laden, Zeitabstände 12 stunden)
% von 22 bis bis 38  Potential Discharge (Zeitabstände eine Stunde zwischen Ergebnissen, außer bei 22 und 38)
% von 39 bis 55 Potential Charge (Zeitabstände eine Stunde zwischen Ergebnissen, außer bei 39 und 55)
% von 56 bis 72 Thermal Discharge (Zeitabstände 12 Stunden zwischen Ergebnissen, außer bei 56)

p = [4, 6, 8, 10, 12]; %"Plot mir die Ergebnisse von der Simulation ...."


%% Plot für 1D Temperaturverlauf über Arrayp
figure('Name','Temperaturverlauf','Color','w')
for m = 1:length(p)
    subplot(1,length(p),m)
    hold on
    %yline(z_Mo) %Begrenzung Matlab oben
    %yline(z_Mu) %Begrenzung Matlab unten
    %yline(0.35) %Inlet Comsol oben
    %yline(-36.35) %Inlet Comsol unten
    rectangle('Position',[40 z_p(end,p(m)) 40 18],'FaceColor',[1 1 1 0.1],'EdgeColor',[0 0 0]) % Kolben
    plot(T_Sys_C(:,p(m)),z,'-','Color',[0.3010 0.7450 0.9330],'LineWidth',1); %Comsol Ergebnisse
    plot(T_M(:,p(m)),z_M,'Color',[0.8500 0.3250 0.0980],'LineWidth',1); %Matlab Ergebnisse
    plot(T_M_fk(:,p(m)),z_M,'--','LineWidth',1); % Matlab Ergebnisse (free convection)
    plot(T_M_nofk(:,p(m)),z_M,'--','LineWidth',1); % Matlab Ergebnisse (free convection)
    grid on
    xticks([40, 50, 60, 70, 80])
    %Achse auf Matlab Koordinaten reduzieren (Comsol hat mehr Höhencoordinaten wegen Halbkugel oben und unten)
    axis([40 80 z_Mu z_Mo])
    % Beschriftung nur am ersten Plot
    if m == 1
        ylabel('System depth [m]','Interpreter','Latex','Fontsize',12)
        xlabel('Temperatur in TPPS [$^{\circ}$C] ','Interpreter','Latex','Fontsize',12)
    end
    if m == length(p)
        legend('$T_{Sys}$ Comsol','$T_{Sys}$ Matlab',...
            '$T_{Sys}$ Matlab (DA FK)$', '$T_{Sys}$ Matlab (DA no FK)$', ...
       'Interpreter','Latex','Fontsize',9)
    end
    sgtitle(['1D Temperature Profil in the 72 m diameter System at Times: ' num2str(t_Sys_M(p)) ] ,'Interpreter','Latex','Fontsize',12)
end




% %% Alternativer Plot mit Temperatur im Ringspalt (Comsol)
% figure('Name','Temperaturverlauf mit Ringspalt','Color','w')
% for m = 1:length(p)
%     subplot(1,length(p),m)
%     hold on
%     %yline(z_Mo) %Begrenzung Matlab oben
%     %yline(z_Mu) %Begrenzung Matlab unten
%     %yline(0.35) %Inlet Comsol oben
%     %yline(-36.35) %Inlet Comsol unten
%     plot(T1_W(:,p(m)),z,'-','Color',[0.3010 0.7450 0.9330],'LineWidth',1); %Comsol Ergebnisse
%     rectangle('Position',[40 z_p(end,p(m)) 40 18],'FaceColor',[1 1 1],'FaceAlpha',0.5,'EdgeColor',[0 0 0]) % Kolben
%     plot(T1_P(2:end,p(m)),z_p(2:end,p(m)),'--','Color',[0.3010 0.7450 0.9330],'LineWidth',1); %Comsol Ergebnisse
%     plot(T_M(:,p(m)),z_M,'Color',[0.8500 0.3250 0.0980],'LineWidth',1); %Matlab Ergebnisse
%     grid on
%     xticks([40, 50, 60, 70, 80])
%     %Achse auf Matlab Koordinaten reduzieren (Comsol hat mehr Höhencoordinaten wegen Halbkugel oben und unten)
%     axis([40 80 z_Mu z_Mo])
%     % Beschriftung nur am ersten Sublot
%     if m == 1
%         ylabel('System depth [m]','Interpreter','Latex','Fontsize',12)
%         xlabel('Temperatur in TPPS [$^{\circ}$C] ','Interpreter','Latex','Fontsize',12)
%     end
%     % Legende nur am letzten Subplot
%     if m == length(p)
%         legend('$T_{Sys}$ Comsol','$T_{Pist}$ Comsol','$T_{Sys}$ Matlab','Interpreter','Latex','Fontsize',10)
%     end
%     sgtitle(['1D Temperature Profil in the 72 m diameter System at Times: ' num2str(t_Sys_M(p)) ] ,'Interpreter','Latex','Fontsize',12)
% end
