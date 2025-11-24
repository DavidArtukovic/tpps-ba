function dTdt = HeatSolid(t,T,Nz,dz,flow,T0init,SW,A,z_RE) % STABILER Geschwindigkeitsbereich keine Korrekturwerte
dTdt = zeros(Nz(10)+Nz(12)*Nz(13),1); % Temperaturgradientenvektor 
Tf = zeros(Nz(13),Nz(12)); % Initalisierung des Temperaturfeldes radiales Erdreich
dTdtf = zeros(Nz(13),Nz(12)); % Initalisierung des Temperaturgradientenfeldes radiales Erdreich
DTd = zeros(Nz(9),1); % Temperatur Dämmung

% DÄMMUNG
% Lufttemperatur Oben
T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4) = T0init(4);
% Temperatur Wasser/Dämmung LADETEMPERATUR

% Kontakttemperatur Wasser Dämmung
T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3) = T0init(6);

% ERDREICH
% Erdtemperatur Unten
T(Nz(3)+1) = T0init(5);

% Kontakttemperatur Unten
T(Nz(3)+Nz(8)) = T0init(6);

%Erdboden
for i = Nz(3)+2:Nz(3)+Nz(8)-1                                                       
    dTdt(i) = SW(2,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2);
end

% KOLBEN
% Kolben Oben
T(Nz(3)) = T0init(6);

% Kolben Unten
% HUK Kolben unten
T(1) = T0init(6);

% Kolben Ortsdiskretisierung mit radialen Erdreich
for i = 2:Nz(3)-1
    dTdt(i) =  SW(2,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2); 
end

% Temperaturvektor in Temperaturfeld umwandeln
for j = 1:Nz(13)
    for i = 1:Nz(12)
        Tf(j,i) = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4+(i-1)*Nz(13)+j); 
    end
end
% Randbedingungen für Radiales Erdreich festlegen
Tf(1,:) = T0init(5); % Erdreich vertikal UNTEN
Tf(end,:) = T0init(4); % Lufttemperatur OBEN
Tf(:,end) = T0init(5); % Erdreich radial Endtemperatur
Tf(:,1) = T0init(6); % Wasser Erdreich Randtemperatur

DTd(:,1) = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3:Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4); % Temperatur Dämmung aus 1D
DTd(:,2) = Tf(Nz(2)+Nz(3)+Nz(4)-2:Nz(2)+Nz(3)+Nz(4)+Nz(9)-3,2); % Temperatur Dämmung zweiter knoten aus 2D
DTd(:,3) = (DTd(:,1)+DTd(:,2))/2; % Kontakttemperatur Dämmung 1D und 2D
Tf(Nz(2)+Nz(3)+Nz(4)-2:Nz(2)+Nz(3)+Nz(4)+Nz(9)-3,1) = DTd(:,3); % % Kontakttemperatur Dämmung 1D Dämmung 2D (Radial)
    
Tf(end-Nz(9)+1,:) = (SW(2,4)*Tf(end-Nz(9),:)+SW(3,4)*Tf(end-Nz(9)+2,:))/(SW(2,4)+SW(3,4)); % Kontakttemperatur Dämmung Erdreich
    
% Temperaturgradienten für radiales Erdreich bis Dämmung berechnenen 
for j = 2:Nz(13)-Nz(9) % Vertikale Komponenten
    for i = 2:Nz(12)-1 % Radiale Komponente
        %dTdtf(j,i) = SW(2,5)*(1/z_RE(1,i)*(Tf(j,i+1)-Tf(j,i-1))/z_RE(3,i)+((Tf(j,i+1)-Tf(j,i))/z_RE(4,i)+(Tf(j,i-1)-Tf(j,i))/z_RE(5,i))*2/z_RE(3,i)); % Ohne Vertikale Komponente
        dTdtf(j,i) = SW(2,5)*((Tf(j+1,i)-2*Tf(j,i)+Tf(j-1,i))/(dz^2)+(1/z_RE(1,i)*(Tf(j,i+1)-Tf(j,i-1))/z_RE(3,i)+((Tf(j,i+1)-Tf(j,i))/z_RE(4,i)+(Tf(j,i-1)-Tf(j,i))/z_RE(5,i))*2/z_RE(3,i))); % Mit Vertikale Komponente
    end
end

% Temperaturgradienten für radiale Dämmung berechnenen 
for j = Nz(13)-Nz(9)-2:Nz(13)-1 % Vertikale Komponenten
    for i = 2:Nz(12)-1 % Radiale Komponente
            %dTdtf(j,i) = SW(2,5)*(1/z_RE(1,i)*(Tf(j,i+1)-Tf(j,i-1))/z_RE(3,i)+((Tf(j,i+1)-Tf(j,i))/z_RE(4,i)+(Tf(j,i-1)-Tf(j,i))/z_RE(5,i))*2/z_RE(3,i)); % Ohne Vertikale Komponente
        dTdtf(j,i) = SW(3,5)*((Tf(j+1,i)-2*Tf(j,i)+Tf(j-1,i))/(dz^2)+(1/z_RE(1,i)*(Tf(j,i+1)-Tf(j,i-1))/z_RE(3,i)+((Tf(j,i+1)-Tf(j,i))/z_RE(4,i)+(Tf(j,i-1)-Tf(j,i))/z_RE(5,i))*2/z_RE(3,i))); % Mit Vertikale Komponente
    end
end

% Temperaturgradientenfeld in Temperaturgradientenvektor umwandeln
for j = 1:Nz(13)
    for i = 1:Nz(12)
        dTdt(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4+(i-1)*Nz(13)+j,1) = dTdtf(j,i); 
    end
end

DT(:,1)=Tf(:,2)-Tf(:,1);

% Temperaturänderung über die Zeit im Dämmmaterial
for i = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2:Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-5
    dTdt(i) = SW(3,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2) + SW(4,1)*DT(i-Nz(8)-Nz(5)-1); % Dämmung inc. Radialer verlustterm
end

end