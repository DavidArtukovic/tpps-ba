function dTdt = HeatFluidSolid(t,T,Nz,dz,flow,T0init,SW,A,z_RE) % STABILER Geschwindigkeitsbereich keine Korrekturwerte
dTdt = zeros(Nz(10)+Nz(12)*Nz(13),1); % Temperaturgradientenvektor 
Tf = zeros(Nz(13),Nz(12)); % Initalisierung des Temperaturfeldes radiales Erdreich
dTdtf = zeros(Nz(13),Nz(12)); % Initalisierung des Temperaturgradientenfeldes radiales Erdreich
DT = zeros(Nz(13),1); % Temperaturdifferenz zwischen Erdboden und Wasser
DTd = zeros(Nz(9),3); % Temperatur Dämmung

% DÄMMUNG
% Lufttemperatur Oben
T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4) = T0init(4);

% Temperatur Wasser/Dämmung LADETEMPERATUR
if flow < 0
    TK_WD = T0init(2);
else
    TK_WD = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-4);
end
% Temperatur Dämmung/Wasser
TK_DW = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2);

% Kontakttemperatur Wasser Dämmung
T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3) = (SW(1,4)*TK_WD+SW(3,4)*TK_DW)/((SW(1,4)+SW(3,4)));

% ERDREICH
% Erdtemperatur Unten
T(Nz(3)+1) = T0init(5);

% Kontakttemperatur Erdreich/Wasser
TK_EW = T(Nz(3)+Nz(8)-1);
% Kontakttemperatur Wasser/Erdreich ENTLADETEMPERATUR
if flow > 0
    TK_WE = T0init(3);
else
    TK_WE = T(Nz(3)+Nz(8)+1);
end

% Kontakttemperatur Unten
T(Nz(3)+Nz(8)) = (SW(1,4)*TK_WE+SW(2,4)*TK_EW)/((SW(1,4)+SW(2,4)));

%Erdboden
for i = Nz(3)+2:Nz(3)+Nz(8)-1                                                       
    dTdt(i) = SW(2,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2);
end

% KOLBEN
% Kolben Oben
if flow >= 0 % thermischer Entladen/Stillstand
    % Wasser/Kolben Oben
    TK_WKO = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)-2);
else % thermisches Laden
    % Wasser/Kolben Oben
    TK_WKO = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)-1);
end
% Kolben/Wasser Oben
TK_KWO = T(Nz(3)-1);
% HUK Kolben oben
T(Nz(3)) = (SW(1,4)*TK_WKO+SW(2,4)*TK_KWO)/((SW(1,4)+SW(2,4)));

% Kolben Unten
if flow <= 0 % thermisches Laden/Stillstand
% Wasser/Kolben Unten
    TK_WKU = T(Nz(3)+Nz(8)+Nz(2)-1);
else % thermisches Entladen
    TK_WKU = T(Nz(3)+Nz(8)+Nz(2)-2);
end
% Kolben/Wasser Unten
TK_KWU = T(2);
% HUK Kolben unten
T(1) = (SW(1,4)*TK_WKU+SW(2,4)*TK_KWU)/((SW(1,4)+SW(2,4)));

% Kolben Ortsdiskretisierung
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
%Tf(:,1) = T0init(6); % Wasser Erdreich Randtemperatur

% % Kontakttemperatur Kolben Erdreich
% for e = Nz(2):Nz(2)+Nz(3)-1
%     Tf(e,1) = (T(1+e-Nz(2))+Tf(e,2))/2;
% end

DTd(:,1) = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3:Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4); % Temperatur Dämmung aus 1D
DTd(:,2) = Tf(Nz(2)+Nz(3)+Nz(4)-2:Nz(2)+Nz(3)+Nz(4)+Nz(9)-3,2); % Temperatur Dämmung zweiter knoten aus 2D
DTd(:,3) = (DTd(:,1)+DTd(:,2))/2; % Kontakttemperatur Dämmung 1D und 2D
Tf(Nz(2)+Nz(3)+Nz(4)-2:Nz(2)+Nz(3)+Nz(4)+Nz(9)-3,1) = DTd(:,3); % % Kontakttemperatur Dämmung 1D Dämmung 2D (Radial)

Tf(end-Nz(9)+1,:) = (SW(2,4)*Tf(end-Nz(9),:)+SW(3,4)*Tf(end-Nz(9)+2,:))/(SW(2,4)+SW(3,4)); % Kontakttemperatur Dämmung Erdreich (vertical)

% Kontakttemperatur Wasser unten Erdreich
for e = 1:Nz(2)
    Tf(e,1) = (SW(1,4)*T(Nz(3)+Nz(8)+e-1)+SW(2,4)*Tf(e,2))/(SW(1,4)+SW(2,4));
end

% Kontakttemperatur Wasser oben Erdreich
for e = Nz(2)+Nz(3)-1:Nz(2)+Nz(3)+Nz(4)-2
    Tf(e,1) = (SW(1,4)*T(Nz(8)+Nz(5)+e-1)+SW(2,4)*Tf(e,2))/(SW(1,4)+SW(2,4));
end


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

% % Temperaturänderung über die Zeit im Dämmmaterial
for i = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2:Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-5
    dTdt(i) = SW(3,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2) + SW(4,1)*DT(i-Nz(8)-Nz(5)-1); % Dämmung inc. Radialer verlustterm
end

% Zeile 145 -150 Integration der Verlustterme, bisher SW(4,3)*DT Eintrag aus dem Erdreich. Hier kommt die Konvektion hinzu. 
% Vegleich Gerle und Micha.
% Wasser
for i = Nz(3)+Nz(8)+1:Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-4
    if i <= Nz(3)+Nz(8)+Nz(2)-1 % Wasser unterhalb des Kolbens
        dTdt(i) = SW(1,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2)-flow*(T(i+1)-T(i-1))/(2*dz) + SW(4,3)*DT(i-(Nz(3)+Nz(8)-1));
    elseif i >= Nz(3)+Nz(8)+Nz(2)+Nz(5)-2 % Wasser oberhalb des Kolbens
        dTdt(i) = SW(1,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2)-flow*(T(i+1)-T(i-1))/(2*dz) + SW(4,3)*DT(i-(Nz(5)+Nz(8)-1));
    else
        dTdt(i) = SW(1,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2)-flow*(T(i+1)-T(i-1))/(2*dz);   
    end
end

%Verlustterme Kolben
dTdt(Nz(3)+Nz(8)+Nz(2)+Nz(5)-2) = dTdt(Nz(3)+Nz(8)+Nz(2)+Nz(5)-2) - (SW(2,1)*A(2))/(A(1)*SW(1,2)*SW(1,3)*dz^2)*(T(Nz(3))-T(Nz(3)-1));
dTdt(Nz(3)+Nz(8)+Nz(2)-1) = dTdt(Nz(3)+Nz(8)+Nz(2)-1)  - (SW(2,1)*A(2))/(SW(1,2)*SW(1,3)*A(1)*dz^2)*(T(1)-T(2));

end