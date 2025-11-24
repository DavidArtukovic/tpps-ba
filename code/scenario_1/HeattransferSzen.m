function [T_Sys,T_REf] = HeattransferSzen(t2,IC_Sys,Nz,dz,flow,T0init,SW,A,z_RE,T_REf,Nt2)

[t2,T_Sys] = ode45(@HeatFluidSolid,t2,IC_Sys,[],Nz,dz,flow,T0init,SW,A,z_RE);
if flow == 0 % thermischer Stillstand
    % HUK Wasser/Kolben oben
    T_Sys(:,Nz(3)) = (SW(1,4)*T_Sys(:,Nz(3)+Nz(8)+Nz(2)+Nz(5)-2)+SW(2,4)*T_Sys(:,Nz(3)-1))/((SW(1,4)+SW(2,4)));
    % HUK Wasser/Kolben unten
    T_Sys(:,1) = (SW(1,4)*T_Sys(:,Nz(3)+Nz(8)+Nz(2)-1)+SW(2,4)*T_Sys(:,2))/((SW(1,4)+SW(2,4)));
elseif flow > 0 % thermisches Entladen
    % HUK Wasser/Kolben oben
    T_Sys(:,Nz(3)) = (SW(1,4)*T_Sys(:,Nz(3)+Nz(8)+Nz(2)+Nz(5)-2)+SW(2,4)*T_Sys(:,Nz(3)-1))/((SW(1,4)+SW(2,4)));
    % HUK Wasser/Kolben unten
    T_Sys(:,1) = (SW(1,4)*T_Sys(:,Nz(3)+Nz(8)+Nz(2)-2)+SW(2,4)*T_Sys(:,2))/((SW(1,4)+SW(2,4)));
            
else % thermisches Laden
    % HUK Wasser/Kolben oben
    T_Sys(:,Nz(3)) = (SW(1,4)*T_Sys(:,Nz(3)+Nz(8)+Nz(2)+Nz(5)-1)+SW(2,4)*T_Sys(:,Nz(3)-1))/((SW(1,4)+SW(2,4)));
    % HUK Wasser/Kolben unten
    T_Sys(:,1) = (SW(1,4)*T_Sys(:,Nz(3)+Nz(8)+Nz(2)-1)+SW(2,4)*T_Sys(:,2))/((SW(1,4)+SW(2,4)));
end

% Lufttemperatur
T_Sys(:,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4) = T0init(4);
% HUK Wasser/Dämmung oben
T_Sys(:,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3) = (SW(1,4)*T_Sys(:,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-4)+SW(2,4)*T_Sys(:,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2))/((SW(1,4)+SW(2,4)));
% HUK Wasser Erdreich
T_Sys(:,Nz(3)+Nz(8)) = (SW(1,4)*T_Sys(:,Nz(3)+Nz(8)+1)+SW(2,4)*T_Sys(:,Nz(3)+Nz(8)-1))/((SW(1,4)+SW(2,4)));


% RADIALES ERDREICH
for tt = 1:1+Nt2 % Zeitschritte Für Radialen Tensor
    for j = 1:Nz(13) % Ortsschritte für vertikales Erdreich
        for i = 1:Nz(12) % Ortsschritte für radiales Erdreich
            T_REf(j,i) = T_Sys(tt,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4+(i-1)*Nz(13)+j); % Erzeugen des Temperaturfeldes (z,r)
        end
    end
        % Recalculation der Randbedingungen 
        T_REf(1,:) = T0init(5); % Recalc Erdreich vertikal UNTEN
        T_REf(end,:) = T0init(4); % Recalc Luft /Erdreich OBEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        T_REf(:,end) = T0init(5); % Recalc Erdreich Radial

    for e = 1:Nz(2) % Randbedingung untere Druckzone
        T_REf(e,1) = (SW(1,4)*T_Sys(tt,Nz(3)+Nz(8)+e-1)+SW(2,4)*T_REf(e,2))/(SW(1,4)+SW(2,4));
    end
    for e = Nz(2)+Nz(3)-1:Nz(2)+Nz(3)+Nz(4)-2 % Randbedingung obere Druckzone
        T_REf(e,1) = (SW(1,4)*T_Sys(tt,Nz(8)+Nz(5)+e-1)+SW(2,4)*T_REf(e,2))/(SW(1,4)+SW(2,4));
    end

    DTd(:,1) = T_Sys(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3:Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4); % Temperatur Dämmung aus 1D
    DTd(:,2) = T_REf(Nz(2)+Nz(3)+Nz(4)-2:Nz(2)+Nz(3)+Nz(4)+Nz(9)-3,2); % Temperatur Dämmung zweiter knoten aus 2D
    DTd(:,3) = (DTd(:,1)+DTd(:,2))/2; % Kontakttemperatur Dämmung 1D und 2D
    T_REf(Nz(2)+Nz(3)+Nz(4)-2:Nz(2)+Nz(3)+Nz(4)+Nz(9)-3,1) = DTd(:,3); % Recalc Kontakttemperatur Dämmung 1D Dämmung 2D (Radial)
    T_REf(end-Nz(9)+1,:) = (SW(2,4)*T_REf(end-Nz(9),:)+SW(3,4)*T_REf(end-Nz(9)+2,:))/(SW(2,4)+SW(3,4)); % Recalc Kontakttemperatur Dämmung Erdreich

% Temperaturfeld in Temperaturvektor Umwandeln
    for j = 1:Nz(13)
        for i = 1:Nz(12)
            T_Sys(tt,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4+(i-1)*Nz(13)+j) = T_REf(j,i); 
        end
    end
end

end