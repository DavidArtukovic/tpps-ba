%% HEATSOLID.m
% ---------------------------------------------------------------
%  SUMMARY:
%   Computes the time derivatives of all temperature states in the TPPS
%   model during the purely conductive initialization phase (no forced
%   water flow through the storage).
%
% DESCRIPTION:
%   - Assembles 1D heat conduction in piston, water, air, insulation and
%     vertical soil.
%   - Discretizes the 1D heat conduction equation in (Häsulein, 2024) with
%     the second-order central difference (for the second derivatives).
%   - Applies the prescribed boundary conditions at piston top and bottom,
%     at the vertical soil boundaries and at the outer radial boundary.
%   - Returns dTdt in a format compatible with ode45 for time integration.
%
% INPUT:
%   t      - current simulation time (s), required by ode45
%   T      - current temperature state vector of the whole system
%   Nz     - number of grid cells per 1D / 2D region
%   dz     - vertical grid spacing (m)
%   flow   - flow flag (kept for interface consistency, no convection here)
%   T0init - boundary / ambient temperatures and reference values
%   SW     - material property matrix (density, heat capacity, etc.)
%   A      - cross-sectional areas of the different regions
%   z_RE   - radial grid spacing (m) for the radial soil discretization (variable spacing)
%
% OUTPUT:
%   dTdt   - time derivative of the temperature state vector (K/s)
%
% NOTES:
%   Visually, the grid and positions can be inspected in the Teilsysteme_TPPS.svg 
%   at the left, the second from the left side (with orange soil) is the accurate 
%   graphic for this function. 
%
% ---------------------------------------------------------------

function dTdt = HeatSolid(t,T,Nz,dz,flow,T0init,SW,A,z_RE) 

% Flaten temperature gradient vector, vertical system + vertical soil system x radial soil system
dTdt = zeros(Nz(10)+Nz(12)*Nz(13),1); 

% Initilization of temperature (gradient) field in radial soil
Tf = zeros(Nz(13),Nz(12)); % rows -> vertical, columns -> radial
dTdtf = zeros(Nz(13),Nz(12)); % Initilization of temperature gradient field in radial soil

DTd = zeros(Nz(9),3); % temperature isolation

%%%------------------------------------------%%%
% 01 Isolation
%%%------------------------------------------%%%

% airtemperature top
T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4) = T0init(4);
% Temperatur Wasser/Dämmung LADETEMPERATUR

% Kontakttemperatur Wasser Dämmung
T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3) = T0init(6);

%%%------------------------------------------%%%
% 02 Soil
%%%------------------------------------------%%%

% Soil temperature at the first node (bottom)
T(Nz(3)+1) = T0init(5);

% Contact temperature between soil and piston
T(Nz(3)+Nz(8)) = T0init(6);

% discretized temperature gradient in the soil
for i = Nz(3)+2:Nz(3)+Nz(8)-1      
    % Discretized one dimensional heat equation (see eq.1 Häuslein 2024)                                                 
    dTdt(i) = SW(2,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2);
end

%%%------------------------------------------%%%
% 03 Piston
%%%------------------------------------------%%%

% Piston temperature at piston top
T(Nz(3)) = T0init(6);

% Piston temperature at the first piston node (bottom)
T(1) = T0init(6);

% discretized temperature gradient in the piston
for i = 2:Nz(3)-1
    % Discretized one dimensional heat equation (see eq.1 Häuslein 2024)   
    dTdt(i) =  SW(2,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2); 
end

% Mapping: 1D-temperature vector → 2D-temperature field for the radial soil 
sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4; % system top index
for j = 1:Nz(13)
    for i = 1:Nz(12)
        Tf(j,i) = T(sys_top_idx+(i-1)*Nz(13)+j); 
    end
end

%%%------------------------------------------%%%
% 04 Boundary Conditions for Radial Soil
%%%------------------------------------------%%%

Tf(1,:) = T0init(5); % Soil temperature at the bottom
Tf(end,:) = T0init(4); % Soil temperature at the top equals air temperature
Tf(:,end) = T0init(5); % Radial soil (last column) outer temperature
Tf(:,1) = T0init(6); % Radial soil (first column -> contact to water) inner temperature

%%%------------------------------------------%%%
% 05 Coupling 1D Insulation - 2D Soil
%%%------------------------------------------%%%
sys_top_no_ins_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3; % system top without insulation
DTd(:,1) = T(sys_top_no_ins_idx:sys_top_idx); % temperature insulation 1D

radial_soil_top_idx = Nz(2)+Nz(3)+Nz(4)+Nz(9)-3; % top index of radial soil (very left graphic in Teilsystem_TPPS.svg)
radial_soil_top_no_ins_idx = Nz(2)+Nz(3)+Nz(4)-2; % top index of radial soil without insulattion (very left graphic in Teilsystem_TPPS.svg)
DTd(:,2) = Tf(radial_soil_top_no_ins_idx:radial_soil_top_idx, 2); % insulation temperature (above radial soil) of second radial grid elements

DTd(:,3) = (DTd(:,1)+DTd(:,2))/2; % contact temperature between system insulation (1D) and radial soil insulation (2D)

% update contact temperature in vertical direction between system insulation (1D) and radial insulation (2D)
Tf(radial_soil_top_no_ins_idx:radial_soil_top_idx,1) = DTd(:,3);
    
% Update the contact temperature between the soil and insulation (2D soil grid)
Tf(end-Nz(9)+1,:) = (SW(2,4)*Tf(end-Nz(9),:)+SW(3,4)*Tf(end-Nz(9)+2,:))/(SW(2,4)+SW(3,4));
    
%%%------------------------------------------%%%
% 06 Temperature Gradients
%%%------------------------------------------%%%

% temperature gradient for radial soil until insulation (2D) 
for j = 2:Nz(13)-Nz(9) % vertical component
    for i = 2:Nz(12)-1 % radial component
        %dTdtf(j,i) = SW(2,5)*(1/z_RE(1,i)*(Tf(j,i+1)-Tf(j,i-1))/z_RE(3,i)+((Tf(j,i+1)-Tf(j,i))/z_RE(4,i)+(Tf(j,i-1)-Tf(j,i))/z_RE(5,i))*2/z_RE(3,i)); % Ohne Vertikale Komponente
        dTdtf(j,i) = SW(2,5)*((Tf(j+1,i)-2*Tf(j,i)+Tf(j-1,i))/(dz^2)...
        + (1/z_RE(1,i)*(Tf(j,i+1)-Tf(j,i-1))/z_RE(3,i)...
        + ((Tf(j,i+1)-Tf(j,i))/z_RE(4,i)...
        + (Tf(j,i-1)-Tf(j,i))/z_RE(5,i))*2/z_RE(3,i))); % with vertical component
    end
end

% temperature gradient for radial insulation (2D)
for j = Nz(13)-Nz(9)-2:Nz(13)-1 % vertical component
    for i = 2:Nz(12)-1 % radial component
            %dTdtf(j,i) = SW(2,5)*(1/z_RE(1,i)*(Tf(j,i+1)-Tf(j,i-1))/z_RE(3,i)+((Tf(j,i+1)-Tf(j,i))/z_RE(4,i)+(Tf(j,i-1)-Tf(j,i))/z_RE(5,i))*2/z_RE(3,i)); % Ohne Vertikale Komponente
        dTdtf(j,i) = SW(3,5)*((Tf(j+1,i)-2*Tf(j,i)+Tf(j-1,i))/(dz^2)...
        + (1/z_RE(1,i)*(Tf(j,i+1)-Tf(j,i-1))/z_RE(3,i)...
        + ((Tf(j,i+1)-Tf(j,i))/z_RE(4,i)...
        + (Tf(j,i-1)-Tf(j,i))/z_RE(5,i))*2/z_RE(3,i))); % with vertical component
    end
end

% Convert temperature gradient field into temperature gradient vector (2D -> 1D)
for j = 1:Nz(13)
    for i = 1:Nz(12)
        dTdt(sys_top_idx+(i-1)*Nz(13)+j,1) = dTdtf(j,i); 
    end
end

% radial temperature difference between first node and second node
DT(:,1)=Tf(:,2)-Tf(:,1);

% temperature gradient in system insulation (1D)
for i = sys_top_no_ins_idx+1:sys_top_idx-1 % Note the +1/-1 to allow for central difference calculation
    dTdt(i) = SW(3,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2)...
    + SW(4,1)*DT(i-Nz(8)-Nz(5)-1); % Isolation inc. radial loss term
end

end