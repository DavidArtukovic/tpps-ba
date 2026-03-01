%% HEATSOLID.m
% ---------------------------------------------------------------
%  SUMMARY:
%   Computes the time derivatives of all temperature states in the TPPS
%   model during the purely conductive initialization phase.
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
%   T0init - boundary / ambient temperatures and reference values
%   SW     - material property matrix (density, heat capacity, etc.)
%   z_RE   - radial grid spacing (m) for the radial soil discretization (variable spacing)
%   z_RP   - radial grid spacing (m) for the radial piston discretization (variable spacing)
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

function dTdt = HeatSolid(t, T, Nz, dz, T0init, SW, z_RE, z_RP) 

    % Flaten temperature gradient vector, vertical system + vertical soil system x radial soil system
    dTdt = zeros(Nz(10)+Nz(12)*Nz(13)+Nz(14)*Nz(3),1); 

    % Initilization of temperature gradient field in radial soil (RE)
    dTdt_REf = zeros(Nz(13),Nz(12)); % rows -> vertical, columns -> radial

    % Initilization of temperature gradient field in piston (P)
    dTdt_Pf = zeros(Nz(3),Nz(14)); % rows -> vertical, columns -> radial
    
    DTd = zeros(Nz(9),3); % temperature isolation

    %%%------------------------------------------%%%
    % 01 Isolation Dirichlet boundary conditions in 1D system
    %%%------------------------------------------%%%

    % air temperature top
    T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4) = T0init(4);

    % contact temperature water - isolation
    T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3) = T0init(6);

    %%%------------------------------------------%%%
    % 02 1D Axial Soil
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
    % 03 2D-Piston
    %%%------------------------------------------%%%
    
    % Piston temperature at piston top
    T(Nz(3)) = T0init(6);

    % Piston temperature at the first piston node (bottom)
    T(1) = T0init(6);

    % First index of 2D piston block
    sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4;      % system top index
    soil_2d_top_idx = sys_top_idx+Nz(12)*Nz(13);              % last index of 2D soil field in system vector
    piston_2d_first_idx = soil_2d_top_idx+1;                  % first index of 2D piston field in system vector

    % Extract 2D piston temperature block
    block_P = T(piston_2d_first_idx + (0:Nz(14)*Nz(3)-1));
    T_Pf    = reshape(block_P, Nz(3), Nz(14));

    % Dirichlet top and bottom
    T_Pf(1,:)   = T0init(6);
    T_Pf(end,:) = T0init(6);

    % Outer radial boundary
    T_Pf(:,end) = T0init(6);

    %--------------------------------------------------------------
    % 3.1 set gradients for 2D piston temperature field inside piston (vectorized)
    %--------------------------------------------------------------
    jIdx = 2:Nz(3)-1;      % interior z
    iIdx = 2:Nz(14)-1;     % interior r

    % vertical second derivative (central diff) for all interior nodes
    d2T_dz2 = (T_Pf(jIdx+1, iIdx) - 2*T_Pf(jIdx, iIdx) + T_Pf(jIdx-1, iIdx)) / (dz^2);

    % radial second derivative with variable spacing
    drsum = z_RP(3, iIdx);     % 1 × NrInt
    dr_f  = z_RP(4, iIdx);     % 1 × NrInt
    dr_b  = z_RP(5, iIdx);     % 1 × NrInt
    ri    = z_RP(1, iIdx);     % 1 × NrInt

    % Broadcast row vectors across all z rows via implicit expansion
    d2T_dr2 = (2 ./ drsum) .* ( ...
                (T_Pf(jIdx, iIdx+1) - T_Pf(jIdx, iIdx)) ./ dr_f  ...  % forward component
              - (T_Pf(jIdx, iIdx)   - T_Pf(jIdx, iIdx-1)) ./ dr_b );  % backward component

    % first radial derivative term (1/r * dT/dr) with variable spacing
    dT_dr = (T_Pf(jIdx, iIdx+1) - T_Pf(jIdx, iIdx-1)) ./ drsum;
    one_over_r_dTdr = (1 ./ ri) .* dT_dr;

    % assemble togehter in the interior of the piston domain
    dTdt_Pf(jIdx, iIdx) = SW(2,5) * ( d2T_dz2 + d2T_dr2 + one_over_r_dTdr );

    %--------------------------------------------------------------
    % 3.2 Axis treatment (r = 0 symmetry): vectorized over z
    %--------------------------------------------------------------
    dr0 = z_RP(5,2);  % distance between radial node 1 and 2

    % second derivative at the axis (central difference)
    d2T_dz2_axis = (T_Pf(jIdx+1,1) - 2*T_Pf(jIdx,1) + T_Pf(jIdx-1,1)) / (dz^2);      

    % radial second derivative at the axis (using ghost node approach for Neumann BC)
    radial_axis  = 4 * (T_Pf(jIdx,2) - T_Pf(jIdx,1)) / (dr0^2);  % = 2*(2*(...)/dr^2)   

    dTdt_Pf(jIdx,1) = SW(2,5) * ( d2T_dz2_axis + radial_axis ); % 1/r * dT/dr term drops out at the axis due to symmetry

    %--------------------------------------------------------------
    % 3.3 Dirichlet boundaries lead to zero temperature gradient 
    % (dT/dt = 0) at the boundary nodes, enforce this explicitly
    %--------------------------------------------------------------
    dTdt_Pf(1,:)   = 0;
    dTdt_Pf(end,:) = 0;
    dTdt_Pf(:,end) = 0;

    %--------------------------------------------------------------
    % 3.4  map 2D piston temperature gradient field into 1D system vector
    %--------------------------------------------------------------
    dTdt(piston_2d_first_idx + (0:Nz(14)*Nz(3)-1)) = dTdt_Pf(:);


    %%%------------------------------------------%%%
    % 04 2D-Radial Soil
    %%%------------------------------------------%%%

    sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4;      % system top index
    soil_2d_first_idx = sys_top_idx+1;                        % first index of 2D soil field in system vector

    % Mapping: 1D-temperature vector → 2D-temperature field for the radial soil 
    block_RE = T(soil_2d_first_idx + (0:Nz(12)*Nz(13)-1));
    T_REf    = reshape(block_RE, Nz(13), Nz(12));

    %--------------------------------------------------------------
    % 4.1 Boundary Conditions for Radial Soil
    %--------------------------------------------------------------
    T_REf(1,:) = T0init(5); % Soil temperature at the bottom
    T_REf(end,:) = T0init(4); % Soil temperature at the top equals air temperature
    T_REf(:,end) = T0init(5); % Radial soil (last column) outer temperature
    T_REf(:,1) = T0init(6); % Radial soil (first column -> contact to water) inner temperature

    %--------------------------------------------------------------
    % 4.2 Coupling 1D Insulation - 2D Soil
    %--------------------------------------------------------------
    sys_top_no_ins_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3; % system top without insulation
    DTd(:,1) = T(sys_top_no_ins_idx:sys_top_idx);           % temperature insulation 1D

    radial_soil_top_idx = Nz(2)+Nz(3)+Nz(4)+Nz(9)-3; % top index of radial soil (very left graphic in Teilsystem_TPPS.svg)
    radial_soil_top_no_ins_idx = Nz(2)+Nz(3)+Nz(4)-2; % top index of radial soil without insulattion (very left graphic in Teilsystem_TPPS.svg)
    DTd(:,2) = T_REf(radial_soil_top_no_ins_idx:radial_soil_top_idx, 2); % insulation temperature (above radial soil) of second radial grid elements

    DTd(:,3) = (DTd(:,1)+DTd(:,2))/2; % contact temperature between system insulation (1D) and radial soil insulation (2D)

    % update contact temperature in vertical direction between system insulation (1D) and radial insulation (2D)
    T_REf(radial_soil_top_no_ins_idx:radial_soil_top_idx,1) = DTd(:,3);
        
    % Update the contact temperature between the soil and insulation (2D soil grid)
    T_REf(end-Nz(9)+1,:) = (SW(2,4)*T_REf(end-Nz(9),:)+SW(3,4)*T_REf(end-Nz(9)+2,:))/(SW(2,4)+SW(3,4));
        
    %--------------------------------------------------------------
    % 4.3 set gradients for 2D radial soil temperature field (vectorized)
    %--------------------------------------------------------------

    % Shared index sets (interior in r, split in z for soil vs insulation)
    iIdx = 2:Nz(12)-1;                   % interior radial nodes
    jSoil = 2:(Nz(13)-Nz(9));            % vertical nodes in soil region (below insulation)
    jInsu = (Nz(13)-Nz(9)-2):(Nz(13)-1); % vertical nodes in insulation region

    % Shared radial geometry terms (row vectors, broadcast over z)
    drsum = z_RE(3, iIdx);            % [m] Δr_{i+1/2} + Δr_{i-1/2}
    dr_f  = z_RE(4, iIdx);            % [m] Δr_{i+1/2}
    dr_b  = z_RE(5, iIdx);            % [m] Δr_{i-1/2}
    ri    = z_RE(1, iIdx);            % [m] local radius r_i

    % Helper: radial operator (cylindrical, non-uniform grid)
    % Note: uses implicit expansion: (NzZ × NrInt) .* (1 × NrInt)
    radial_term = @(Tfld, jIdx) (...
        (1 ./ ri) .* (Tfld(jIdx, iIdx+1) - Tfld(jIdx, iIdx-1)) ./ drsum ... % 1/r * dT/dr term
      + (2 ./ drsum) .* ( ...
          (Tfld(jIdx, iIdx+1) - Tfld(jIdx, iIdx)) ./ dr_f ...   % forward diff component of d2T/dr2
        - (Tfld(jIdx, iIdx) - Tfld(jIdx, iIdx-1)) ./ dr_b ...   % backward diff component of d2T/dr2
        ) ...
    );

    % Soil region (material SW(2,*))
    d2T_dz2_soil = (T_REf(jSoil+1, iIdx) - 2*T_REf(jSoil, iIdx) + T_REf(jSoil-1, iIdx)) / (dz^2);
    dTdt_REf(jSoil, iIdx) = SW(2,5) * ( d2T_dz2_soil + radial_term(T_REf, jSoil) );

    % Insulation region (material SW(3,*))
    d2T_dz2_insu = (T_REf(jInsu+1, iIdx) - 2*T_REf(jInsu, iIdx) + T_REf(jInsu-1, iIdx)) / (dz^2);
    dTdt_REf(jInsu, iIdx) = SW(3,5) * ( d2T_dz2_insu + radial_term(T_REf, jInsu) );

    %--------------------------------------------------------------
    % 4.4 Dirichlet boundaries in 2D radial soil (keep fixed)
    %--------------------------------------------------------------
    dTdt_REf(1,:)   = 0;              % bottom boundary
    dTdt_REf(end,:) = 0;              % top boundary
    dTdt_REf(:,1)   = 0;              % inner radial boundary (contact column)
    dTdt_REf(:,end) = 0;              % outer radial boundary

    %--------------------------------------------------------------
    % 4.5  map 2D piston temperature gradient field into 1D system vector
    %--------------------------------------------------------------
    dTdt(soil_2d_first_idx + (0:Nz(12)*Nz(13)-1)) = dTdt_REf(:);


    %%%------------------------------------------%%%
    % 05 1D System Insulation
    %%%------------------------------------------%%%

    % radial temperature difference between first node and second node
    DT(:,1)=dTdt_REf(:,2)-dTdt_REf(:,1);

    for i = sys_top_no_ins_idx+1:sys_top_idx-1
        dTdt(i) = SW(3,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2)...
        + SW(4,1)*DT(i-Nz(8)-Nz(5)-1); % Isolation inc. radial loss term
    end

end