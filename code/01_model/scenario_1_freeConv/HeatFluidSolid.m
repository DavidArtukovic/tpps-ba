%% HEATFLUIDSOLID.m
% ---------------------------------------------------------------
%  SUMMARY:
%   Computes the time derivatives of all temperature states in the TPPS
%   model during operation with flowing water in the storage.
%
% DESCRIPTION:
%   - Assembles 1D heat conduction in piston, water, air, insulation and
%     vertical soil.
%   - Adds 1D axial forced convection in the water column based on the flow flag
%     "flow".
%   - Assembles 2D heat conduction in the surrounding radial soil and
%     radial insulation.
%   - Couples the 1D system (piston / water / air / insulation) to the
%     2D radial soil via mixed (Robin-type) boundary conditions.
%   - Returns dTdt in a format compatible with ode45 for time integration.
%
% INPUT:
%   t      - time (s), required by ODE interface but not used explicitly
%   T      - temperature state vector (K), structured as:
%            [piston(1D); vertical soil(1D); water(1D); air(1D);
%             insulation(1D); radial soil/insulation(2D -> flattened)]
%   Nz     - vector with number of grid points in each sub-domain
%   dz     - vertical grid spacing (m)
%   bypass_indices - index of the lower and upper bypass inlet in the 1D system and membrane
%   flow   - signed water flow rate in the 1D water domain
%            (m/s or consistent internal units). The sign determines
%            charging / discharging direction.
%   T0init - boundary / ambient temperatures and reference values
%   SW     - material property matrix (density, heat capacity, etc.)
%   A      - cross-sectional areas of the different regions
%   z_RE   - geometric information for the radial grid (variable spacing)
%
% OUTPUT:
%   dTdt   - time derivative of the temperature state vector (K/s)
%
% NOTES:
%   Visually, the grid and positions can be inspected in the
%   Teilsysteme_TPPS.svg. The present function corresponds to the
%   configuration with flowing water in the storage.
%
% ---------------------------------------------------------------


function dTdt = HeatFluidSolid(t, T, Nz, dz, bypass_indices, flow, T0init, SW, A, z_RE, fklog)
    % Pre-allocate temperature time derivative vector:
    %  - Nz(10) = number of 1D vertical nodes in piston, soil, water, air,
    %             and 1D insulation
    %  - Nz(12)*Nz(13) = number of nodes in the 2D radial soil / insulation
    dTdt = zeros(Nz(10)+Nz(12)*Nz(13),1); 

    % Initialize radial soil temperature field (2D) and its time derivative.
    % rows    -> vertical direction
    % columns -> radial direction
    Tf = zeros(Nz(13),Nz(12)); % temperature field
    dTdtf = zeros(Nz(13),Nz(12)); % temperature time derivative field

    % Temperature difference vectors used for coupling terms
    DT = zeros(Nz(13),1); % soil - water temperature difference (radial)
    DTd = zeros(Nz(9),3); % insulation temperatures (1D / 2D coupling)

    % top of the 1D system (air node)
    sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4;

    %%%------------------------------------------%%%
    % 01 Vertical Insulation and Air (1D)
    %%%------------------------------------------%%%

    % Prescribed air temperature at the very top node
    T(sys_top_idx) = T0init(4);

    % Contact temperature between water and vertical insulation on the water side.
    % When charging (flow < 0), the water temperature is prescribed by the
    % charging inlet temperature; otherwise we take the local water node.
    if flow < 0
        TK_WD = T0init(2); % charging temperature at water / insulation
    else
        TK_WD = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-4); % water node next to insulation
    end

    % Contact temperature between insulation and water on the insulation side
    TK_DW = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2);

    % Mixed contact temperature water <-> vertical insulation (1D)
    T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3) = ...
        (SW(1,4)*TK_WD + SW(3,4)*TK_DW) / (SW(1,4) + SW(3,4));


    %%%------------------------------------------%%%
    % 02 Vertical soil below the water (1D) and water/soil contact
    %%%------------------------------------------%%%

    % Prescribed soil temperature at the very bottom of the vertical soil column
    T(Nz(3)+1) = T0init(5);

    % Contact temperature between soil and water at the bottom interface
    TK_EW = T(Nz(3)+Nz(8)-1);

    % Contact temperature between water and soil (water side).
    % When discharging (flow > 0) we impose an outlet/inlet temperature,
    % otherwise we use the local water node.
    if flow > 0
        TK_WE = T0init(3);
    else
        TK_WE = T(Nz(3)+Nz(8)+1);
    end

    % Mixed contact temperature water <-> vertical soil at the bottom
    T(Nz(3)+Nz(8)) = ...
        (SW(1,4)*TK_WE + SW(2,4)*TK_EW) / (SW(1,4) + SW(2,4));

    % 1D heat conduction in the vertical soil column (between bottom and water)
    for i = Nz(3)+2 : Nz(3)+Nz(8)-1
        % Discretized one dimensional heat equation
        dTdt(i) = SW(2,5) * (T(i+1) - 2*T(i) + T(i-1)) / dz^2;
    end

    %%%------------------------------------------%%%
    % 03 Piston (1D) including coupling to water
    %%%------------------------------------------%%%

    %--------------------------------------------------------------
    % 3.1 Piston top contact (to upper water region)
    %--------------------------------------------------------------

    % Water/piston contact temperature at the piston top (water side)
    if flow >= 0
        % thermal discharging / standstill
        TK_WKO = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)-2);
    else
        % thermal charging
        TK_WKO = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)-1);
    end
    % Piston/water contact temperature at piston top (piston side)
    TK_KWO = T(Nz(3)-1);

    % Mixed contact temperature at piston top
    T(Nz(3)) = (SW(1,4)*TK_WKO+SW(2,4)*TK_KWO)/((SW(1,4)+SW(2,4)));

    %--------------------------------------------------------------
    % 3.2 Piston bottom contact (to lower water region) 
    %--------------------------------------------------------------

    % Water/piston contact temperature at piston bottom (water side)
    if flow <= 0
        % thermal charging / standstill
        TK_WKU = T(Nz(3)+Nz(8)+Nz(2)-1);
    else 
        % thermal discharging
        TK_WKU = T(Nz(3)+Nz(8)+Nz(2)-2);
    end
    % Piston/water contact temperature at piston bottom (piston side)
    TK_KWU = T(2);

    % Mixed contact temperature at piston bottom
    T(1) = (SW(1,4)*TK_WKU+SW(2,4)*TK_KWU)/((SW(1,4)+SW(2,4)));

    % 1D heat conduction inside the piston
    for i = 2:Nz(3)-1
        dTdt(i) =  SW(2,5)*(T(i+1)-2*T(i)+T(i-1))/(dz^2); 
    end

    %%%------------------------------------------%%%
    % 04 Map 1D temperature vector to 2D radial soil field
    %%%------------------------------------------%%%

    % Temperaturvektor in Temperaturfeld umwandeln
    for j = 1:Nz(13)
        for i = 1:Nz(12)
            Tf(j,i) = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4+(i-1)*Nz(13)+j); 
        end
    end

    %%%------------------------------------------%%%
    % 05 Boundary conditions for radial soil / insulation (2D)
    %%%------------------------------------------%%%

    % Vertical boundaries of the radial soil
    Tf(1,:) = T0init(5); % bottom soil (same as vertical soil bottom)
    Tf(end,:) = T0init(4); % top soil equals air temperature

    % Outer radial boundary of the soil
    Tf(:,end) = T0init(5);

    %--------------------------------------------------------------
    % 5.1 Coupling 1D insulation with 2D radial insulation
    %--------------------------------------------------------------

    % Extract vertical temperature profile in 1D insulation at system top
    sys_top_no_ins_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3;
    DTd(:,1) = T(sys_top_no_ins_idx :sys_top_idx);

    % Indices for the top of the radial soil column (first radial ring)
    radial_soil_top_idx        = Nz(2)+Nz(3)+Nz(4)+Nz(9)-3;
    radial_soil_top_no_ins_idx = Nz(2)+Nz(3)+Nz(4)-2;

    % Vertical temperature profile in the second radial column (just inside insulation)
    DTd(:,2) = Tf(radial_soil_top_no_ins_idx:radial_soil_top_idx,2); 

    % Average contact temperature between 1D insulation and 2D radial insulation
    DTd(:,3) = (DTd(:,1)+DTd(:,2))/2;

    % Update contact temperatures in the first radial column of the soil / insulation
    Tf(Nz(2)+Nz(3)+Nz(4)-2:Nz(2)+Nz(3)+Nz(4)+Nz(9)-3,1) = DTd(:,3); 


    % Mixed contact between vertical soil and radial insulation at the very top
    Tf(end-Nz(9)+1,:) = ...
        (SW(2,4)*Tf(end-Nz(9),:) + SW(3,4)*Tf(end-Nz(9)+2,:)) / (SW(2,4)+SW(3,4));

    %--------------------------------------------------------------
    % 5.2 Coupling 1D water with 2D radial soil at the inner radius
    %--------------------------------------------------------------

    % Coupling 1D water (lower part) with 2D radial soil at inner radius
    for e = 1:Nz(2)
        Tf(e,1) = ...
            (SW(1,4)*T(Nz(3)+Nz(8)+e-1) + SW(2,4)*Tf(e,2)) / (SW(1,4)+SW(2,4));
    end

    % Coupling 1D water (upper part) with 2D radial soil at inner radius
    for e = Nz(2)+Nz(3)-1 : Nz(2)+Nz(3)+Nz(4)-2
        Tf(e,1) = ...
            (SW(1,4)*T(Nz(8)+Nz(5)+e-1) + SW(2,4)*Tf(e,2)) / (SW(1,4)+SW(2,4));
    end

    %%%------------------------------------------%%%
    % 06 Temperature gradients in radial soil / insulation (2D)
    %%%------------------------------------------%%%

    % 6.1 Soil region until radial insulation
    for j = 2 : Nz(13)-Nz(9)          % vertical index
        for i = 2 : Nz(12)-1          % radial index
            dTdtf(j,i) = SW(2,5) * ( ...
                (Tf(j+1,i) - 2*Tf(j,i) + Tf(j-1,i)) / dz^2 + ...
                (Tf(j,i+1) - Tf(j,i-1)) ./ (z_RE(1,i) .* z_RE(3,i)) + ...
                ( (Tf(j,i+1) - Tf(j,i)) ./ z_RE(4,i) + ...
                (Tf(j,i-1) - Tf(j,i)) ./ z_RE(5,i) ) .* 2 ./ z_RE(3,i) );
        end
    end

    % 6.2 Radial insulation region
    for j = Nz(13)-Nz(9)-2 : Nz(13)-1 % vertical index
        for i = 2 : Nz(12)-1          % radial index
            dTdtf(j,i) = SW(3,5) * ( ...
                (Tf(j+1,i) - 2*Tf(j,i) + Tf(j-1,i)) / dz^2 + ...
                (Tf(j,i+1) - Tf(j,i-1)) ./ (z_RE(1,i) .* z_RE(3,i)) + ...
                ( (Tf(j,i+1) - Tf(j,i)) ./ z_RE(4,i) + ...
                (Tf(j,i-1) - Tf(j,i)) ./ z_RE(5,i) ) .* 2 ./ z_RE(3,i) );
        end
    end

    % Map 2D radial temperature gradients back into the global dTdt vector
    for j = 1 : Nz(13)
        for i = 1 : Nz(12)
            dTdt(sys_top_idx + (i-1)*Nz(13) + j,1) = dTdtf(j,i);
        end
    end

    %%%------------------------------------------%%%
    % 07 Vertical insulation (1D) including radial losses
    %%%------------------------------------------%%%

    % radial temperature difference between first and second radial nodes
    DT(:,1) = Tf(:,2) - Tf(:,1);

    % Temperature gradient in the 1D system insulation including radial loss
    for i = sys_top_no_ins_idx+1 : sys_top_idx-1
        dTdt(i) = SW(3,5) * (T(i+1) - 2*T(i) + T(i-1)) / dz^2 + ...
                SW(4,1) * DT(i-Nz(8)-Nz(5)-1);
    end


    %%%------------------------------------------%%%
    % 08 Water region (1D) with axial conduction, convection and radial losses
    %%%------------------------------------------%%%

    idx_w_top    = sys_top_no_ins_idx - 1;      % topmost interior water node
    idx_w_bottom = Nz(3)+Nz(8) + 1;             % bottommost interior water node

    % Bypass indices for the current time step
    idx_bypass_lower = bypass_indices(1);
    idx_bypass_upper = bypass_indices(2);
    idx_membrane = bypass_indices(3);

    for i = Nz(3)+Nz(8)+1 : idx_w_top

        is_adv_segment = (i >= idx_bypass_upper && i <= idx_w_top) || ...
                        (i >= idx_w_bottom    && i <= idx_bypass_lower);

        if ~is_adv_segment || flow == 0
            adv = 0;
        elseif flow > 0
            adv = - flow * (T(i) - T(i-1)) / dz;
        else % flow < 0
            adv = - flow * (T(i+1) - T(i)) / dz;
        end
        
        if i <= Nz(3)+Nz(8)+Nz(2)-1
            % Water segment below< the piston
            radial_idx = i - (Nz(3)+Nz(8)-1);
            dTdt(i) = SW(1,5) * (T(i+1) - 2*T(i) + T(i-1)) / dz^2 ...
                    + adv ...
                    + SW(4,3) * DT(radial_idx);

        elseif i >= Nz(3)+Nz(8)+Nz(2)+Nz(5)-2
            % Water segment above the piston
            radial_idx = i - (Nz(5)+Nz(8)-1);
            dTdt(i) = SW(1,5) * (T(i+1) - 2*T(i) + T(i-1)) / dz^2 ...
                    + adv ...
                    + SW(4,3) * DT(radial_idx);
                
        else 
            % Now we are in the ring gap / piston region
            % Check wether we are at the membrane position
            if i == idx_membrane
                % Membrane node: do not exchange axially with neighbors (adiabatic)
                dTdt(i) = 0; 
                continue
            elseif i == idx_membrane - 1
                % work with ghost cells to remove coupling to i+1 (which would be membrane)
                diffusion = 2*(T(i-1) - T(i)) / dz^2;
            elseif i == idx_membrane + 1
                % Upper neighbor: remove coupling to i-1 (which would be membrane), work with ghost cells
                diffusion = 2*(T(i+1) - T(i)) / dz^2;
            else
                % Standard Laplacian
                diffusion = (T(i+1) - 2*T(i) + T(i-1)) / dz^2;
            end
            % Water segment inside the piston region (no direct radial coupling)
            % update with axial conduction + advection only
            % between bypass indices only difussion is active
            % at mebrane node dTdt = 0 (adiabatic)
            dTdt(i) = SW(1,5) * diffusion ...
                    + adv;
        end
    end

    %%%------------------------------------------%%%
    % 09 Additional piston loss terms into the water nodes
    %%%------------------------------------------%%%

    % Heat flux from piston top into upper water node
    alpha_pw = (SW(2,1) * A(2)) / (A(1) * SW(1,2) * SW(1,3) * dz^2);

    % Upper water node adjacent to piston
    idx_water_top_piston = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;
    dTdt(idx_water_top_piston) = dTdt(idx_water_top_piston) ...
        - alpha_pw * (T(Nz(3)) - T(Nz(3)-1));

    % Lower water node adjacent to piston
    idx_water_bottom_piston = Nz(3)+Nz(8)+Nz(2)-1;
    dTdt(idx_water_bottom_piston) = dTdt(idx_water_bottom_piston) ...
        - alpha_pw * (T(1) - T(2));
end