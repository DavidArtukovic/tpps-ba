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
%   z_RP   - geometric information for the radial grid in the piston region
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


function dTdt = HeatFluidSolid(t, T, Nz, dz, bypass_indices, flow, T0init, SW, A, z_RE, z_RP)
    % Pre-allocate temperature time derivative vector:
    %  - Nz(10) = number of 1D vertical nodes in piston, soil, water, air,
    %             and 1D insulation
    %  - Nz(12)*Nz(13) = number of nodes in the 2D radial soil / insulation
    %  - Nz(14)*Nz(3) = number of nodes in the 2D piston
    dTdt = zeros(Nz(10)+Nz(12)*Nz(13)+Nz(14)*Nz(3),1); 


    % indices in the 1D system vector for coupling with 2D fields
    sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4;

    soil_2d_first_idx = sys_top_idx+1;                                  % first index of 2D soil field in system vector
    soil_2d_top_idx = soil_2d_first_idx + Nz(12)*Nz(13) - 1;            % last index of 2D soil field in system vector

    piston_2d_first_idx = soil_2d_top_idx+1;                      % first index of 2D piston field in system vector
    piston_2d_top_idx = piston_2d_first_idx + Nz(14)*Nz(3) - 1;   % last index of 2D piston field in system vector


    % Mapping: 1D-temperature vector → 2D-temperature field for the radial soil 
    block_RE = T(soil_2d_first_idx:soil_2d_top_idx);
    T_REf    = reshape(block_RE, Nz(13), Nz(12));

    % Initialize radial soil temperature field (2D) and its time derivative.
    % rows    -> vertical direction
    % columns -> radial direction
    dTdt_REf = zeros(Nz(13), Nz(12)); % temperature time derivative field

    % Mapping: 1D-temperature vector → 2D-temperature field for piston 
    block_P = T(piston_2d_first_idx:piston_2d_top_idx);
    T_Pf    = reshape(block_P, Nz(3), Nz(14));

    % Initilization of temperature and gradient field in piston (P)
    dTdt_Pf = zeros(Nz(3), Nz(14)); % rows -> vertical, columns -> radial


    % Temperature difference vectors used for coupling terms
    DT = zeros(Nz(13),1); % soil - water temperature difference (radial)
    DTd = zeros(Nz(9),3); % insulation temperatures (1D / 2D coupling)

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
    % 03 Piston (2D) including coupling to water
    %%%------------------------------------------%%%

    % set legacy 1D piston temperature time derivatives to zero
    dTdt(1:Nz(3)) = 0;

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
    TK_KWO_f = T_Pf(end-1, :); % radial temperature distribution at piston top

    % Mixed contact temperature at piston top
    T_Pf(end,:) = (SW(1,4)*TK_WKO + SW(2,4)*TK_KWO_f) / (SW(1,4)+SW(2,4));

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
    TK_KWU_f = T_Pf(2, :); % radial temperature distribution at piston bottom

    % Mixed contact temperature at piston bottom
    T_Pf(1,:)   = (SW(1,4)*TK_WKU + SW(2,4)*TK_KWU_f) / (SW(1,4)+SW(2,4));


    %--------------------------------------------------------------
    % 3.3 Piston radial contact (to ring-gap water) 
    %--------------------------------------------------------------

        N_piston = Nz(3);
        N_ring = Nz(5);

        for e = 1:Nz(3)

            % relative position in soil column (0 ... 1)
            rel = (e - 1) / (N_piston-1);

            % corresponding ring-gap index (1 ... Nz(5))
            j = floor(rel * (N_ring-1)) + 1;

            T_Pf(e,end) = ...
                (SW(1,4) * T(Nz(3)+Nz(8)+Nz(2)-1 + j) ...
                + SW(2,4) * T_Pf(e,end-1)) ...
                / (SW(1,4) + SW(2,4));
        end



    %--------------------------------------------------------------
    % 3.4 Map 2D piston temperature field back to 1D system vector
    %--------------------------------------------------------------

    T(piston_2d_first_idx:piston_2d_top_idx) = T_Pf(:);

    %--------------------------------------------------------------
    % 3.5 Temperature gradients in the piston (2D) with coupling to water at top and bottom
    %--------------------------------------------------------------

    jIdx = 2:Nz(3)-1;      % interior z
    iIdx = 2:Nz(14)-1;     % interior r

    % vertical second derivative
    d2T_dz2 = (T_Pf(jIdx+1,iIdx) - 2*T_Pf(jIdx,iIdx) + T_Pf(jIdx-1,iIdx)) / dz^2;

    % radial geometry terms
    drsum = z_RP(3,iIdx);
    dr_f  = z_RP(4,iIdx);
    dr_b  = z_RP(5,iIdx);
    ri    = z_RP(1,iIdx);

    % radial second derivative (non-uniform cylindrical)
    d2T_dr2 = (2 ./ drsum) .* ( ...
                (T_Pf(jIdx,iIdx+1) - T_Pf(jIdx,iIdx)) ./ dr_f ...
              - (T_Pf(jIdx,iIdx)   - T_Pf(jIdx,iIdx-1)) ./ dr_b );

    % 1/r * dT/dr term
    dT_dr = (T_Pf(jIdx,iIdx+1) - T_Pf(jIdx,iIdx-1)) ./ drsum;
    one_over_r_dTdr = (1 ./ ri) .* dT_dr;

    % assemble
    dTdt_Pf(jIdx,iIdx) = SW(2,5) * (d2T_dz2 + d2T_dr2 + one_over_r_dTdr);
    
    %--------------------------------------------------------------
    % 3.6 Axis treatment (r = 0 symmetry): vectorized over z
    %--------------------------------------------------------------

    dr0 = z_RP(5,2);

    d2T_dz2_axis = (T_Pf(jIdx+1,1) - 2*T_Pf(jIdx,1) + T_Pf(jIdx-1,1)) / dz^2;
    radial_axis = 4 * (T_Pf(jIdx,2) - T_Pf(jIdx,1)) / dr0^2;

    dTdt_Pf(jIdx,1) = SW(2,5) * (d2T_dz2_axis + radial_axis);

    DT_p = T_Pf(:,end-1) - T_Pf(:,end); % gradient at piston outer radius (ring gap)

    %--------------------------------------------------------------
    % 3.7 set all boundary gradients to zero (for now)
    %--------------------------------------------------------------
    dTdt_Pf(:,end) = 0; % piston - ring-gap water
    dTdt_Pf(1,:) = 0;   % piston - lower water volume
    dTdt_Pf(end,:) = 0; % piston - upper water volume
 
    % Bring bakc in 1D system vector
    dTdt(piston_2d_first_idx: piston_2d_top_idx) = dTdt_Pf(:);

    %%%------------------------------------------%%%
    % 04 Map 1D temperature vector to 2D radial soil field
    %%%------------------------------------------%%%

    % Temperaturvektor in Temperaturfeld umwandeln
    for j = 1:Nz(13)
        for i = 1:Nz(12)
            T_REf(j,i) = T(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4+(i-1)*Nz(13)+j); 
        end
    end

    %%%------------------------------------------%%%
    % 05 Boundary conditions for radial soil / insulation (2D)
    %%%------------------------------------------%%%

    % Vertical boundaries of the radial soil
    T_REf(1,:) = T0init(5); % bottom soil (same as vertical soil bottom)
    T_REf(end,:) = T0init(4); % top soil equals air temperature

    % Outer radial boundary of the soil
    T_REf(:,end) = T0init(5);

    %--------------------------------------------------------------
    % 5.1 Coupling 1D insulation with 2D radial insulation
    %--------------------------------------------------------------

    % Extract vertical temperature profile in 1D insulation at system top
    sys_top_no_ins_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3;
    DTd(:,1) = T(sys_top_no_ins_idx:sys_top_idx);

    % Indices for the top of the radial soil column (first radial ring)
    radial_soil_top_idx        = Nz(2)+Nz(3)+Nz(4)+Nz(9)-3;
    radial_soil_top_no_ins_idx = Nz(2)+Nz(3)+Nz(4)-2;

    % Vertical temperature profile in the second radial column (just inside insulation)
    DTd(:,2) = T_REf(radial_soil_top_no_ins_idx:radial_soil_top_idx,2); 

    % Average contact temperature between 1D insulation and 2D radial insulation
    DTd(:,3) = (DTd(:,1)+DTd(:,2))/2;

    % Update contact temperatures in the first radial column of the soil / insulation
    T_REf(Nz(2)+Nz(3)+Nz(4)-2:Nz(2)+Nz(3)+Nz(4)+Nz(9)-3,1) = DTd(:,3); 


    % Mixed contact between vertical soil and radial insulation at the very top
    T_REf(end-Nz(9)+1,:) = ...
        (SW(2,4)*T_REf(end-Nz(9),:) + SW(3,4)*T_REf(end-Nz(9)+2,:)) / (SW(2,4)+SW(3,4));

    %--------------------------------------------------------------
    % 5.2 Coupling 1D water with 2D radial soil at the inner radius
    %--------------------------------------------------------------

    % Coupling 1D water (lower part) with 2D radial soil at inner radius
    for e = 1:Nz(2)
        T_REf(e,1) = ...
            (SW(1,4)*T(Nz(3)+Nz(8)+e-1) + SW(2,4)*T_REf(e,2)) / (SW(1,4)+SW(2,4));
    end

    % Coupling 1D water (upper part) with 2D radial soil at inner radius
    for e = Nz(2)+Nz(3)-1 : Nz(2)+Nz(3)+Nz(4)-2
        T_REf(e,1) = ...
            (SW(1,4)*T(Nz(8)+Nz(5)+e-1) + SW(2,4)*T_REf(e,2)) / (SW(1,4)+SW(2,4));
    end

    % Coupling 1D water (ring-gap) with 2D radial soil
    e_start = Nz(2);
    e_end   = Nz(2) + Nz(3) - 1;

    N_soil = Nz(3);
    N_ring = Nz(5);

    for e = e_start:e_end

        % relative position in soil column (0 ... 1)
        rel = (e - e_start) / (N_soil-1);

        % corresponding ring-gap index (1 ... Nz(5))
        j = floor(rel * (N_ring-1)) + 1;

        T_REf(e,1) = ...
            (SW(1,4) * T(Nz(3)+Nz(8)+Nz(2)-1 + j) ...
            + SW(2,4) * T_REf(e,2)) ...
            / (SW(1,4) + SW(2,4));
    end


    %%%------------------------------------------%%%
    % 06 Temperature gradients in radial soil / insulation (2D)
    %%%------------------------------------------%%%

    % 6.1 Soil region until radial insulation
    for j = 2 : Nz(13)-Nz(9)          % vertical index
        for i = 2 : Nz(12)-1          % radial index
            dTdt_REf(j,i) = SW(2,5) * ( ...
                (T_REf(j+1,i) - 2*T_REf(j,i) + T_REf(j-1,i)) / dz^2 + ...
                (T_REf(j,i+1) - T_REf(j,i-1)) ./ (z_RE(1,i) .* z_RE(3,i)) + ...
                ( (T_REf(j,i+1) - T_REf(j,i)) ./ z_RE(4,i) + ...
                (T_REf(j,i-1) - T_REf(j,i)) ./ z_RE(5,i) ) .* 2 ./ z_RE(3,i) );
        end
    end

    % 6.2 Radial insulation region
    for j = Nz(13)-Nz(9)-2 : Nz(13)-1 % vertical index
        for i = 2 : Nz(12)-1          % radial index
            dTdt_REf(j,i) = SW(3,5) * ( ...
                (T_REf(j+1,i) - 2*T_REf(j,i) + T_REf(j-1,i)) / dz^2 + ...
                (T_REf(j,i+1) - T_REf(j,i-1)) ./ (z_RE(1,i) .* z_RE(3,i)) + ...
                ( (T_REf(j,i+1) - T_REf(j,i)) ./ z_RE(4,i) + ...
                (T_REf(j,i-1) - T_REf(j,i)) ./ z_RE(5,i) ) .* 2 ./ z_RE(3,i) );
        end
    end

    % Map 2D radial temperature gradients back into the global dTdt vector
    for j = 1 : Nz(13)
        for i = 1 : Nz(12)
            dTdt(sys_top_idx + (i-1)*Nz(13) + j,1) = dTdt_REf(j,i);
        end
    end

    %%%------------------------------------------%%%
    % 07 Vertical insulation (1D) including radial losses
    %%%------------------------------------------%%%

    % radial temperature difference between first and second radial nodes
    DT(:,1) = T_REf(:,2) - T_REf(:,1);

    % Temperature gradient in the 1D system insulation including radial loss
    for i = sys_top_no_ins_idx+1 : sys_top_idx-1
        dTdt(i) = SW(3,5) * (T(i+1) - 2*T(i) + T(i-1)) / dz^2 + ...
                SW(4,1) * DT(i-Nz(8)-Nz(5)-1);
    end


    %%%------------------------------------------%%%
    % 08 Water region (1D) with axial conduction, convection and radial losses
    %%%------------------------------------------%%%

    idx_w_top    = sys_top_no_ins_idx - 1;      % topmost interior water node
    idx_w_bottom = Nz(3) + Nz(8);             % bottommost interior water node

    % Bypass indices for the current time step
    idx_bypass_lower = bypass_indices(1);
    idx_bypass_upper = bypass_indices(2);
    idx_membrane = bypass_indices(3);

    % begin of piston for radial soil
    e_start_soil = Nz(2);
    e_end_soil   = Nz(2)+Nz(3)-1;

    % begin of piston for radial soil
    e_start_piston = 1;
    e_end_piston   = Nz(3);

    for i = idx_w_bottom : idx_w_top

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
                diffusion = 0;
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

            % Ring gap: scale axial diffusivity due to replacement-volume compression
            s_geom      = A(3) / A(1);              % = H(5)/H(3) because A1*H5 = A3*H3
            alpha_ring  = SW(1,5) * s_geom^2;       % preserve diffusion time scale tau ~ L^2/alpha

            % ----------------------------------------------------------
            % Extra axial mixing below/above outlet (eddy diffusion)
            % ----------------------------------------------------------
            alpha_max = 3e-5;   % [m^2/s] tune
            L_mix = max(abs(idx_bypass_upper - idx_membrane) * dz, dz); % [m]

            % Define mixing zone: between membrane and upper bypass
            in_mix_zone = (i >= idx_membrane) && (i <= idx_bypass_upper);

            if in_mix_zone
                d  = abs(i - idx_bypass_upper) * dz;      % [m], 0 at stub
                alpha_turb = alpha_max * (1 - d/(2*L_mix))*s_geom^2; % strong near inlet
            else
                alpha_turb = 0;
            end

            alpha_ring_eff = alpha_ring + alpha_turb;

            % Compute mean radial gradient based on current node position in the ring gap
            % ring-gap local index (1 ... Nz(5))
            j = i - (Nz(3)+Nz(8)+Nz(2)-1);

            % determine corresponding soil index interval
            e_low_soil  = e_start_soil + floor((j-1)/Nz(5) * Nz(3));
            e_high_soil = e_start_soil + floor(j/Nz(5)     * Nz(3)) - 1;

            % determine corresponding piston index interval
            e_low_piston  = 1 + floor((j-1)/Nz(5) * Nz(3));
            e_high_piston = 1 + floor(j/Nz(5)     * Nz(3)) - 1;

            % indices in soil relevant for water node in ring-gap
            e_low_soil  = max(e_low_soil,  e_start_soil);
            e_high_soil = min(e_high_soil, e_end_soil);

            % mean radial gradient to soil
            mean_DT = mean( DT(e_low_soil:e_high_soil,1) );

            % indices in soil relevant for water node in ring-gap
            e_low_piston  = max(e_low_piston,  e_start_piston);
            e_high_piston = min(e_high_piston, e_end_piston);

            % mean radial gradient to piston (for nodes adjacent to piston)
            mean_DT_p = mean( DT_p(e_low_piston:e_high_piston) );

            % Water segment inside the piston region (no direct radial coupling)
            % update with advection + scaled diffusion and  radial loss to outer soil
            % at mebrane node dTdt = 0 (adiabatic)
            dTdt(i) = alpha_ring_eff * diffusion ...
                    + adv ...
                    + SW(4,4) * mean_DT ...          % soil coupling
                    + SW(4,5) * mean_DT_p;           % piston coupling
        end
    end

    %%%------------------------------------------%%%
    % 09 Additional piston loss terms into the water nodes (distributed)
    % For numerical smoothening
    %%%------------------------------------------%%%

    alpha_pw = (SW(2,1) * A(2)) / (A(1) * SW(1,2) * SW(1,3) * dz^2);

    % Water nodes energy is distributed to
    Nspread = 20;                  % tune
    ell     = 4*dz;                % [m] decay length scale (controls smoothness)

    % Exponential weights, normalized to sum=1
    k  = (0:Nspread-1).';
    wk = exp(-(k*dz)/ell);
    wk = wk / sum(wk);  % normalize weights to sum to 1

    % normalize by sum of ZR_P(1,:) to get proper weights
    gradient_weights = z_RP(2,:) / sum(z_RP(2,:));

    % weighted mean temperature gradient at piston top
    dT_WKO_mean = sum((T_Pf(end-1,:)-T_Pf(end,:)).* gradient_weights);
    
    % weighted mean temperature gradient at piston bottom
    dT_WKU_mean = sum((T_Pf(1,:)-T_Pf(2,:)).* gradient_weights);
    

    % Upper water nodes adjacent to piston
    idx_water_top_piston = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;
    for j = 0:Nspread-1
        idx = idx_water_top_piston + j;            
        dTdt(idx) = dTdt(idx) - wk(j+1) * alpha_pw * dT_WKO_mean;
    end

    % Lower water nodes adjacent to piston
    idx_water_bottom_piston = Nz(3)+Nz(8)+Nz(2)-1;
    for j = 0:Nspread-1
        idx = idx_water_bottom_piston - j;    
        dTdt(idx) = dTdt(idx) - wk(j+1) * alpha_pw * dT_WKU_mean;
    end

end