%% HEATTRANSFERSZEN.m
% ---------------------------------------------------------------
%  SUMMARY:
%   Computes the transient TPPS system evolution for a given flow
%   scenario (charging, discharging or standstill) and couples the
%   1D system model to the 2D radial soil model.
%
% DESCRIPTION:
%   - Calls HEATFLUIDSOLID.m to integrate the full 1D TPPS model in
%     piston, water, air, insulation and vertical soil.
%   - Applies Robin-type interface conditions between water and piston
%     at top and bottom, depending on the flow direction.
%   - Enforces air temperature at the top of the system.
%   - Applies Robin-type coupling between water and vertical soil at the
%     water/soil interface.
%   - Maps the 1D temperature vector of the radial soil into a 2D
%     (z,r)-temperature field for each stored time step.
%   - Re-applies boundary conditions in the 2D radial soil domain,
%     including coupling to water at the inner radial boundary.
%   - Couples the 1D insulation temperatures to the 2D radial insulation
%     by averaging both sides at the contact interface.
%   - Maps the updated 2D temperature field back into the 1D system
%     temperature vector at every stored time step.
%
% INPUT:
%   t2     - time vector for the scenario simulation
%   IC_Sys - initial condition vector for the 1D system states
%   Nz     - vector with number of grid points for each subdomain
%   dz     - vertical grid spacing for the 1D domains
%   bypass_indices - indeces of the lower/upper bypass inlet and mebrane in the 1D system
%   flow   - flow flag:
%            = 0  : no flow (thermal standstill)
%            > 0  : discharging
%            < 0  : charging
%   T0init - vector of prescribed boundary / reference temperatures
%   SW     - parameter matrix with material and model switches
%   A      - cross-sectional areas and geometric parameters
%   z_RE   - geometric parameters of the radial soil grid
%   z_RP   - geometric parameters of the radial piston grid
%   T_REf  - 2D temperature field in the radial soil (z,r)
%   Nt2    - number of stored time steps during the scenario
%
% OUTPUT:
%   T_Sys  - time history of the 1D system temperature vector
%            Each row in the solution array T_Sys corresponds to 
%            a value returned in column vector t2: 
%            T_Sys  ∈  ℝ^(Nt × Nstates)
%   T_REf  - updated 2D radial soil temperature field (z,r)
%
% NOTES:
%   - Indexing conventions follow HEATTRANSFERINIT.m to keep the
%     coupling logic between 1D and 2D domains consistent.
% ---------------------------------------------------------------


function [T_Sys, T_REf, fk_code] = HeattransferSzen(t2, IC_Sys, Nz, dz, bypass_indices, flow,T0init,SW,A,z_RE,z_RP,T_REf,Nt2,fklog)

    %--------------------------------------------------------------
    % Preallocation for contact temperatures between system insulation (1D)
    % and radial insulation (2D)
    %--------------------------------------------------------------
    DTd = zeros(Nz(9),3);  % [1D insulation, 2D insulation, contact]

    %%%------------------------------------------%%%
    % 01: Build effective idle diffusivity profile and integrate ODE
    %%%------------------------------------------%%%

    alpha_w_eff = SW(1,5) * ones(numel(IC_Sys),1);

    if flow == 0
        % Empirical idle-mixing parameters for inverse thermocline treatment
        C_fk     = 2;                   % tune against COMSOL
        k_fk     = 0.5;                 % literature-inspired exponent
        grad_tol = .2;                  % [K/m] activation threshold
        alpha_cap = 10 * SW(1,5);       % upper bound for effective diffusivity

        alpha_w_eff = build_idle_diffusivity_profile(IC_Sys, Nz, dz, SW, C_fk, k_fk, grad_tol, alpha_cap);
    end

    [t2, T_Sys] = ode45(@HeatFluidSolid, t2, IC_Sys, [], Nz, dz, bypass_indices, flow, T0init, SW, A, z_RE, z_RP, alpha_w_eff);


    %%%------------------------------------------%%%
    % 02: Free convection as DISCRETE mixing operator (post-ODE)
    %%%------------------------------------------%%%

    % dt_mix = t2(2) - t2(1); % time step size for mixing operator [s]
    dt_mix = 900; % time step size for mixing operator [s]

    [T_Sys, fk_code] = apply_free_convection(T_Sys, dt_mix, flow, Nz, dz, SW, A, bypass_indices, fklog);
    
    %%%------------------------------------------%%%
    % 03: Define frequently used indices
    %%%------------------------------------------%%%

    % Air node at the very top of the system
    sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4;

    % Water / insulation interface at the system top (1D)
    sys_top_no_ins_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3;

    % Water+soil index (1D)
    water_soil_idx = Nz(3)+Nz(8);
    
    % Radial soil indices in 2D field
    radial_soil_top_idx        = Nz(2)+Nz(3)+Nz(4)+Nz(9)-3; % with insulation
    radial_soil_top_no_ins_idx = Nz(2)+Nz(3)+Nz(4)-2;       % without insulation


    soil_2d_first_idx = sys_top_idx+1;                                  % first index of 2D soil field in system vector
    soil_2d_top_idx = soil_2d_first_idx + Nz(12)*Nz(13) - 1;            % last index of 2D soil field in system vector

    piston_2d_first_idx = soil_2d_top_idx+1;                      % first index of 2D piston field in system vector
    piston_2d_top_idx = piston_2d_first_idx + Nz(14)*Nz(3) - 1;   % last index of 2D piston field in system vector

    % Mapping: 1D-temperature vector → 2D-temperature field for piston 
    block_P = T_Sys(end, piston_2d_first_idx:piston_2d_top_idx);
    T_Pf    = reshape(block_P, Nz(3), Nz(14));

    % Air temperature at system top (Dirichlet BC)
    T_Sys(:,sys_top_idx) = T0init(4);
    
    % Water / insulation coupling at top (Robin BC)
    T_Sys(:,sys_top_no_ins_idx) = (SW(1,4)*T_Sys(:,sys_top_no_ins_idx-1)+SW(2,4)...
                                  *T_Sys(:,sys_top_no_ins_idx+1))/((SW(1,4)+SW(2,4)));

    % Water / soil coupling condition (Robin BC)
    T_Sys(:,water_soil_idx) = (SW(1,4)*T_Sys(:,water_soil_idx+1)+SW(2,4)*T_Sys(:,water_soil_idx-1))/((SW(1,4)+SW(2,4)));


    % Set new membrane index temperature
    T_Sys(:,bypass_indices(3)) = (T_Sys(:,bypass_indices(3)-1)+T_Sys(:,bypass_indices(3)+1))/2;

    %%%------------------------------------------%%%
    % 04: Update water/piston contact temperatures at piston top and bottom based on flow direction
    %%%------------------------------------------%%%

    % Water/piston contact temperature at the piston top (water side)
    if flow >= 0
        % thermal discharging / standstill
        TK_WKO = T_Sys(end,Nz(3)+Nz(8)+Nz(2)+Nz(5)-2);
    else
        % thermal charging
        TK_WKO = T_Sys(end,Nz(3)+Nz(8)+Nz(2)+Nz(5)-1);
    end
    % Piston/water contact temperature at piston top (piston side)
    TK_KWO_f = T_Pf(end-1, :); % radial temperature distribution at piston top

    % Mixed contact temperature at piston top
    T_Pf(end,:) = (SW(1,4)*TK_WKO + SW(2,4)*TK_KWO_f) / (SW(1,4)+SW(2,4));


    % Water/piston contact temperature at piston bottom (water side)
    if flow <= 0
        % thermal charging / standstill
        TK_WKU = T_Sys(end,Nz(3)+Nz(8)+Nz(2)-1);
    else 
        % thermal discharging
        TK_WKU = T_Sys(end,Nz(3)+Nz(8)+Nz(2)-2);
    end
    % Piston/water contact temperature at piston bottom (piston side)
    TK_KWU_f = T_Pf(2, :); % radial temperature distribution at piston bottom

    % Mixed contact temperature at piston bottom
    T_Pf(1,:)   = (SW(1,4)*TK_WKU + SW(2,4)*TK_KWU_f) / (SW(1,4)+SW(2,4));

    % Update 1D system vector with new piston top and bottom temperatures
    T_Sys(end,piston_2d_first_idx:piston_2d_top_idx) = T_Pf(:).';


    %%%------------------------------------------%%%
    % 05: Mapping 1D system state <-> 2D radial soil field
    %%%------------------------------------------%%%

    % Time loop over all stored time steps for radial soil coupling
    for tt = 1:1+Nt2

        %----------------------------------------------------------
        % 5.1: Build 2D radial soil temperature field from 1D state
        %----------------------------------------------------------

        block = T_Sys(tt, sys_top_idx+1 : sys_top_idx + Nz(12)*Nz(13));
        T_REf = reshape(block, Nz(13), Nz(12));

        %----------------------------------------------------------
        % 5.2: Re-apply boundary conditions in the 2D radial soil
        %----------------------------------------------------------

        % Vertical boundaries
        T_REf(1,:) = T0init(5); % soil bottom (vertical)
        T_REf(end,:) = T0init(4); % soil bottom (vertical) 

        % Outer radial boundary
        T_REf(:,end) = T0init(5); % outer soil temperature

        % Inner radial boundary in contact with water / piston:
        % lower pressure zone
        for e = 1:Nz(2) 
            T_REf(e,1) = (SW(1,4)*T_Sys(tt,water_soil_idx+e-1)...
                         +SW(2,4)*T_REf(e,2))/(SW(1,4)+SW(2,4));
        end

        % upper pressure zone
        for e = Nz(2)+Nz(3)-1:radial_soil_top_no_ins_idx
            T_REf(e,1) = (SW(1,4)*T_Sys(tt,Nz(8)+Nz(5)+e-1)...
                         +SW(2,4)*T_REf(e,2))/(SW(1,4)+SW(2,4));
        end

        %----------------------------------------------------------
        % 5.3: Couple system insulation (1D) with radial insulation (2D)
        %----------------------------------------------------------

        % Temperature in system insulation (1D) at system top
        DTd(:,1) = T_Sys(tt, sys_top_no_ins_idx:sys_top_idx);

        % Temperature in radial insulation (2D) at second radial node
        DTd(:,2) = T_REf(radial_soil_top_no_ins_idx:radial_soil_top_idx,2); 
        
        % Average to obtain contact temperature between 1D and 2D insulation
        DTd(:,3) = (DTd(:,1)+DTd(:,2))/2;

        % Update contact column in 2D field (first radial node)
        T_REf(radial_soil_top_no_ins_idx:radial_soil_top_idx,1) = DTd(:,3);

        % Bottom contact between radial insulation and soil
        T_REf(end-Nz(9)+1,:) = (SW(2,4)*T_REf(end-Nz(9),:)+SW(3,4)*T_REf(end-Nz(9)+2,:))...
                               /(SW(2,4)+SW(3,4));

        %----------------------------------------------------------
        % 5.4: Map updated 2D temperature field back into system vector
        %----------------------------------------------------------
        
        T_Sys(tt, sys_top_idx+1 : sys_top_idx + Nz(12)*Nz(13)) = T_REf(:).';
    end

end