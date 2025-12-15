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
%   flow   - flow flag:
%            = 0  : no flow (thermal standstill)
%            > 0  : discharging
%            < 0  : charging
%   T0init - vector of prescribed boundary / reference temperatures
%   SW     - parameter matrix with material and model switches
%   A      - cross-sectional areas and geometric parameters
%   z_RE   - geometric parameters of the radial soil grid
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


function [T_Sys,T_REf] = HeattransferSzen(t2,IC_Sys,Nz,dz,flow,T0init,SW,A,z_RE,T_REf,Nt2)

    %--------------------------------------------------------------
    % Preallocation for contact temperatures between system insulation (1D)
    % and radial insulation (2D)
    %--------------------------------------------------------------
    DTd = zeros(Nz(9),3);  % [1D insulation, 2D insulation, contact]

    %%%------------------------------------------%%%
    % 01: Integrate transient 1D TPPS model with flow (fluid + solid)
    %%%------------------------------------------%%%

    [t2,T_Sys] = ode45(@HeatFluidSolid,t2,IC_Sys,[],Nz,dz,flow,T0init,SW,A,z_RE);
    %%%------------------------------------------%%%
    % 02: Free convection as DISCRETE mixing operator (post-ODE)
    %%%------------------------------------------%%%

    dt_mix = diff(t2);
    disp(dt_mix);
    tol_t  = 1e-9;                % time tolerance for event detection

    % Indices 
    replacement_top_idx    = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;
    replacement_bottom_idx = Nz(3)+Nz(8)+Nz(2)-1;
    replacement_mid_idx    = round((replacement_top_idx + replacement_bottom_idx)/2, 0);

    z_inlet_upper_idx = replacement_mid_idx;
    z_inlet_lower_idx = Nz(3)+Nz(8) + round(Nz(2)*2/3, 0);

    % lower water bounds
    z_w_lower_start_idx = Nz(3)+Nz(8)+1;
    z_w_lower_end_idx   = Nz(3)+Nz(8)+Nz(2)-1;

    % upper water bounds
    z_w_upper_start_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;              % first upper water node
    z_w_upper_end_idx   = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3;        % last water node before insulation interface

    % --- Parameters (same as HeatFluidSolid.m) ---
    rho_w = SW(1,2);      % [kg/m^3]
    c_w   = SW(1,3);      % [J/(kg K)]
    A_hws = A(3);         % [m^2]

    v_flow = abs(flow);   % [m/s]
    Vdot   = v_flow * A_hws;
    mdot   = rho_w * Vdot;

    for tt = 2:size(T_Sys,1)

        % Apply mixing sequentially along stored time points
        t_next_mix = t2(1);

        % Use the already-updated previous state as baseline (sequential operator)
        Tcur = T_Sys(tt,:).';

        T_w_in_upper = Tcur(z_inlet_upper_idx);
        T_w_in_lower = Tcur(z_inlet_lower_idx);

        if flow < 0
            %----------------------------------------------------------
            % Charging: free convection in LOWER water volume
            % Use upper inlet temperature as "inflow" temperature (as in old code)
            %----------------------------------------------------------
            Tin = T_w_in_upper;
            z_mix_idx = z_inlet_lower_idx;

            if T_w_in_upper > T_w_in_lower
                % upward direction in lower volume (increasing index)
                while (Tin > Tcur(z_mix_idx)) && (z_mix_idx < z_w_lower_end_idx)
                    z_mix_idx = z_mix_idx + 1;
                end
                ids_mix = z_inlet_lower_idx:z_mix_idx;
                Tmix = Tcur(ids_mix);
                I_mix = sum((Tin - Tmix) * dz);          % [K*m]
            else
                % downward direction in lower volume (decreasing index)
                while (Tin < Tcur(z_mix_idx)) && (z_mix_idx > z_w_lower_start_idx)
                    z_mix_idx = z_mix_idx - 1;
                end
                ids_mix = z_inlet_lower_idx:-1:z_mix_idx;
                Tmix = Tcur(ids_mix);
                I_mix = sum((Tin - Tmix) * dz * (-1));   % [K*m] sign fix for reversed bounds
            end

            Qdot = mdot * c_w * (T_w_in_upper - T_w_in_lower);  % [W]

        elseif flow > 0
            %----------------------------------------------------------
            % Discharging: free convection in UPPER water volume
            % Use lower inlet temperature as "inflow" temperature (as in old code)
            %----------------------------------------------------------
            Tin = T_w_in_lower;
            z_mix_idx = z_inlet_upper_idx;

            if T_w_in_lower > T_w_in_upper
                % upward direction in upper volume (increasing index)
                while (Tin > Tcur(z_mix_idx)) && (z_mix_idx < z_w_upper_end_idx)
                    z_mix_idx = z_mix_idx + 1;
                end
                ids_mix = z_inlet_upper_idx:z_mix_idx;
                Tmix = Tcur(ids_mix);
                I_mix = sum((Tin - Tmix) * dz);          % [K*m]
            else
                % downward direction in upper volume (decreasing index)
                while (Tin < Tcur(z_mix_idx)) && (z_mix_idx > z_w_upper_start_idx)
                    z_mix_idx = z_mix_idx - 1;
                end
                % NOTE: fixed vs your old code: use -1 step to avoid empty ids
                ids_mix = (z_inlet_upper_idx-1):-1:z_mix_idx;
                Tmix = Tcur(ids_mix);
                I_mix = sum((Tin - Tmix) * dz * (-1));   % [K*m]
            end

            Qdot = mdot * c_w * (T_w_in_lower - T_w_in_upper);  % [W]

        else
            % No flow -> no free convection mixing
            break;
        end

        %--------------------------------------------------------------
        % Apply the CLOSED-FORM update
        % T_new = Tin + geom_factor * ( I_mix - Qdot*dt_mix/(rho*c*A) )
        %--------------------------------------------------------------
        dz_mix = max((numel(ids_mix)-1) * dz, eps);       % [m]
        z_dist = (ids_mix - z_mix_idx).' * dz;            % [m]
        geom_factor = 2 * z_dist / (dz_mix^2);            % [1/m]

        Tnew_mix = Tin + geom_factor .* ( I_mix - (Qdot * dt_mix)/(rho_w*c_w*A_hws) );

        Tcur(ids_mix) = Tnew_mix;

        T_Sys(tt,:) = Tcur.';
    end


    %%%------------------------------------------%%%
    % 03: Define frequently used indices
    %%%------------------------------------------%%%

    % Air node at the very top of the system
    sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4;

    % Water / insulation interface at the system top (1D)
    sys_top_no_ins_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3;

    % Water+soil index (1D)
    water_soil_idx = Nz(3)+Nz(8);

    % replacement volume index (1D)
    replacement_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;
    replacement_bottom_idx = Nz(3)+Nz(8)+Nz(2)-1;

    % Radial soil indices in 2D field
    radial_soil_top_idx        = Nz(2)+Nz(3)+Nz(4)+Nz(9)-3; % with insulation
    radial_soil_top_no_ins_idx = Nz(2)+Nz(3)+Nz(4)-2;       % without insulation


    %%%------------------------------------------%%%
    % 04: Coupling conditions between 1D domains
    %%%------------------------------------------%%%

    % 3.1: Water (Volumen Saugrohr) / piston coupling at top and bottom (Dirichlet BC)
    if flow == 0
        % Thermal standstill
        % Top: piston / upper water replacement volume
        T_Sys(:,Nz(3)) = (SW(1,4)*T_Sys(:,replacement_top_idx)...
                         +SW(2,4)*T_Sys(:,Nz(3)-1))/((SW(1,4)+SW(2,4)));
        % Bottom: piston / lower water replacement volume
        T_Sys(:,1) = (SW(1,4)*T_Sys(:,replacement_bottom_idx)...
                     +SW(2,4)*T_Sys(:,2))/((SW(1,4)+SW(2,4)));
    elseif flow > 0
        % Discharging
        % Top: piston / upper water replacement volume
        T_Sys(:,Nz(3)) = (SW(1,4)*T_Sys(:,replacement_top_idx)...
                         +SW(2,4)*T_Sys(:,Nz(3)-1))/((SW(1,4)+SW(2,4)));
        % Bottom: piston / lower water replacement volume
        T_Sys(:,1) = (SW(1,4)*T_Sys(:,replacement_bottom_idx-1)+SW(2,4)*T_Sys(:,2))/((SW(1,4)+SW(2,4)));
                
    else
        % Charging
        % Top: piston / upper water replacement volume
        T_Sys(:,Nz(3)) = (SW(1,4)*T_Sys(:,replacement_top_idx+1)+SW(2,4)*T_Sys(:,Nz(3)-1))/((SW(1,4)+SW(2,4)));
          % Bottom: piston / lower water replacement volume
        T_Sys(:,1) = (SW(1,4)*T_Sys(:,replacement_bottom_idx)+SW(2,4)*T_Sys(:,2))/((SW(1,4)+SW(2,4)));
    end

    % Air temperature at system top (Dirichlet BC)
    T_Sys(:,sys_top_idx) = T0init(4);
    
    % Water / insulation coupling at top (Robin BC)
    T_Sys(:,sys_top_no_ins_idx) = (SW(1,4)*T_Sys(:,sys_top_no_ins_idx-1)+SW(2,4)...
                                  *T_Sys(:,sys_top_no_ins_idx+1))/((SW(1,4)+SW(2,4)));

    % Water / soil coupling condition (Robin BC)
    T_Sys(:,water_soil_idx) = (SW(1,4)*T_Sys(:,water_soil_idx+1)+SW(2,4)*T_Sys(:,water_soil_idx-1))/((SW(1,4)+SW(2,4)));

    %%%------------------------------------------%%%
    % 05: Mapping 1D system state <-> 2D radial soil field
    %%%------------------------------------------%%%

    % Time loop over all stored time steps for radial soil coupling
    for tt = 1:1+Nt2

        %----------------------------------------------------------
        % 5.1: Build 2D radial soil temperature field from 1D state
        %----------------------------------------------------------

        for j = 1:Nz(13) % vertical direction (z)
            for i = 1:Nz(12) % radial direction (r)
                T_REf(j,i) = T_Sys(tt,sys_top_idx+(i-1)*Nz(13)+j);
            end
        end

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
        
        for j = 1:Nz(13)
            for i = 1:Nz(12)
                T_Sys(tt,sys_top_idx +(i-1)*Nz(13)+j) = T_REf(j,i); 
            end
        end
    end

end