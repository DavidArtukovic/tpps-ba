
function T_Sys = apply_free_convection(T_Sys, dt_mix, flow, Nz, dz, SW, A)
% Applies discrete free convection mixing operator
% Sequential in time, conservative by construction
    if flow==0
        % No flow -> no free convection mixing
        return;
    end
    
    % Indices 
    replacement_top_idx    = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;
    replacement_bottom_idx = Nz(3)+Nz(8)+Nz(2)-1;
    replacement_mid_idx    = round((replacement_top_idx + replacement_bottom_idx)/2, 0);

    z_inlet_upper_idx = replacement_mid_idx; % approx. middle of upper water height
    z_inlet_lower_idx = Nz(3)+Nz(8) + round(Nz(2)*2/3, 0);  % approx. 2/3 of lower water height

    % lower water bounds
    z_w_lower_start_idx = Nz(3)+Nz(8)+1;
    z_w_lower_end_idx   = Nz(3)+Nz(8)+Nz(2)-1;

    % upper water bounds
    z_w_upper_start_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;              % first upper water node
    z_w_upper_end_idx   = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3;        % last water node before insulation

    % Parameters
    rho_w = SW(1,2);      % [kg/m^3]
    c_w   = SW(1,3);      % [J/(kg K)]
    A_hws = A(1);         % [m^2]

    v_flow = abs(flow);   % [m/s]
    Vdot   = v_flow * A_hws;
    mdot   = rho_w * Vdot;

    for tt = 2:size(T_Sys,1)
        
        % Use the already-updated previous state as baseline (sequential operator)
        Tcur = T_Sys(tt,:).';

        T_w_in_upper = Tcur(z_inlet_upper_idx);
        T_w_in_lower = Tcur(z_inlet_lower_idx);

        if flow < 0
            %--------------------------------------------------------------
            % Thermal Charging (warm water enters at top of upper volume): 
            % free convection in LOWER water volume since surplus volume needs
            % to be transfered to lower water volume.
            % Use upper inlet temperature as "inflow" temperature.
            %--------------------------------------------------------------
            Tin = T_w_in_upper;
            z_mix_idx = z_inlet_lower_idx; % set end of mix-zone to lower inlet

            if Tin > T_w_in_lower % water entering at lower inlet is WARMER than current lower inlet temp
                % upward direction in lower volume (increasing index)
                while (Tin > Tcur(z_mix_idx)) && (z_mix_idx < z_w_lower_end_idx)
                    z_mix_idx = z_mix_idx + 1;
                end
                ids_mix = z_inlet_lower_idx:z_mix_idx; % nodes where mixing occurs (upward)

            else % water entering at lower inlet is COLDER than current lower inlet temp
                % downward direction in lower volume (decreasing index)
                while (Tin < Tcur(z_mix_idx)) && (z_mix_idx > z_w_lower_start_idx)
                    z_mix_idx = z_mix_idx - 1;
                end
                ids_mix = z_inlet_lower_idx:-1:z_mix_idx; % nodes where mixing occurs (downward)
            end
            disp('I am here!')

        elseif flow > 0
            %----------------------------------------------------------
            % Thermal Discharging: (warm water exit's at top of upper volume): 
            % free convection in UPPER water volume since shortage volume needs to be replaced
            % Use lower inlet temperature as "inflow" temperature
            %----------------------------------------------------------
            Tin = T_w_in_lower; % inflow temperature at upper inlet is the temperature at lower inlet
            z_mix_idx = z_inlet_upper_idx; % set end of mix-zone to upper inlet

            if Tin > T_w_in_upper % water entering at upper inlet is WARMER than current upper inlet temp
                % upward direction in upper volume (increasing index)
                while (Tin > Tcur(z_mix_idx)) && (z_mix_idx < z_w_upper_end_idx)
                    z_mix_idx = z_mix_idx + 1;
                end
                ids_mix = z_inlet_upper_idx:z_mix_idx;

            else
                % downward direction in upper volume (decreasing index)
                while (Tin < Tcur(z_mix_idx)) && (z_mix_idx > z_w_upper_start_idx)
                    z_mix_idx = z_mix_idx - 1;
                end
                % NOTE: fixed vs your old code: use -1 step to avoid empty ids
                ids_mix = (z_inlet_upper_idx-1):-1:z_mix_idx;
            end
        end

        Tmix = Tcur(ids_mix);                               % current temperatures in the mixing zone
        I_mix = abs(sum((Tin - Tmix) * dz));                % [K*m] why abs()? - check out the integral
                                                            % in model_overview.md or 01_code_info.md

        Qdot_mix = mdot * c_w * abs(T_w_in_upper - T_w_in_lower);   % [W] % the abs term here is theoretically NOT UNDERSTOOD!
                                                                    % but seems to be necessary to get phyisically meaningful results
        

        % Bracket term in derivation (equation 15 model_overview.md)
        bracket = I_mix - (Qdot_mix * dt_mix) / (rho_w * c_w * A_hws);   % [K*m]
        z_dist = (ids_mix - z_mix_idx) * dz;   % [m]
        % Geometric factor
        L = (numel(ids_mix)-1)*dz;  % [m]
        geom = 2 * z_dist / (L^2);    % [1/m]

        %--------------------------------------------------------------
        % Apply the CLOSED-FORM update
        % T_new = Tin + geom_factor * ( I_mix - Qdot*dt_mix/(rho*c*A) )
        %--------------------------------------------------------------

        % Updated temperature profile
        Tcur(ids_mix) = Tin + geom .* bracket;
        T_Sys(tt,:) = Tcur.';
        disp('---------------------------------------------------');
        disp(['bracket            = ' num2str(bracket)])
        disp(['Tin        = ' num2str(Tin)])
        disp(['Qdot       = ' num2str(Qdot_mix)])
        disp(['mdot       = ' num2str(mdot)])
        disp(['dz         = ' num2str(dz)])
        disp(['T_w_in_upper      = ' num2str(T_w_in_upper)])
        disp(['T_w_in_lower      = ' num2str(T_w_in_lower)])
        disp(['T_w_in_upper - T_w_in_lower = ' num2str(abs(T_w_in_upper - T_w_in_lower))])
        disp(['ids_mix first = ' num2str(ids_mix(1))])
        disp(['ids_mix last  = '  num2str(ids_mix(end))])
        % disp(['abs(T_w_in_lower - T_w_in_upper) = ' num2str(abs(T_w_in_lower - T_w_in_upper))])
        disp(['Mixing zone nodes: ' num2str(numel(ids_mix))]);
        disp(['T_mix_end_before mixing        = ' num2str(Tmix(end))])
        disp(['T_mix_min        = ' num2str(Tcur(ids_mix(1)))])
        disp(['T_mix_max        = ' num2str(Tcur(ids_mix(end)))])
        disp(['z_dist_min       = ' num2str(z_dist(1))]);
        disp(['z_dist_max       = ' num2str(z_dist(end))]);
        disp('---------------------------------------------------');
        % disp('Tcur(ids_mix) =');
        % disp(Tcur(ids_mix));
        % disp(['Tnew_max        = ' num2str(max(Tnew_mix))])
        % disp(['Tnew_min        = ' num2str(min(Tnew_mix))])

        % disp(['T min mix zone: '  num2str(min(Tcur(z_w_lower_start_idx:z_w_lower_end_idx)))])
        % disp([ 'T max mix zone: ' num2str(max(Tcur(z_w_upper_start_idx:z_w_upper_end_idx)))]);
    end
end