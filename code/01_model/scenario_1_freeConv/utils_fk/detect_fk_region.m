
function [ids_mix, Tin, z_mix_idx, T_w_in_upper, T_w_in_lower, fk_code] = detect_fk_region(Tcur, flow, params, fklog)
% ---------------------------------------------------------------
% SUMMARY:
%   Detects the free-convection mixing region and inlet conditions.
%
% DESCRIPTION:
%   Determines the active FK region, inlet temperature and numerical
%   mixing boundary based on flow direction and stratification.
%
% INPUT:
%   Tcur    - current system temperature vector
%   flow    - flow flag
%   params  - FK parameter struct
%   fklog   - logging function handle
%
% OUTPUT:
%   ids_mix        - indices of the FK mixing region
%   Tin            - inlet temperature driving free convection
%   z_mix_idx      - index of numerical mixing boundary
%   T_w_in_upper   - inlet-adjacent upper water temperature
%   T_w_in_lower   - inlet-adjacent lower water temperature
%   fk_code        - FK regime code identifier
% ---------------------------------------------------------------
    fk_code = 0;
    T_w_in_upper = Tcur(params.z_inlet_upper_idx);
    T_w_in_lower = Tcur(params.z_inlet_lower_idx);

    ids_mix = [];

    if flow < 0
        %--------------------------------------------------------------
        % Thermal Charging (warm water enters at top of upper volume): 
        % free convection in LOWER water volume since surplus volume needs
        % to be transfered to lower water volume.
        % Use upper inlet temperature as "inflow" temperature.
        %--------------------------------------------------------------
        fklog('***********************************');
        fklog('Thermal Charging detected (flow < 0)');
        Tin = T_w_in_upper;
        z_mix_idx = params.z_inlet_lower_idx; % set end of mix-zone to lower inlet
        fklog(['z_inlet_lower_idx = ' num2str(params.z_inlet_lower_idx)]);

        if Tin > T_w_in_lower
            % upward direction in lower volume (increasing index)
            while (Tin > Tcur(z_mix_idx)) && (z_mix_idx < params.z_membrane_idx)
                z_mix_idx = z_mix_idx + 1;
            end
            ids_mix = params.z_inlet_lower_idx:z_mix_idx;
            fklog('Upward mixing in lower volume detected');
            fk_code = -2;
        else
            % downward direction in lower volume (decreasing index)
            while (Tin < Tcur(z_mix_idx)) && (z_mix_idx > params.z_w_lower_start_idx)
                z_mix_idx = z_mix_idx - 1;
            end
            ids_mix = params.z_inlet_lower_idx:-1:z_mix_idx;
            fklog('Downward mixing in lower volume detected');
            fk_code = -1;
        end

    elseif flow > 0
        %----------------------------------------------------------
        % Thermal Discharging: (warm water exit's at top of upper volume): 
        % free convection in UPPER water volume since shortage volume needs to be replaced
        % Use lower inlet temperature as "inflow" temperature
        %----------------------------------------------------------
        fklog('***********************************');
        fklog('Thermal Discharging detected (flow > 0)');
        Tin = T_w_in_lower;
        z_mix_idx = params.z_inlet_upper_idx;
        fklog(['z_inlet_upper_idx = ' num2str(params.z_inlet_upper_idx)]);
        if Tin > T_w_in_upper
            % upward direction in upper volume (increasing index)
            while (Tin > Tcur(z_mix_idx)) && (z_mix_idx < params.z_w_upper_end_idx)
                z_mix_idx = z_mix_idx + 1;
            end
            ids_mix = params.z_inlet_upper_idx:z_mix_idx;
            fklog('Upward mixing in upper volume detected');
            fk_code = +2;
        else
            % downward direction in upper volume (decreasing index)
            while (Tin < Tcur(z_mix_idx)) && (z_mix_idx > params.z_membrane_idx)
                z_mix_idx = z_mix_idx - 1;
            end
            ids_mix = params.z_inlet_upper_idx:-1:z_mix_idx;
            fklog('Downward mixing in upper volume detected');
            fk_code = +1;
        end
     
    end
    fklog(['z_mix_zone_idx= '  num2str(z_mix_idx)]);
    fklog(['fk_code    = ' num2str(fk_code)]);
    fklog(['Tin        = ' num2str(Tin)]);
    fklog(['T_w_in_upper      = ' num2str(T_w_in_upper)]);
    fklog(['T_w_in_lower      = ' num2str(T_w_in_lower)]);
    fklog('***********************************');

end