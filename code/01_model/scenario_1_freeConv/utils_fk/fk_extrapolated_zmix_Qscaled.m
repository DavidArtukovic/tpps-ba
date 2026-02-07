
function Tcur = fk_extrapolated_zmix_Qscaled(Tcur, ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower, dt_mix, params, fklog)
% ---------------------------------------------------------------
% SUMMARY:
%   Applies FK update with extrapolated mixing height and energy scaling.
%
% DESCRIPTION:
%   Uses a virtual mixing boundary outside the physical domain and
%   rescales the linear temperature reconstruction to conserve energy
%   over the actual numerical mixing zone.
%
% INPUT / OUTPUT:
%   Tcur    - system temperature vector
%
% ADDITIONAL INPUT:
%   ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower,
%   dt_mix, params, fklog
% ---------------------------------------------------------------
    
    %--------------------------------------------------------------
    % 1) Current temperatures in numerical mixing zone
    %--------------------------------------------------------------
    Tmix = Tcur(ids_mix);
    %--------------------------------------------------------------
    % 2) Numerical integral over ACTUAL mixing zone
    %--------------------------------------------------------------
    dT_node = Tin - Tmix;  % nodal difference
    % Use trapezoid sum for integral
    I_num = abs( sum( 0.5 * (dT_node(1:end-1) + dT_node(2:end)) ) * params.dz );  % [K*m]

    %--------------------------------------------------------------
    % 3) Determine extrapolated z_mix* via GLOBAL linear regression
    %    over the current numerical mixing zone
    %--------------------------------------------------------------

    % Coordinates in meters relative to mixing boundary
    mix_height = (numel(ids_mix)-1)*params.dz;
    dTdz = (Tmix(end) - Tmix(1)) / mix_height;    % geometric mean absolute slope
   
    if abs(dTdz) < 1e-6
        return; % later: fully_mixed
        fklog('Warning: near-zero gradient in extrapolated_zmix case');
    end

    % Distance needed to reach Tin by extrapolation (positive in upward, negative in downward direction)
    dz_star = (Tin - Tmix(end)) / abs(dTdz);   % [m]

    %--------------------------------------------------------------
    % 4) Geometric integral over extrapolated mixing zone
    %--------------------------------------------------------------
    z_in  = 0;
    z_mix = (numel(ids_mix)-1) * params.dz;
    z_mix_star = z_mix + abs(dz_star);

    G = 0.5*(z_mix^2 - z_in^2) - z_mix_star*(z_mix - z_in);

    %--------------------------------------------------------------
    % 5) Calculate Inlet Energy rate Qdot_in
    %--------------------------------------------------------------
    Qdot_in = params.mdot * params.c_w * abs(T_w_in_upper - T_w_in_lower);

    E = (Qdot_in * dt_mix) / (params.rho_w * params.c_w * params.A_hws);

    %--------------------------------------------------------------
    % 6) Caculate new Slope a
    %--------------------------------------------------------------
    a = (E - I_num) / G;

    %--------------------------------------------------------------
    % 7) Apply closed-form update
    %--------------------------------------------------------------
    z = (ids_mix - z_mix_idx) * params.dz;   % physical z
    Tcur(ids_mix) = Tin + a * (z - dz_star);

    fklog('---------------------------------------------------');
    fklog('case: extrapolated_zmix');         
    fklog(['a                 = ' num2str(a)]);
    fklog(['I_num             = ' num2str(I_num)]);
    fklog(['G                 = ' num2str(G)]);
    fklog(['E                 = ' num2str(E)]);
    fklog(['dz_star           = ' num2str(dz_star)]);
    fklog(['mix_height       = ' num2str(mix_height)]);
    fklog(['T_in         = ' num2str(Tin)]);
    fklog(['T_mix_end_before_mixing        = ' num2str(Tmix(end))]);
    fklog(['T_mix_begin_before_mixing        = ' num2str(Tmix(1))]);
    fklog(['T_mix_end_after_mixing        = ' num2str(Tcur(ids_mix(end)))]);
    fklog(['T_mix_begin_after_mixing        = ' num2str(Tcur(ids_mix(1)))]);
    fklog('---------------------------------------------------');

end
