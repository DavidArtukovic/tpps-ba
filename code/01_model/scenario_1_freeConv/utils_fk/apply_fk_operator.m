
function Tcur = apply_fk_operator(Tcur, ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower, fk_case, dt_mix, params,fklog)
% ---------------------------------------------------------------
% SUMMARY:
%   Applies the FK temperature update for the classified regime.
%
% DESCRIPTION:
%   Dispatches to the appropriate closed-form FK operator depending on
%   the detected numerical free-convection case.
%
% INPUT:
%   Tcur            - current system temperature vector
%   ids_mix         - mixing region indices
%   z_mix_idx       - numerical mixing boundary index
%   Tin             - inlet temperature
%   T_w_in_upper    - upper inlet-adjacent temperature
%   T_w_in_lower    - lower inlet-adjacent temperature
%   fk_case         - classified FK regime
%   dt_mix          - FK time step size [s]
%   params          - FK parameter struct
%   fklog           - logging function handle
%
% OUTPUT:
%   Tcur            - updated temperature vector after FK operator
% ---------------------------------------------------------------
    switch fk_case

        case 'internal_zmix'
            %--------------------------------------------------------------
            % Internal mixing zone:
            % The inlet temperature lies within the temperature range of
            % the mixing zone.
            %
            % Physical interpretation:
            % Free convection develops only up to the height where the
            % local storage temperature equals the inlet temperature.
            % Above this point, the stratification remains unaffected.
            %
            % Modeling assumption:
            % The thermal energy introduced at the inlet is redistributed
            % within the physically bounded mixing zone using a linear
            % temperature reconstruction anchored at z_mix.
            % Energy conservation is enforced via integral scaling.
            %--------------------------------------------------------------
            Tmix = Tcur(ids_mix);                               % current temperatures in the mixing zone

            dT_node = Tin - Tmix;  % nodal difference
            % Use trapezoid sum for integral
            I_mix = abs( sum( 0.5 * (dT_node(1:end-1) + dT_node(2:end)) ) * params.dz );  % [K*m]

            Qdot_mix = params.mdot * params.c_w * abs(T_w_in_upper - T_w_in_lower);        % [W]

            % Bracket term in derivation (equation 15 model_overview.md)
            bracket = I_mix - (Qdot_mix * dt_mix) / (params.rho_w * params.c_w * params.A_hws);

            z_dist = (ids_mix - z_mix_idx) * params.dz;         % [m]
            L = (numel(ids_mix)-1) * params.dz;                 % [m]
            geom = 2 * z_dist / (L^2);                           % [1/m]

            %--------------------------------------------------------------
            % Apply the CLOSED-FORM update
            % T_new = Tin + geom_factor * ( I_mix - Qdot*dt_mix/(rho*c*A) )
            %--------------------------------------------------------------
            Tcur(ids_mix) = Tin + geom .* bracket;
            fklog('---------------------------------------------------');
            fklog('case: internal_zmix');
            fklog(['bracket            = ' num2str(bracket)]);
            fklog(['Qdot       = ' num2str(Qdot_mix)]);
            fklog(['last mixing node = '  num2str(z_dist(end))]);
            fklog(['T_mix_end_before_mixing        = ' num2str(Tmix(end))]);
            fklog(['T_mix_begin_before_mixing        = ' num2str(Tmix(1))]);
            fklog(['T_mix_end_after_mixing        = ' num2str(Tcur(ids_mix(end)))]);
            fklog(['T_mix_begin_after_mixing        = ' num2str(Tcur(ids_mix(1)))]);
            fklog('---------------------------------------------------');

        case 'extrapolated_zmix'
            %--------------------------------------------------------------
            % Extrapolated mixing zone:
            % Inlet temperature lies OUTSIDE the temperature range of the
            % mixing zone. The theoretical mixing height lies outside the
            % physical domain and is extrapolated.
            %
            % Physical interpretation:
            % The entire FK region participates in the convective process,
            % but the linear temperature reconstruction is anchored at the
            % VIRTUAL mixing height.
            %
            % Energy conservation is ensured by scaling the reconstructed
            % temperature profile using the integrated energy content.
            %--------------------------------------------------------------
            Tcur = fk_extrapolated_zmix_Qscaled(Tcur, ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower, dt_mix, params, fklog);

        case 'fully_mixed'
            %--------------------------------------------------------------
            % Fully mixed case:
            % The temperature variation within the mixing zone is BELOW
            % a defined threshold. No meaningful stratification exists.
            %
            % Physical interpretation:
            % Free convection faces an homogenized FK region.
            % Any additional thermal energy introduced at the inlet
            % leads to a UNIFORM temperature shift of the mixing zone.
            %
            % Modeling assumption:
            % The total thermal energy associated with the inlet mass flow
            % (relative to the inlet-adjacent temperature) is distributed
            % linearly over the entire mixing zone.
            %--------------------------------------------------------------
            fklog('---------------------------------------------------');
            fklog('case: fully_mixed');
            fklog(['T_mix_end_before_mixing        = ' num2str(Tcur(ids_mix(end)))]);
            fklog(['T_mix_begin_before_mixing        = ' num2str(Tcur(ids_mix(1)))]);        

            % Reference temperature at inlet height
            T_ref = Tcur(ids_mix(1));

            % Thermal energy rate introduced by inlet
            Qdot_mix = params.mdot * params.c_w * (Tin - T_ref);

            % Total energy added over dt
            Q_add = Qdot_mix * dt_mix;

            % Volume / mass of mixing zone
            M_mix = params.rho_w * params.A_hws * params.dz * (numel(ids_mix)-1);

            % Uniform temperature shift
            dT = Q_add / (M_mix * params.c_w);

            % Apply uniform shift (positive or negative)
            Tcur(ids_mix) = Tcur(ids_mix) + dT;
            
            fklog(['dT uniform shift   = ' num2str(dT)]);
            fklog(['T_mix_end_after_mixing        = ' num2str(Tcur(ids_mix(end)))]);
            fklog(['T_mix_begin_after_mixing        = ' num2str(Tcur(ids_mix(1)))]);
            fklog('---------------------------------------------------');
    end
end