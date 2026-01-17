
function T_Sys = apply_free_convection(T_Sys, dt_mix, flow, Nz, dz, SW, A, fklog)
% Applies discrete free convection mixing operator
% Sequential in time, conservative by construction

    if flow==0
        % No flow -> no free convection mixing
        return;
    end

    %--------------------------------------------------------------
    % Precompute geometric indices and physical parameters
    %--------------------------------------------------------------
    params = init_fk_parameters(Nz, dz, SW, A, flow);
        
    % Use the already-updated previous state as baseline (sequential operator)
    Tcur = T_Sys(end,:).';

    %----------------------------------------------------------
    % 1) Detect free-convection region and inflow temperature
    %----------------------------------------------------------
    [ids_mix, Tin, z_mix_idx, T_w_in_upper, T_w_in_lower] = detect_fk_region(Tcur, flow, params, fklog);

    if isempty(ids_mix)
        return  ; % no mixing region detected
    end

    %----------------------------------------------------------
    % 2) Classify numerical free-convection case
    %----------------------------------------------------------
    fk_case = classify_fk_case(Tcur(ids_mix), Tin, 1);  % last number indicates threshold for fully mixed [K]

    %----------------------------------------------------------
    % 3) Apply appropriate free-convection operator
    %----------------------------------------------------------
    Tcur_temp = apply_fk_operator(Tcur, ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower, fk_case, dt_mix, params, fklog);

    %----------------------------------------------------------
    % 4) FK boundary check + lambda hybrid (if needed)
    %----------------------------------------------------------
    Tcur_new = apply_fk_lambda_hybrid(Tcur, Tcur_temp, Tin, ids_mix, params, fk_case, dt_mix, fklog);

    T_Sys(end,:) = Tcur_new.';

end

function params = init_fk_parameters(Nz, dz, SW, A, flow)

    % Indices 
    params.replacement_top_idx    = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;
    params.replacement_bottom_idx = Nz(3)+Nz(8)+Nz(2)-1;
    params.replacement_mid_idx    = round((params.replacement_top_idx + ...
                                           params.replacement_bottom_idx)/2,0);
    % static case okay, dynamic case: mid index changes with water levels

    params.z_inlet_upper_idx = params.replacement_mid_idx+1;          % approx. middle of upper water height
    params.z_inlet_lower_idx = Nz(3)+Nz(8)+round(Nz(2)*9/10,0);        % approx. 9/10 of lower water height


    % lower water bounds
    params.z_w_lower_start_idx = Nz(3)+Nz(8)+1;
    % params.z_w_lower_end_idx   = Nz(3)+Nz(8)+Nz(2)-1;
    params.z_w_lower_end_idx = params.replacement_mid_idx-1;

    % upper water bounds
    params.z_w_upper_start_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;
    params.z_w_upper_end_idx   = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3;

    % Parameters
    params.rho_w = SW(1,2);      % [kg/m^3]
    params.c_w   = SW(1,3);      % [J/(kg K)]
    params.A_hws = A(1);         % [m^2]

    v_flow = abs(flow);          % [m/s]
    Vdot   = v_flow * params.A_hws;
    params.mdot = params.rho_w * Vdot;

    params.dz = dz;
end

function [ids_mix, Tin, z_mix_idx, T_w_in_upper, T_w_in_lower] = detect_fk_region(Tcur, flow, params, fklog)

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

        if Tin > T_w_in_lower
            % upward direction in lower volume (increasing index)
            while (Tin > Tcur(z_mix_idx)) && (z_mix_idx < params.z_w_lower_end_idx)
                z_mix_idx = z_mix_idx + 1;
            end
            ids_mix = params.z_inlet_lower_idx:z_mix_idx;
            fklog('Upward mixing in lower volume detected');

        else
            % downward direction in lower volume (decreasing index)
            while (Tin < Tcur(z_mix_idx)) && (z_mix_idx > params.z_w_lower_start_idx)
                z_mix_idx = z_mix_idx - 1;
            end
            ids_mix = params.z_inlet_lower_idx:-1:z_mix_idx;
            fklog('Downward mixing in lower volume detected');

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

        if Tin > T_w_in_upper
            % upward direction in upper volume (increasing index)
            while (Tin > Tcur(z_mix_idx)) && (z_mix_idx < params.z_w_upper_end_idx)
                z_mix_idx = z_mix_idx + 1;
            end
            ids_mix = params.z_inlet_upper_idx:z_mix_idx;
            fklog('Upward mixing in upper volume detected');
        else
            % downward direction in upper volume (decreasing index)
            while (Tin < Tcur(z_mix_idx)) && (z_mix_idx > params.z_w_upper_start_idx)
                z_mix_idx = z_mix_idx - 1;
            end
            ids_mix = (params.z_inlet_upper_idx-1):-1:z_mix_idx;
            fklog('Downward mixing in upper volume detected');
        end
        
        
    end
   
    fklog(['Tin        = ' num2str(Tin)]);
    fklog(['T_w_in_upper      = ' num2str(T_w_in_upper)]);
    fklog(['T_w_in_lower      = ' num2str(T_w_in_lower)]);
    fklog('***********************************');

end

function fk_case = classify_fk_case(Tmix, Tin, dT_fully_mixed)
    % Classifies free convection regime
    %
    % Tmix            : temperature vector in mixing region
    % Tin             : inlet temperature
    % dT_fully_mixed  : threshold for fully mixed assumption [K]

    Tmin = min(Tmix);
    Tmax = max(Tmix);

    % Fully mixed: almost uniform temperature profile
    if (Tmax - Tmin) <= dT_fully_mixed
        fk_case = 'fully_mixed';

    % Inlet temperature lies within temperature range
    elseif  Tin >= min(Tmix) && Tin <= max(Tmix)
        fk_case = 'internal_zmix';

    % Inlet temperature outside temperature range
    else
        fk_case = 'extrapolated_zmix';
    end
end

function Tcur = apply_fk_operator(Tcur, ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower, fk_case, dt_mix, params,fklog)
    % Applies free convection operator depending on classified FK case
    %
    % fk_case : 'internal_zmix' | 'extrapolated_zmix' | 'fully_mixed'
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
            I_mix = abs(sum((Tin - Tmix) * params.dz));         % [K*m]

            Qdot_mix = params.mdot * params.c_w * ...
                    abs(T_w_in_upper - T_w_in_lower);        % [W]

            % Bracket term in derivation (equation 15 model_overview.md)
            bracket = I_mix - (Qdot_mix * dt_mix) / ...
                    (params.rho_w * params.c_w * params.A_hws);

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
            Tcur = fk_extrapolated_zmix_Qscaled(Tcur, ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower, dt_mix, params,fklog);

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

function Tcur = fk_extrapolated_zmix_Qscaled(Tcur, ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower, dt_mix, params, fklog)
    % Free convection with extrapolated z_mix and Q-scaling
    % Ensures that full Q_in*dt is allocated in the actual numerical mixing zone
    
    %--------------------------------------------------------------
    % 1) Current temperatures in numerical mixing zone
    %--------------------------------------------------------------
    Tmix = Tcur(ids_mix);

    %--------------------------------------------------------------
    % 2) Numerical integral over ACTUAL mixing zone
    %--------------------------------------------------------------
    I_num = abs(sum((Tin - Tmix) * params.dz));   % [K*m]

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
    dz_star = abs((Tin - Tmix(end)) / dTdz);   % [m]

    %--------------------------------------------------------------
    % 4) Geometric integral over extrapolated mixing zone
    %--------------------------------------------------------------
    z_in  = 0;
    z_mix = (numel(ids_mix)-1) * params.dz;
    z_mix_star = z_mix + dz_star;

    G = 0.5*(z_mix^2 - z_in^2) - z_mix_star*(z_mix - z_in);

    %--------------------------------------------------------------
    % 5) Calculate Inlet Energy rate Qdot_in
    %--------------------------------------------------------------
    Qdot_in = params.mdot * params.c_w * ...
          (T_w_in_upper - T_w_in_lower);

    E = (Qdot_in * dt_mix) / ...
        (params.rho_w * params.c_w * params.A_hws);

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

function Tnew = apply_fk_lambda_hybrid(Told, Tfk, Tin, ids_mix, params, fk_case, dt_mix, fklog)
    % Applies lambda-weighted FK/FM hybrid if FK violates boundary conditions

    Tnew = Tfk;   % default: accept FK result

    % --- only for relevant linear approaches cases ---
    if ~ismember(fk_case, {'internal_zmix','extrapolated_zmix'})   
        return
    end

    % --- indices ---
    idx_zin  = ids_mix(1);
    idx_zmix = ids_mix(end);

    % --- boundary values ---
    Told_in  = Told(idx_zin);
    Told_mix = Told(idx_zmix);

    Tfk_in   = Tfk(idx_zin);
    Tfk_mix  = Tfk(idx_zmix);

    % --- check violation ---
    tol = 0.2;   % temperature tolerance [°C]

    violates = (Tfk_in  < Told_in  - tol) || ...
            (Tfk_mix < Told_mix - tol);

    if ~violates
        return
    end

    fklog('--- FK boundary violation detected: applying FK/FM hybrid ---');
    % Thermal energy rate introduced by inlet
    Qdot_mix = params.mdot * params.c_w * (Tin - Told_in);
    
    % Total energy added over dt
    Q_add = Qdot_mix * dt_mix;

    % Volume / mass of mixing zone
    M_mix = params.rho_w * params.A_hws * params.dz * (numel(ids_mix)-1);

    % Uniform temperature shift
    dT_fm = Q_add / (M_mix * params.c_w);

    Tfm = Told;
    Tfm(ids_mix) = Told(ids_mix) + dT_fm;

    % --- lambda from boundary constraints ---
    epsT = 1e-12;
    lambda_in  = (Tfm(idx_zin)  - Told_in ) / ...
                (Tfm(idx_zin)  - Tfk_in + epsT);

    lambda_mix = (Tfm(idx_zmix) - Told_mix) / ...
                (Tfm(idx_zmix) - Tfk_mix + epsT);

    lambda = min([lambda_in, lambda_mix, 1]);
    lambda = max(lambda,0);

    % --- apply lambda blend in mixing zone ---
    Tnew(ids_mix) = ...
        lambda * Tfk(ids_mix) + ...
        (1-lambda) * Tfm(ids_mix);

    % --- logging ---
    fklog(sprintf('FK limited by lambda = %.3f (case: %s)', lambda, fk_case));
    fklog(['T_in_before_hybrid       = ' num2str(Tfk_in)]);
    fklog(['T_in_after_hybrid        = ' num2str(Tnew(idx_zin))]);
    fklog(['T_mix_before_hybrid      = ' num2str(Tfk_mix)]);
    fklog(['T_mix_after_hybrid       = ' num2str(Tnew(idx_zmix))]);
    fklog(sprintf('FK hybrid applied: lambda=%.3f, T_in %.3f→%.3f, T_mix %.3f→%.3f', ...
    lambda, Tfk_in, Tnew(idx_zin), Tfk_mix, Tnew(idx_zmix)));

end




