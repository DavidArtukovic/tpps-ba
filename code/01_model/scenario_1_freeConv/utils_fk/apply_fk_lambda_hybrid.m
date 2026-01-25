
function Tnew = apply_fk_lambda_hybrid(Told, Tfk, Tin, ids_mix, params, fk_case, dt_mix, fklog)
% ---------------------------------------------------------------
% SUMMARY:
%   Enforces physical boundary constraints via FK/FM lambda hybrid.
%
% DESCRIPTION:
%   Detects boundary violations in linear FK updates and blends the
%   FK result with a fully mixed solution using a lambda limiter.
%
% INPUT:
%   Told     - temperature vector before FK update
%   Tfk      - temperature vector after FK update
%   Tin      - inlet temperature
%   ids_mix  - mixing region indices
%   params   - FK parameter struct
%   fk_case  - classified FK regime
%   dt_mix   - FK time step size [s]
%   fklog    - logging function handle
%
% OUTPUT:
%   Tnew     - physically consistent hybrid temperature vector
% ---------------------------------------------------------------

    Tnew = Tfk;   % default: accept FK result

    % --- only for relevant linear approaches cases ---
    if ~ismember(fk_case, {'internal_zmix', 'extrapolated_zmix'})   
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