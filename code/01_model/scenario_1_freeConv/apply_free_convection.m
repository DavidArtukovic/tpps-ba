function [T_Sys, fk_code] = apply_free_convection(T_Sys, dt_mix, flow, Nz, dz, SW, A, bypass_indices, fklog)
% ---------------------------------------------------------------
% SUMMARY:
%   Applies the discrete free-convection (FK) mixing operator to the
%   1D system temperature state in a sequential and energy-conservative way.
%
% DESCRIPTION:
%   - Detects the free-convection mixing region depending on flow direction.
%   - Classifies the numerical FK regime (internal, extrapolated, fully mixed).
%   - Applies the corresponding FK temperature update.
%   - Enforces physical boundary consistency via FK/FM lambda hybrid.
%
% INPUT:
%   T_Sys   - system temperature history [Nt x Nstates]
%   dt_mix  - FK time step size [s]
%   flow    - flow flag (0: none, >0: discharge, <0: charge)
%   Nz      - number of grid points per subdomain
%   dz      - vertical grid spacing [m]
%   SW      - material and model parameter matrix
%   A       - cross-sectional and geometric parameters
%   bypass_indices - indices of the lower/upper bypass inlet and membrane in the 1D system state vector
%   fklog   - logging function handle
%
% OUTPUT:
%   T_Sys   - updated system temperature history
%   fk_code - FK regime code applied in this step
% ---------------------------------------------------------------
    fk_code = 0;

    if flow==0
        % No flow -> no free convection mixing
        return;
    end

    %--------------------------------------------------------------
    % Precompute geometric indices and physical parameters
    %--------------------------------------------------------------
    params = init_fk_parameters(Nz, dz, SW, A, bypass_indices, flow);
        
    % Use the already-updated previous state as baseline (sequential operator)
    Tcur = T_Sys(end,:).';

    %----------------------------------------------------------
    % 1) Detect free-convection region and inflow temperature
    %----------------------------------------------------------
    [ids_mix, Tin, z_mix_idx, T_w_in_upper, T_w_in_lower, fk_code] = detect_fk_region(Tcur, flow, params, fklog);

    if numel(ids_mix) < 3
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
    % 4) FK boundary check + lambda hybrid (if not needed, returns original FK result)
    %----------------------------------------------------------
    Tcur_new = apply_fk_lambda_hybrid(Tcur, Tcur_temp, Tin, ids_mix, params, fk_case, dt_mix, fklog);

    T_Sys(end,:) = Tcur_new.';

end



