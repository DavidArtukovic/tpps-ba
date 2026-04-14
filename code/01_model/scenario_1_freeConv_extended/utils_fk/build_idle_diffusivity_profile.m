
function alpha_w_eff = build_idle_diffusivity_profile(T, Nz, dz, SW, C_fk, k_fk, grad_tol, alpha_cap)
% ---------------------------------------------------------------
% SUMMARY:
%   Builds a nodal effective axial diffusivity profile for idle water
%   regions based on inverse thermocline detection.
%
% INPUT:
%   T         - current state vector
%   Nz        - number of grid points per subdomain
%   dz        - axial grid spacing [m]
%   SW        - material parameter matrix
%   C_fk      - empirical coefficient in lambda_enh = C * |dT/dz|^k
%   k_fk      - exponent for enhanced conductivity [-]
%   grad_tol  - activation threshold for inverse gradient [K/m]
%   alpha_cap - upper bound for effective diffusivity [m^2/s]
%
% OUTPUT:
%   alpha_w_eff - nodal diffusivity vector on full state vector [m^2/s]
% ---------------------------------------------------------------

    alpha_w_eff = SW(1,5) * ones(numel(T),1);
    lambda_w = SW(1,1);
    rho_w    = SW(1,2);
    cp_w     = SW(1,3);

    [mask_full, grad_full] = detect_inverse_thermocline(T, Nz, dz, grad_tol);

    active_idx = find(mask_full);
    for m = 1:numel(active_idx)
        i = active_idx(m);

        lambda_enh = C_fk * abs(grad_full(i))^k_fk;
        lambda_eff = max(lambda_w, lambda_enh);
        alpha_eff  = lambda_eff / (rho_w * cp_w);
        alpha_w_eff(i) = min(alpha_eff, alpha_cap);
    end
end
