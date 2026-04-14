
function [mask_full, grad_full] = detect_inverse_thermocline(T, Nz, dz, grad_tol)
% ---------------------------------------------------------------
% SUMMARY:
%   Detects inverse thermocline regions in the lower and upper water
%   volumes using a gradient tolerance.
%
% INPUT:
%   T        - current 1D/extended temperature state vector
%   Nz       - number of grid points per subdomain
%   dz       - axial grid spacing [m]
%   grad_tol - activation threshold for inverse gradient [K/m]
%
% OUTPUT:
%   mask_full - logical mask on full state vector, true where inverse
%               thermocline is detected in lower/upper water
%   grad_full - axial temperature gradient dT/dz on the same entries [K/m]
% ---------------------------------------------------------------

    n_states = numel(T);
    mask_full = false(n_states,1);
    grad_full = zeros(n_states,1);

    idx_lower = (Nz(3)+Nz(8)+1):(Nz(3)+Nz(8)+Nz(2)-1);
    idx_upper = (Nz(3)+Nz(8)+Nz(2)+Nz(5)-2):(Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3);

    local_fill(idx_lower);
    local_fill(idx_upper);

    function local_fill(idx_range)
        if numel(idx_range) < 3
            return;
        end

        for k = 2:numel(idx_range)-1
            i = idx_range(k);
            grad_i = (T(i+1) - T(i-1)) / (2*dz);
            grad_full(i) = grad_i;

            % Inverse thermocline: temperature increases downward
            if grad_i < -grad_tol
                mask_full(i) = true;
            end
        end

        % Edge nodes use one-sided gradients
        i = idx_range(1);
        grad_full(i) = (T(i+1) - T(i)) / dz;
        if grad_full(i) < -grad_tol
            mask_full(i) = true;
        end

        i = idx_range(end);
        grad_full(i) = (T(i) - T(i-1)) / dz;
        if grad_full(i) < -grad_tol
            mask_full(i) = true;
        end
    end
end
