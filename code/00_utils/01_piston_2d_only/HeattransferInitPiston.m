function [T_RPf] = HeattransferInitPiston(tspan, T_RPf0, Nz, dz, T0init, SW, z_RP)
% -------------------------------------------------------------------------
% SUMMARY:
%   Time integration wrapper for 2D piston heat conduction.
%
% DESCRIPTION:
%   - Flattens 2D piston field to vector form
%   - Solves transient heat conduction using ode45
%   - Reshapes solution back to 2D field
%   - Applies Dirichlet boundary conditions after integration
%
% INPUT:
%   tspan   [1x2]   time interval [s]
%   T_RPf0  [Nz(3) x Nz(14)] initial piston temperature field [°C]
%   Nz      struct  grid sizes
%   dz      scalar  axial grid spacing [m]
%   T0init  vector  boundary temperatures [°C]
%   SW      matrix  material properties
%   z_RP    matrix  radial geometry
%
% OUTPUT:
%   T_RPf   [Nz(3) x Nz(14)] updated piston temperature field [°C]
% -------------------------------------------------------------------------   

    % ODE options with sparse Jacobian structure
    options = odeset('RelTol',1e-3, ...
                    'AbsTol',1e-5, ...
                    'MaxStep',900,...
                    'InitialStep', 10);

    T0_vec = T_RPf0(:);

    [~, T_sol] = ode45(@HeatSolidPiston, tspan, T0_vec, options, Nz, dz, T0init, SW, z_RP);

    T_RPf = reshape(T_sol(end,:), Nz(3), Nz(14));
    %--------------------------------------------------------------
    % 1: Set boundary conditions in the 2D piston domain
    %--------------------------------------------------------------
    T_RPf(1,:) = T0init(6);   % piston temperature at bottom
    T_RPf(end,:) = T0init(6); % piston temperature at top
    T_RPf(:,end) = T0init(6); % piston temperature at outer radial boundary

end
