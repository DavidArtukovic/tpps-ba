
function params = init_fk_parameters(Nz, dz, SW, A, bypass_indices, flow)
% ---------------------------------------------------------------
% SUMMARY:
%   Initializes geometric indices and physical parameters for FK modeling.
%
% INPUT:
%   Nz    - number of grid points per subdomain
%   dz    - vertical grid spacing [m]
%   SW    - material parameters
%   A     - geometric parameters

%   flow  - flow flag
%
% OUTPUT:
%   params - struct containing FK-related indices and physical constants
% ---------------------------------------------------------------

    % Indices 

    params.z_inlet_lower_idx = bypass_indices(1);  
    params.z_inlet_upper_idx = bypass_indices(2);
    params.z_membrane_idx    = bypass_indices(3);      


    % lower water bounds
    params.z_w_lower_start_idx = Nz(3)+Nz(8)+1;
    params.z_w_lower_end_idx   = Nz(3)+Nz(8)+Nz(2)-1;

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
