function dTdt = HeatSolidPiston(t, T_vec, Nz, dz, T0init, SW, z_RP)
% -------------------------------------------------------------------------
% SUMMARY:
%   Computes time derivative of 2D piston temperature field.
%
% DESCRIPTION:
%   - Implements heat conduction in cylindrical coordinates (z, r)
%   - Includes axial diffusion and radial diffusion + (1/r dT/dr) term
%   - Applies symmetry condition at r = 0
%   - Applies Dirichlet boundary conditions at outer boundaries
%
% GOVERNING EQUATION:
%   dT/dt = α ( d²T/dz² + d²T/dr² + (1/r) dT/dr )
%
% INPUT:
%   T_vec   vectorized piston temperature field
%   Nz      grid sizes
%   dz      axial grid spacing [m]
%   T0init  boundary temperatures [°C]
%   SW      material properties (SW(2,5) = α_piston)
%   z_RP    radial geometry
%
% OUTPUT:
%   dTdt    time derivative (vectorized)
% -------------------------------------------------------------------------

    T_Pf = reshape(T_vec, Nz(3), Nz(14));
    % Initilization of temperature gradient field in piston (P)
    dTdt_Pf = zeros(Nz(3),Nz(14)); % rows -> vertical, columns -> radial
    


    %%%------------------------------------------%%%
    % 01 2D-Piston
    %%%------------------------------------------%%%
    

    % Dirichlet top and bottom
    T_Pf(1,:)   = T0init(6);
    T_Pf(end,:) = T0init(6);

    % Outer radial boundary
    T_Pf(:,end) = T0init(6);

    %--------------------------------------------------------------
    % 3.1 set gradients for 2D piston temperature field inside piston (vectorized)
    %--------------------------------------------------------------
    jIdx = 2:Nz(3)-1;      % interior z
    iIdx = 2:Nz(14)-1;     % interior r

    % vertical second derivative (central diff) for all interior nodes
    d2T_dz2 = (T_Pf(jIdx+1, iIdx) - 2*T_Pf(jIdx, iIdx) + T_Pf(jIdx-1, iIdx)) / (dz^2);

    % radial second derivative with variable spacing
    drsum = z_RP(3, iIdx);     % 1 × NrInt
    dr_f  = z_RP(4, iIdx);     % 1 × NrInt
    dr_b  = z_RP(5, iIdx);     % 1 × NrInt
    ri    = z_RP(1, iIdx);     % 1 × NrInt

    % Broadcast row vectors across all z rows via implicit expansion
    d2T_dr2 = (2 ./ drsum) .* ( ...
                (T_Pf(jIdx, iIdx+1) - T_Pf(jIdx, iIdx)) ./ dr_f  ...  % forward component
              - (T_Pf(jIdx, iIdx)   - T_Pf(jIdx, iIdx-1)) ./ dr_b );  % backward component

    % first radial derivative term (1/r * dT/dr) with variable spacing
    dT_dr = (T_Pf(jIdx, iIdx+1) - T_Pf(jIdx, iIdx-1)) ./ drsum;
    one_over_r_dTdr = (1 ./ ri) .* dT_dr;

    % assemble togehter in the interior of the piston domain
    dTdt_Pf(jIdx, iIdx) = SW(2,5) * ( d2T_dz2 + d2T_dr2 + one_over_r_dTdr );

    %--------------------------------------------------------------
    % 3.2 Axis treatment (r = 0 symmetry): vectorized over z
    %--------------------------------------------------------------
    dr0 = z_RP(5,2);  % distance between radial node 1 and 2

    % second derivative at the axis (central difference)
    d2T_dz2_axis = (T_Pf(jIdx+1,1) - 2*T_Pf(jIdx,1) + T_Pf(jIdx-1,1)) / (dz^2);      

    % radial second derivative at the axis (using ghost node approach for Neumann BC)
    radial_axis  = 4 * (T_Pf(jIdx,2) - T_Pf(jIdx,1)) / (dr0^2);  % = 2*(2*(...)/dr^2)   

    dTdt_Pf(jIdx,1) = SW(2,5) * ( d2T_dz2_axis + radial_axis ); % 1/r * dT/dr term drops out at the axis due to symmetry

    %--------------------------------------------------------------
    % 3.3 Dirichlet boundaries lead to zero temperature gradient 
    % (dT/dt = 0) at the boundary nodes, enforce this explicitly
    %--------------------------------------------------------------
    dTdt_Pf(1,:)   = 0;
    dTdt_Pf(end,:) = 0;
    dTdt_Pf(:,end) = 0;

    % Map to 1D vector format for ode45
    dTdt = dTdt_Pf(:);

end