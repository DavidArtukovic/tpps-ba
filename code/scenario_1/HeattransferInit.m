%% HEATTRANSFERINIT.m
% ---------------------------------------------------------------
%  SUMMARY:
%   Computes the purely conductive initialization of the TPPS model and
%   couples the 1D system model to the 2D radial soil model.
%
% DESCRIPTION:
%   - Calls HEATSOLID.m to integrate the 1D conduction problem in piston,
%     water, air, insulation and vertical soil.
%   - Enforces time-independent Dirichlet boundary conditions for piston,
%     water and air temperatures in the system vector.
%   - Maps the 1D temperature vector of the radial soil into a 2D
%     (z,r)-temperature field.
%   - Re-applies boundary conditions in the 2D radial soil domain.
%   - Couples the 1D insulation temperatures to the 2D radial insulation
%     by averaging both sides at the contact interface.
%   - Maps the updated 2D temperature field back into the 1D system
%     temperature vector at every stored time step.
%
% INPUT:
%   t2      - time vector for the purely conductive initialization
%   IC_Sys  - initial temperature vector of the TPPS system
%   Nz      - vector with number of grid points in each subdomain
%   dz      - vertical grid spacing for the 1D domains
%   flow    - flow flag (not used in the purely conductive phase, but
%             passed to HEATSOLID for interface compatibility)
%   T0init  - vector of prescribed initial / boundary temperatures
%   SW      - parameter matrix with material and model switches
%   A       - cross-sectional areas and geometric parameters
%   z_RE    - geometric parameters of the radial soil grid
%   T_REf   - 2D temperature field in the radial soil (z,r)
%   Nt2     - number of stored time steps during initialization
%
% OUTPUT:
%   T_Sys   - time history of the 1D system temperature vector
%             Each row in the solution array T_Sys corresponds to 
%             a value returned in column vector t: T_Sys  ∈  ℝ^(Nt × Nstates)
%   T_REf   - updated 2D radial soil temperature field at final time step
%
% ---------------------------------------------------------------

function [T_Sys,T_REf] = HeattransferInit(t2,IC_Sys,Nz,dz,flow,T0init,SW,A,z_RE,T_REf,Nt2)
    
% Preallocation for contact temperatures between system insulation (1D)
% and radial insulation (2D)

DTd = zeros(Nz(9),3); % temperature isolation

% Integrate purely conductive phase of the 1D TPPS model
[~,T_Sys] = ode45(@HeatSolid,t2,IC_Sys,[],Nz,dz,flow,T0init,SW,A,z_RE);


%%%------------------------------------------%%%
% 01 Enforce Dirichlet boundary conditions in 1D system
%%%------------------------------------------%%%

% Piston temperature at top and bottom nodes (adapt for all time steps)
T_Sys(:,Nz(3)) = T0init(6); % piston top
T_Sys(:,1) = T0init(6);     % piston bottom
 
% Air temperature at top of system (above water and insulation)
sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4;
T_Sys(:,sys_top_idx) = T0init(4);

% Water / insulation temperature at system top (without insulation)
sys_top_no_ins_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3;
T_Sys(:,sys_top_no_ins_idx) = T0init(6);
T_Sys(:,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3) = T0init(6);

% Water temperature at interface to vertical soil
water_soil_idx = Nz(3)+Nz(8);
T_Sys(:,water_soil_idx) = T0init(6);

%%%------------------------------------------%%%
% 02 Mapping 1D system state <-> 2D radial soil field
%%%------------------------------------------%%%

radial_soil_top_idx        = Nz(2)+Nz(3)+Nz(4)+Nz(9)-3; % with insulation
radial_soil_top_no_ins_idx = Nz(2)+Nz(3)+Nz(4)-2;       % without insulation

% Time loop over all stored initialization time steps
% Radial Soil
for tt = 1:1+Nt2
    %--------------------------------------------------------------
    % 2.1: Build 2D radial soil temperature field from 1D system state
    %--------------------------------------------------------------
    for j = 1:Nz(13) % vertical direction (z)
        for i = 1:Nz(12) % radial direction (r)
            T_REf(j,i) = T_Sys(tt,sys_top_idx+(i-1)*Nz(13)+j);
        end
    end
    %--------------------------------------------------------------
    % 2.2: Re-apply boundary conditions in the 2D radial soil domain
    %--------------------------------------------------------------
    T_REf(1,:) = T0init(5); % soil temperature at bottom (vertical) (outer soil temp)
    T_REf(end,:) = T0init(4); % air / soil interface at the top (vertical) (air temp)
    T_REf(:,end) = T0init(5); % soil temperature at outer radial boundary (outer soil temp)
    T_REf(1:end-Nz(9),1) = T0init(6);   % Left radial boundary in contact with water/piston-soil

    %--------------------------------------------------------------
    % 2.3: Couple system insulation (1D) with radial insulation (2D)
    %--------------------------------------------------------------

    % Temperature in system insulation (1D) at system top
    DTd(:,1) = T_Sys(sys_top_no_ins_idx:sys_top_idx); 

    % Temperature in radial insulation (2D) at second radial node
    DTd(:,2) = T_REf(radial_soil_top_no_ins_idx:radial_soil_top_idx,2); 

    % Average to obtain contact temperature between 1D and 2D insulation
    DTd(:,3) = (DTd(:,1)+DTd(:,2))/2;

    % Update contact column in 2D field (first radial node)
    T_REf(radial_soil_top_no_ins_idx:radial_soil_top_idx,1) = DTd(:,3);

    % Bottom contact between radial insulation and soil
    T_REf(end-Nz(9)+1,:) = ...
        (SW(2,4)*T_REf(end-Nz(9),:) + SW(3,4)*T_REf(end-Nz(9)+2,:))...
        /(SW(2,4)+SW(3,4));
    
    %--------------------------------------------------------------
    % 2.4: Map updated 2D temperature field back into system vector
    %--------------------------------------------------------------
    for j = 1:Nz(13)
        for i = 1:Nz(12)
            T_Sys(tt,sys_top_idx+(i-1)*Nz(13)+j) = T_REf(j,i); 
        end
    end
end
end