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
%   - Maps the updated 2D temperature for soild and piston field back into 
%     the 1D system temperature vector at every stored time step.
%
% INPUT:
%   tspan   - time vector for the purely conductive initialization
%   IC_Sys  - initial temperature vector of the TPPS system
%   Nz      - vector with number of grid points in each subdomain
%   dz      - vertical grid spacing for the 1D domains
%   T0init  - vector of prescribed initial / boundary temperatures
%   SW      - parameter matrix with material and model switches
%   z_RE    - geometric parameters of the radial soil grid
%   z_RP    - geometric parameters of the radial piston grid
%
% OUTPUT:
%   T_Sys   - time history of the 1D system temperature vector
%             Each row in the solution array T_Sys corresponds to 
%             a value returned in column vector t: T_Sys  ∈  ℝ^(Nt × Nstates)
%   T_REf   - updated 2D radial soil temperature field at final time step
%   T_RPf   - updated 2D radial piston temperature field at final time step
%
% ---------------------------------------------------------------

function [T_Sys, T_REf, T_RPf] = HeattransferInit(tspan, IC_Sys, Nz, dz, T0init, SW, z_RE, z_RP)
    
    % Preallocation for contact temperatures between system insulation (1D)
    % and radial insulation (2D)

    DTd = zeros(Nz(9),3); % temperature isolation


    % ODE options with sparse Jacobian structure
    options = odeset('RelTol',1e-3, ...
                    'AbsTol',1e-5, ...
                    'MaxStep',900,...
                    'InitialStep', 10);

    [~,T_Sys] = ode45(@HeatSolid, tspan, IC_Sys, options, Nz, dz, T0init, SW, z_RE, z_RP);


    %%%------------------------------------------%%%
    % 01 Enforce Dirichlet boundary conditions in 1D system
    %%%------------------------------------------%%%
    
    % Air temperature at top of system (above water and insulation)
    sys_top_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)+Nz(9)-4;
    T_Sys(:,sys_top_idx) = T0init(4);

    % Water / insulation temperature at system top (without insulation)
    sys_top_no_ins_idx = Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3;
    T_Sys(:,sys_top_no_ins_idx) = T0init(6);

    % Water temperature at interface to vertical soil
    water_soil_idx = Nz(3)+Nz(8);
    T_Sys(:,water_soil_idx) = T0init(6);

    %%%------------------------------------------%%%
    % 02 Mapping 1D system state <-> 2D radial soil field
    %%%------------------------------------------%%%

    T_Sys900 = T_Sys(end,:); % Extract final time step of temperature vector for mapping

    radial_soil_top_idx        = Nz(2)+Nz(3)+Nz(4)+Nz(9)-3;  % with insulation
    radial_soil_top_no_ins_idx = Nz(2)+Nz(3)+Nz(4)-2;        % without insulation
    soil_2d_first_idx          = Nz(10)+1;                   % first index of 2D soil field in system vector
    piston_2d_first_idx        = Nz(10)+Nz(12)*Nz(13)+1;     % first index of 2D piston field in system vector

    %--------------------------------------------------------------
    % 2.1: Build 2D radial soil temperature field from 1D system state
    %--------------------------------------------------------------
    block = T_Sys900(soil_2d_first_idx + (0:Nz(12)*Nz(13)-1));
    T_REf = reshape(block, Nz(13), Nz(12));

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
    DTd(:,1) = T_Sys900(sys_top_no_ins_idx:sys_top_idx); 

    % Temperature in radial insulation (2D) at second radial node
    DTd(:,2) = T_REf(radial_soil_top_no_ins_idx:radial_soil_top_idx,2); 

    % Average to obtain contact temperature between 1D and 2D insulation
    DTd(:,3) = (DTd(:,1)+DTd(:,2))/2;

    % Update contact column in 2D field (first radial node)
    T_REf(radial_soil_top_no_ins_idx:radial_soil_top_idx, 1) = DTd(:,3);

    % Bottom contact between radial insulation and soil
    T_REf(end-Nz(9)+1,:) = ...
        (SW(2,4)*T_REf(end-Nz(9),:) + SW(3,4)*T_REf(end-Nz(9)+2,:))...
        /(SW(2,4)+SW(3,4));
    
    %--------------------------------------------------------------
    % 2.4: Map updated 2D-soil temperature field back into system vector
    %--------------------------------------------------------------
    T_Sys900(soil_2d_first_idx + (0:Nz(12)*Nz(13)-1)) = T_REf(:);

    %--------------------------------------------------------------
    % 2.5: Build 2D radial piston temperature field from 1D system state
    %--------------------------------------------------------------
    block = T_Sys900(piston_2d_first_idx + (0:Nz(14)*Nz(3)-1));
    T_RPf = reshape(block, Nz(3), Nz(14));

    %--------------------------------------------------------------
    % 2.6: Re-apply boundary conditions in the 2D piston domain
    %--------------------------------------------------------------
    T_RPf(1,:) = T0init(6);   % piston temperature at bottom
    T_RPf(end,:) = T0init(6); % piston temperature at top
    T_RPf(:,end) = T0init(6); % piston temperature at outer radial boundary

    %--------------------------------------------------------------
    % 2.7: Map updated 2D-piston temperature field back into system vector
    %--------------------------------------------------------------
    T_Sys900(piston_2d_first_idx + (0:Nz(14)*Nz(3)-1)) = T_RPf(:);

    T_Sys(end,:) = T_Sys900; % Update final time step of system temperature vector with updated 2D fields
end
