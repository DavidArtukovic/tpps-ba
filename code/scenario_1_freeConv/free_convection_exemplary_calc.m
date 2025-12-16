%% Free convection toy example (operator update)
clear; clc; close all;

%% Spatial grid
z = linspace(-10,10,400)';   % [m]

%% Initial temperature profile

T0 = 50+z;                   % arbitrary linear profile

T0 = 50 + z;                         % arbitrary linear profile
%% Concave initial temperature profile with fixed boundary values
c = 0.08;     % curvature strength (>0 => concave)

T0 = 50 + z + c * z .* (10 - z);
disp(T0(z==0.0025));
disp(T0(z==10));
%% Dirichlet boundary conditions (constant temperatures)
T_left  = T0(1);
T_right = T0(end);

%% -----------------------------
%  Water material parameters (approx. at 25°C, 1 bar)
%% -----------------------------
rho_w = 997.6;        % [kg/m^3]
cp_w  = 4.177;         % [kJ/(kg*K)]
lambda_w = 0.6081;    % [W/(m*K)]
alpha = lambda_w / (rho_w * cp_w);   % [m^2/s]


%% Time steps to visualize
dt_list = [0,1,2,3,4,5,6]*7200;   % [s]

%% Plot
figure; hold on; grid on;

% plot(T0, z, 'k--', 'LineWidth', 1.5);

for dt = dt_list
    T = T0;


    % apply update only in mixing region z ∈ [0,10]

    % ------------------------------------------------------------
    % (1) Free convection operator update (Schäfer-consistent)
    % ------------------------------------------------------------

    ids = (z >= 0 & z <= 10);
    
    Tin  = 60;      % [°C] inlet temperature
    z_in = 0;       % [m]
    z_mix = 10;     % [m]
    L = z_mix - z_in;
    
    % Current profile in mixing zone
    Tmix = T(ids);
    zloc = z(ids);


    T(ids) = 60 + ...
             (z(ids) - 10)* (1 - 9.2e-6 * dt);

    
    % Mixing length
    L = z_mix - z_in;
    
    % Integral I = ∫(Tin - T(z)) dz  (energy content of profile)
    I = sum((Tin - Tmix)) * dz;      % [K*m]

    % Mixing power term (choose a toy value or derive it)
    Qdot_mix = 577;                  % [W] set nonzero if you want
    A_hws = 300;                     % [m^2] toy cross-section (set your value)
    rho = rho_w;                   % [kg/m^3]
    cp  = cp_w;                    % [J/(kg*K)]

    
    % Time increment used by the operator (your Δt)
    dt_mix = dt;                   % [s] in your loop dt is the operator interval
    
    % Bracket term in Schäfer/Gerle formula
    bracket = I - (Qdot_mix * dt_mix) / (rho * cp * A_hws);   % [K*m]
    
    % Geometric factor
    geom = 2 * (zloc - z_mix) / (L^2);    % [1/m]
    
    % Updated temperature profile
    T(ids) = Tin + geom .* bracket;

    % ------------------------------------------------------------
    % (2) Optional diffusion step: dT/dt = alpha * d²T/dz²
    % ------------------------------------------------------------
    if enableDiffusion
        T = diffuse_dirichlet_explicit(T, alpha, dz, dt, T_left, T_right);
    end


    plot(T, z, 'LineWidth', 1.5);
end
plot(T0, z, 'k--', 'LineWidth', 1.5);

xlabel('T', 'FontSize', 18);
ylabel('z', 'FontSize', 18);
title('Closed-form free convection update', 'FontSize', 20);
legend(["Δt = " + string(dt_list) + " s",'initial'], 'Location','best');
title('Closed-form free convection update');
legend(["\Deltat = " + string(dt_list) + " s",'initial'], ...
       'Location','best', 'FontSize', 18);

set(gca, 'FontSize', 14);   % tick labels