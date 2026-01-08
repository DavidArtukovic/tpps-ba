%% Free convection + optional diffusion toy example (operator splitting)
clear; clc; close all;

%%% -----------------------------
%  Toggle diffusion on/off
%%% -----------------------------
enableDiffusion = true;   % <-- switch here

%% Spatial grid
z  = linspace(0,20,4000)';          % [m]
dz = z(2) - z(1);                    % [m]

%% Initial temperature profile
T0 = 50+z;                   % arbitrary linear profile

%% Concave initial temperature profile with fixed boundary values
% c = 0.08;     % curvature strength (>0 => concave)
% 
% T0 = 50 + z + c * z .* (10 - z);
% disp(T0(z==0.0025));
% disp(T0(z==10));
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
dt_list = [0,1,2,3,4,5,6]*7200;       % [s]

%% Plot
figure; hold on; grid on;

for dt = dt_list
    T = T0;

    % ------------------------------------------------------------
    % (1) Free convection operator update (Schäfer-consistent)
    % ------------------------------------------------------------

    ids = (z >= 10 & z <= 20); % upward mixing
    Tin  = 70;      % [°C] inlet temperature
    z_in = 10;       % [m]
    z_mix = 20;     % [m]
      
    % ids = (z >= 0 & z <= 10); % downward mixing
    % Tin  = 50;      % [°C] inlet temperature
    % z_in = 10;       % [m]
    % z_mix = 0;     % [m]
    
    % Mixing length
    L = z_mix - z_in;

    Tmix = T(ids);
    zloc = z(ids);
    
    % Integral I = ∫(Tin - T(z)) dz  (energy content of profile)
    I = abs(sum((Tin - Tmix)) * dz);      % [K*m]

    % Mixing power term (choose a toy value or derive it)
    Qdot_mix = -577;                  % [W]
    A_hws = 300;                     % [m^2]
    rho = rho_w;                   % [kg/m^3]
    cp  = cp_w;                    % [J/(kg*K)]

    
    % Time increment used by the operator (your Δt)
    dt_mix = dt;                   % [s] in your loop dt is the operator interval
    
    % Bracket term in Schäfer/Gerle formula
    bracket = I - (abs(Qdot_mix) * dt_mix) / (rho * cp * A_hws);   % [K*m]
    
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

if enableDiffusion
    title('Free convection operator + diffusion (Dirichlet BC)', 'FontSize', 18);
else
    title('Free convection operator only (diffusion OFF)', 'FontSize', 18);
end

legend(["\Deltat = " + string(dt_list) + " s",'initial'], ...
       'Location','best', 'FontSize', 12);
set(gca, 'FontSize', 14);

%% ------------------------------------------------------------
% Local function: explicit diffusion with Dirichlet BC
% ------------------------------------------------------------
function T = diffuse_dirichlet_explicit(T, alpha, dz, dt_total, T_left, T_right)
%DIFFUSE_DIRICHLET_EXPLICIT Explicit FD for 1D heat equation with Dirichlet BC.
%   dT/dt = alpha * d2T/dz2
%   T(1) and T(end) are kept constant.

    if dt_total <= 0
        T(1)   = T_left;
        T(end) = T_right;
        return;
    end

    % Stability condition for explicit scheme:
    % dt <= dz^2 / (2*alpha). Use safety factor.
    dt_stable = 0.45 * dz^2 / (2*alpha);

    nSteps = max(1, ceil(dt_total / dt_stable));
    dt     = dt_total / nSteps;

    N = numel(T);
    disp(nSteps);
    for k = 1:nSteps
        Tnew = T;

        % Second derivative (central differences) on interior nodes
        Tzz = (T(3:N) - 2*T(2:N-1) + T(1:N-2)) / dz^2;

        % Explicit update
        Tnew(2:N-1) = T(2:N-1) + alpha * dt * Tzz;

        % Dirichlet BC
        Tnew(1)   = T_left;
        Tnew(end) = T_right;

        T = Tnew;
    end
end
