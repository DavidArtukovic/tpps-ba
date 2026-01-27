clc; clear; close all;

%% Grid
Nz = 200;
L  = 1.0;
dz = L/(Nz-1);
z  = linspace(0,L,Nz)';

%% Parameters
u  = 1.0;        % advection velocity
dt = 0.4*dz/u;   % CFL < 1
Nt = 10;

%% Initial condition: sharp front
T0 = zeros(Nz,1);
T0(z < 0.4) = 1.0;

%% Allocate
T_central = T0;
T_upwind  = T0;

%% Time loop
for n = 1:Nt
    Tc = T_central;
    Tu = T_upwind;

    for i = 2:Nz-1
        % Central difference
        T_central(i) = Tc(i) ...
            - dt * u * (Tc(i+1) - Tc(i-1)) / (2*dz);

        % Upwind
        T_upwind(i) = Tu(i) ...
            - dt * u * (Tu(i) - Tu(i-1)) / dz;
    end
end

%% Plot
figure; hold on; box on;
plot(z, T0,        'k--', 'LineWidth',1.5)
plot(z, T_central,'r',  'LineWidth',1.5)
plot(z, T_upwind, 'b',  'LineWidth',1.5)
legend('Initial','Central','Upwind','Location','best')
xlabel('z'); ylabel('T')
title('Central vs Upwind Advection')
