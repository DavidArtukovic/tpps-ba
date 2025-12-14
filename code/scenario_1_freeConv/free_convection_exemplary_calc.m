%% Free convection toy example (operator update)
clear; clc; close all;

%% Spatial grid
z = linspace(-10,10,400)';   % [m]

%% Initial temperature profile
T0 = 50+z;                   % arbitrary linear profile

%% Time steps to visualize
dt_list = [0, 900, 1800, 2700, 3600];   % [s]

%% Plot
figure; hold on; grid on;

plot(T0, z, 'k--', 'LineWidth', 1.5);

for dt = dt_list
    T = T0;

    % apply update only in mixing region z ∈ [0,10]
    ids = (z >= 0 & z <= 10);

    T(ids) = 60 + ...
             (z(ids) - 10)* (1 - 9.2e-6 * dt);

    plot(T, z, 'LineWidth', 1.5);
end

xlabel('T');
ylabel('z');
legend(['initial', "Δt = " + string(dt_list) + " s"], 'Location','best');
title('Operator-type free convection update');
