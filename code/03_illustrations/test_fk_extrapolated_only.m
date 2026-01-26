%% Test script: Extrapolated free convection only
% ---------------------------------------------------------------
% Tests FK extrapolation for:
%  (A) Tin < min(Tmix)
%  (B) Tin > max(Tmix)
%
% No internal mixing, no fully mixed case.
% Geometry and energy scaling only.
% ---------------------------------------------------------------
clear; clc; close all;

%% Spatial grid
z  = linspace(0,10,2000)';     % [m]
dz = z(2) - z(1);

%% Initial stratified temperature profile (monotonic)
T0 = 40 + z;               % 40°C → 48°C

%% Mixing region definition (entire column)
ids_mix = (z >= 0 & z <= 10);
z_mix_idx = find(z == 0, 1, 'first');   % numerical inlet at bottom

%% Physical parameters (toy-consistent)
params.dz    = dz;
params.rho_w = 1000;           % [kg/m^3]
params.c_w   = 4180;           % [J/(kg K)]
params.A_hws = 1.0;            % [m^2]
params.mdot  = 0.15;           % [kg/s]

dt_mix = 2*3600;               % 2 h FK step

%% Two extrapolated inlet temperatures
Tin_list = [ ...
    min(T0(ids_mix)) - 5;   % colder than entire profile
    max(T0(ids_mix)) + 5];  % hotter than entire profile

labels = ["Tin < Tmin(Tmix)", "Tin > Tmax(Tmix)"];
colors = lines(2);

figure; hold on; grid on;

%% Loop over extrapolated cases
for k = 1:2

    T = T0;
    Tin = Tin_list(k);

    % inlet-adjacent temperatures (only magnitude matters here)
    T_w_in_upper = T(ids_mix(end));
    T_w_in_lower = T(ids_mix(1));

    % ---- Apply extrapolated FK operator ----
    T = fk_extrapolated_zmix_Qscaled( ...
        T, ids_mix, z_mix_idx, Tin, ...
        T_w_in_upper, T_w_in_lower, ...
        dt_mix, params, @(s)[] );

    plot(T, z, 'LineWidth', 2, 'Color', colors(k,:));
end

%% Reference plot
plot(T0, z, 'k--', 'LineWidth', 1.5);

xlabel('Temperature [°C]');
ylabel('z [m]');
title('FK extrapolated mixing test (Tin outside temperature range)');
legend([labels, "initial"], 'Location','best');
set(gca,'FontSize',14);



function Tcur = fk_extrapolated_zmix_Qscaled(Tcur, ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower, dt_mix, params, ~)
% ---------------------------------------------------------------
% SUMMARY:
%   Applies FK update with extrapolated mixing height and energy scaling.
%
% DESCRIPTION:
%   Uses a virtual mixing boundary outside the physical domain and
%   rescales the linear temperature reconstruction to conserve energy
%   over the actual numerical mixing zone.
%
% INPUT / OUTPUT:
%   Tcur    - system temperature vector
%
% ADDITIONAL INPUT:
%   ids_mix, z_mix_idx, Tin, T_w_in_upper, T_w_in_lower,
%   dt_mix, params, fklog
% ---------------------------------------------------------------
    
    %--------------------------------------------------------------
    % 1) Current temperatures in numerical mixing zone
    %--------------------------------------------------------------
    Tmix = Tcur(ids_mix);

    %--------------------------------------------------------------
    % 2) Numerical integral over ACTUAL mixing zone
    %--------------------------------------------------------------
    I_num = abs(sum((Tin - Tmix) * params.dz));   % [K*m]

    %--------------------------------------------------------------
    % 3) Determine extrapolated z_mix* via GLOBAL linear regression
    %    over the current numerical mixing zone
    %--------------------------------------------------------------

    % Coordinates in meters relative to mixing boundary
    mix_height = (numel(ids_mix)-1)*params.dz;
    dTdz = (Tmix(end) - Tmix(1)) / mix_height;    % geometric mean absolute slope
   
    if abs(dTdz) < 1e-6
        return; % later: fully_mixed
        fklog('Warning: near-zero gradient in extrapolated_zmix case');
    end

    % Distance needed to reach Tin by extrapolation (positive in upward, negative in downward direction)
    dz_star = abs((Tin - Tmix(end)) / dTdz);   % [m]

    %--------------------------------------------------------------
    % 4) Geometric integral over extrapolated mixing zone
    %--------------------------------------------------------------
    z_in  = 0;
    z_mix = (numel(ids_mix)-1) * params.dz;
    z_mix_star = z_mix + dz_star;

    G = 0.5*(z_mix^2 - z_in^2) - z_mix_star*(z_mix - z_in);

    %--------------------------------------------------------------
    % 5) Calculate Inlet Energy rate Qdot_in
    %--------------------------------------------------------------
    Qdot_in = params.mdot * params.c_w * abs(T_w_in_upper - T_w_in_lower);

    E = (Qdot_in * dt_mix) / (params.rho_w * params.c_w * params.A_hws);

    %--------------------------------------------------------------
    % 6) Caculate new Slope a
    %--------------------------------------------------------------
    a = (E - I_num) / G;

    %--------------------------------------------------------------
    % 7) Apply closed-form update
    %--------------------------------------------------------------
    z = (ids_mix - z_mix_idx) * params.dz;   % physical z
    Tcur(ids_mix) = Tin + a * (z - dz_star);


end