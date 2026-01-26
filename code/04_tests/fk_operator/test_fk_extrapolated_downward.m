%% TEST: extrapolated_zmix – downward FK
% ---------------------------------------------------------------
% Tests the FK operator for an EXTRAPOLATED mixing height with
% DOWNWARD free convection.
%
% Focus:
%   - apply_fk_operator
%   - fk_case = 'extrapolated_zmix'
%   - fk_extrapolated_zmix_Qscaled is exercised internally
%
% No classification, no wrapper logic.
% ---------------------------------------------------------------

run(fullfile('..', '..', '..', 'configs', 'paths_relative.m'));


%% --------------------------------------------------------------
% 1) Setup: geometry & temperature profile
% ---------------------------------------------------------------

dz = 0.005;                    % [m]
H  = 2.0;                      % [m]
Nz = round(H/dz) + 1;

z = (0:Nz-1).' * dz;

% Monotone stratified profile (warmer at top)
Tcur = 40 + 2.0*z;             % 2 °C / m

T_old = Tcur;

z_inlet_idx = round(Nz/2);     % inlet at half height
z_mix_idx   = 1;               % mixing extends downward to bottom
ids_mix     = z_inlet_idx:-1:z_mix_idx;

% Inlet temperature BELOW entire profile → extrapolated zmix
Tin = min(Tcur) - 1.0;         % [°C]

% Inlet-adjacent temperatures (only magnitude matters)
T_w_in_upper = Tin;
T_w_in_lower = Tcur(z_inlet_idx);

%% --------------------------------------------------------------
% 2) FK parameters (hard-coded, consistent)
% ---------------------------------------------------------------

params.dz    = dz;
params.rho_w = 1000;           % [kg/m^3]
params.c_w   = 4180;           % [J/(kg K)]
params.A_hws = 1.0;            % [m^2]
params.mdot  = 0.1;            % [kg/s]

dt_mix = 3600;                 % [s]

fk_case = 'extrapolated_zmix';
fklog   = @(s) [];

%% --------------------------------------------------------------
% 3) Action: apply FK operator
% ---------------------------------------------------------------

T_new = apply_fk_operator( ...
    Tcur, ids_mix, z_mix_idx, Tin, ...
    T_w_in_upper, T_w_in_lower, ...
    fk_case, dt_mix, params, fklog);

%% --------------------------------------------------------------
% 4) Calculate Energies before and after FK
% ---------------------------------------------------------------

T_ref = T_old(z_inlet_idx);
Qdot_mix = params.mdot * params.c_w * (Tin - T_ref);
Q_add    = Qdot_mix * dt_mix;

% --- cell-averaged temperatures in mixing zone ---
Tcell_old = 0.5 * (T_old(ids_mix(1:end-1)) + T_old(ids_mix(2:end)));
Tcell_new = 0.5 * (T_new(ids_mix(1:end-1)) + T_new(ids_mix(2:end)));

% --- thermal energy (integral over cells) ---
E_old = params.rho_w * params.c_w * params.A_hws * params.dz * ...
        sum(Tcell_old);

E_new = params.rho_w * params.c_w * params.A_hws * params.dz * ...
        sum(Tcell_new);


%% --------------------------------------------------------------
% 5) Plausibility checks
% ---------------------------------------------------------------


tol = 1e-6;

% --- 4.1 Lower boundary values in mix-zone ---
assert(T_new(z_inlet_idx) <= T_old(z_inlet_idx) - tol, ...
    'Inlet temperature increased for extrapolated downward FK');

assert(T_new(1) <= T_old(1) - tol, ...
    'Bottom temperature increased for extrapolated downward FK');

% --- 4.2 Monotonicity preserved ---
assert(all(diff(T_new(ids_mix)) <= -tol), ...
    'Temperature profile lost monotonicity');

% --- 4.3 Extrapolated bound check (Tin must not be undershot) ---
assert(min(T_new(ids_mix)) >= Tin - tol, ...
    'Extrapolated FK undershot inlet temperature Tin');
    
% --- 4.4 Check energy conservation ---
rel_err = abs((E_new - E_old) - Q_add) / max(abs(Q_add), tol);
assert(rel_err < tol, ...
    'Energy is not conserved in extrapolated downward FK');
    
disp('✓ TEST PASSED: extrapolated_zmix / downward FK');

%% --------------------------------------------------------------
% 6) Visual Sanity Check
% ---------------------------------------------------------------
figure; hold on; grid on;
plot(T_old, z, 'k--', 'LineWidth', 1.2);
plot(T_new, z, 'b-',  'LineWidth', 1.5);
xlabel('Temperature [°C]');
ylabel('z [m]');
legend('old profile','FK updated profile','Location','best');
title('Extrapolated FK (downward) – profile comparison');