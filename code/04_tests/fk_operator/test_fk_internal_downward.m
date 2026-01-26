%% TEST: internal_zmix – downward FK
% ---------------------------------------------------------------
% Internal mixing height, downward free convection.
% Inlet temperature lies WITHIN the temperature range.
% ---------------------------------------------------------------

run(fullfile('..', '..', '..', 'configs', 'paths_relative.m'));

%% --------------------------------------------------------------
% 1) Setup: geometry & temperature profile
% ---------------------------------------------------------------

dz = 0.005;                    % [m]
H  = 2.0;                      % [m]
Nz = round(H/dz) + 1;

z = (0:Nz-1).' * dz;

% Monotone stratified profile
Tcur = 40 + 2.0*z;             % 40°C .. 44°C

T_old = Tcur;

% Inlet at half height
z_inlet_idx = round(Nz/2);

% Mixing zone BELOW inlet (downward FK)
z_mix_idx = 1;
ids_mix   = z_inlet_idx:-1:z_mix_idx;

% Internal inlet temperature (within profile range)
Tin = 40;                    % [°C]

% Single-volume proxy temperatures
T_w_in_lower = Tcur(z_inlet_idx);
T_w_in_upper = Tin;

%% --------------------------------------------------------------
% 2) FK parameters
% ---------------------------------------------------------------

params.dz    = dz;
params.rho_w = 1000;           % [kg/m^3]
params.c_w   = 4180;           % [J/(kg K)]
params.A_hws = 1.0;            % [m^2]
params.mdot  = 0.1;            % [kg/s]

dt_mix = 3600;                 % [s]

fk_case = 'internal_zmix';
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

% --- 4.1 No new extrema (internal mix!) ---
assert(min(T_new(ids_mix)) >= min(T_old(ids_mix)) - tol, ...
    'New minimum created in internal downward FK');

assert(max(T_new(ids_mix)) <= max(T_old(ids_mix)) + tol, ...
    'New maximum created in internal downward FK');

% --- 4.2 Cooling at inlet ---
assert(T_new(z_inlet_idx) <= T_old(z_inlet_idx) - tol, ...
    'Inlet temperature increased for internal downward FK');

% --- 4.3 Monotonicity preserved ---
assert(all(diff(T_new(ids_mix)) <= -tol), ...
    'Temperature profile lost monotonicity');

% --- 4.4 Check energy conservation ---
rel_err = abs((E_new - E_old) - Q_add) / max(abs(Q_add), tol);
assert(rel_err < tol, ...
    'Energy is not conserved in internal downward FK');

disp('✓ TEST PASSED: internal_zmix / downward FK');

%% --------------------------------------------------------------
% 6) Visual sanity check
% ---------------------------------------------------------------

figure; hold on; grid on;
plot(T_old, z, 'k--', 'LineWidth', 1.2);
plot(T_new, z, 'b-',  'LineWidth', 1.5);
xlabel('Temperature [°C]');
ylabel('z [m]');
legend('old profile','FK updated profile','Location','best');
title('Internal FK (downward) – profile comparison');
