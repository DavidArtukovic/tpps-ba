%% TEST: internal_zmix – downward FULLY MIXED
% ---------------------------------------------------------------
% Internal mixing height, downward fully mixed case.
% Entire mixing zone is shifted downward, cooler water enters.
% ---------------------------------------------------------------

run(fullfile('..', '..', '..', 'configs', 'paths_relative.m'));

%% --------------------------------------------------------------
% 1) Setup: geometry & temperature profile
% --------------------------------------------------------------

dz = 0.005;                    % [m]
H  = 2.0;                      % [m]
Nz = round(H/dz) + 1;

z = (0:Nz-1).' * dz;

% Monotone stratified profile (same as FK test)
Tcur = 40 + 2.0*z;             % 40°C .. 44°C
T_old = Tcur;

% Inlet at half height
z_inlet_idx = round(Nz/2);

% Mixing zone BELOW inlet (downward)
z_mix_idx = 1;
ids_mix   = z_inlet_idx:-1:z_mix_idx;

% Cooler inlet water
Tin = 40;                      % [°C]

% Proxy inlet temperatures
T_w_in_lower = Tin;
T_w_in_upper = Tcur(z_inlet_idx);

%% --------------------------------------------------------------
% 2) Parameters
% --------------------------------------------------------------

params.dz    = dz;
params.rho_w = 1000;           % [kg/m^3]
params.c_w   = 4180;           % [J/(kg K)]
params.A_hws = 1.0;            % [m^2]
params.mdot  = 0.1;            % [kg/s]

dt_mix = 3600;                 % [s]

fk_case = 'fully_mixed';
fklog   = @(s) [];

%% --------------------------------------------------------------
% 3) Action: apply fully mixed operator
% --------------------------------------------------------------

T_new = apply_fk_operator( ...
    Tcur, ids_mix, z_mix_idx, Tin, ...
    T_w_in_upper, T_w_in_lower, ...
    fk_case, dt_mix, params, fklog);

%% --------------------------------------------------------------
% 4) Energy balance
% --------------------------------------------------------------

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
% --------------------------------------------------------------

tol = 1e-6;


% --- 5.1 Cooling anywhere in mixing zone ---
assert(all(T_new(ids_mix) <= T_old(ids_mix) + tol), ...
    'Heating occurred in fully mixed downward case');


% --- 5.2 Energy conservation ---
rel_err = abs((E_new - E_old) - Q_add) / max(abs(Q_add), tol);

assert(rel_err < tol, ...
    'Energy is not conserved in downward fully mixed test');

disp('✓ TEST PASSED: downward FK - FULLY MIXED');

%% --------------------------------------------------------------
% 6) Visual sanity check
% --------------------------------------------------------------

figure; hold on; grid on;
plot(T_old, z, 'k--', 'LineWidth', 1.2);
plot(T_new, z, 'r-',  'LineWidth', 1.5);
xlabel('Temperature [°C]');
ylabel('z [m]');
legend('old profile','fully mixed profile','Location','best');
title('Internal downward fully mixed – profile comparison');
