%% TEST: lambda_hybrid – extrapolated_zmix / upward / violation
% ---------------------------------------------------------------
% Tests the FK/FM lambda hybrid in an EXTRAPOLATED mixing case
% with UPWARD free convection where FK violates boundary conditions.
%
% Focus:
%   - apply_fk_lambda_hybrid
%   - violation at inlet and mixing boundary
%   - hybrid must restore physical admissibility
% ---------------------------------------------------------------

run(fullfile('..', '..', '..', 'configs', 'paths_relative.m'));

%% --------------------------------------------------------------
% 1) Setup: geometry & temperature profile
% ---------------------------------------------------------------

dz = 0.005;
H  = 2.0;
Nz = round(H/dz) + 1;

z = (0:Nz-1).' * dz;

p = 6;                 % 0<p<1  (higher => more "extremely" convex)
x = z./H;

T_old = 40 + 4*(x.^p);    % T(0)=40, T(H)=44, monotone increasing, concave

z_inlet_idx = round(Nz/2);
z_mix_idx   = Nz;
ids_mix     = z_inlet_idx:1:z_mix_idx;

Tin = max(T_old)+8;
Tcur = T_old;

% Single-volume proxy temperatures
T_w_in_lower = Tcur(z_inlet_idx);
T_w_in_upper = Tin;

%% --------------------------------------------------------------
% 2) FK parameters
% ---------------------------------------------------------------

params.dz    = dz;
params.rho_w = 1000;
params.c_w   = 4180;
params.A_hws = 1.0;
params.mdot  = 0.03;

dt_mix = 1800;  % [s]
fk_case = 'extrapolated_zmix';
fklog = @(s) [];

%% --------------------------------------------------------------
% 3) Action: call  fk operator and lambda hybrid
% ---------------------------------------------------------------
T_fk = apply_fk_operator( ...
    Tcur, ids_mix, z_mix_idx, Tin, ...
    T_w_in_upper, T_w_in_lower, ...
    fk_case, dt_mix, params, fklog);

T_new = apply_fk_lambda_hybrid( ...
    T_old, T_fk, Tin, ids_mix, ...
    params, fk_case, dt_mix, fklog);


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

% --- 4.1 Boundary violations fixed by lambda hybrid ---
assert(T_new(z_inlet_idx) >= T_old(z_inlet_idx) - tol, ...
    'Lambda hybrid did not fix inlet boundary violation (upward FK)');

assert(T_new(end) >= T_old(end) - tol, ...
    'Lambda hybrid did not fix top boundary violation (upward FK)');

% --- 4.2 Monotonicity preserved in mixing zone ---
assert(all(diff(T_new(ids_mix)) >= tol), ...
    'Temperature profile lost monotonicity after lambda hybrid');

% --- 4.3 Extrapolated bound check (Tin must not be undershot) ---
assert(max(T_new(ids_mix)) <= Tin - tol, ...
    'Hybrid result overshoot inlet temperature Tin');

% --- 4.4 Check energy conservation ---
rel_err = abs(E_new - (E_old + Q_add)) / max(abs(Q_add), tol);
assert(rel_err < tol, ...
    'Energy is not conserved in lambda extrapolated upward violation FK');

disp('✓ TEST PASSED: extrapolated_zmix / upward FK + lambda hybrid');


%% --------------------------------------------------------------
% 6) Visual Sanity Check
% ---------------------------------------------------------------
figure; hold on; grid on;

plot(T_old, z, 'k--', 'LineWidth', 1.2);
plot(T_fk,  z, 'r-',  'LineWidth', 1.2);
plot(T_new, z, 'b-',  'LineWidth', 1.5);

xlabel('Temperature [°C]');
ylabel('z [m]');
legend('old profile','FK (violating)','hybrid result','Location','best');
title('Extrapolated FK (upward) – lambda hybrid repair');
