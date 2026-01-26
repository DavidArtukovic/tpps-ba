%% TEST: lambda_hybrid – extrapolated_zmix / downward / no violation
% ---------------------------------------------------------------
% Tests the FK/FM lambda hybrid in an EXTRAPOLATED mixing case
% with DOWNWARD free convection where NO boundary violation occurs.
%
% Focus:
%   - apply_fk_lambda_hybrid
%   - fk_case = 'extrapolated_zmix'
%   - FK result should be accepted unchanged (lambda = 1)
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
T_old = 40 + 2.0*z;            % 2 °C / m

z_inlet_idx = round(Nz/2);
z_mix_idx   = 1;
ids_mix     = z_inlet_idx:-1:z_mix_idx;

Tin = min(T_old) - 1.0;        % extrapolated (cold inlet)


Tcur = T_old;
% Single-volume proxy temperatures
T_w_in_lower = T_old(z_inlet_idx);
T_w_in_upper = Tin;

%% --------------------------------------------------------------
% 2) FK parameters (hard-coded)
% ---------------------------------------------------------------

params.dz    = dz;
params.rho_w = 1000;
params.c_w   = 4180;
params.A_hws = 1.0;
params.mdot  = 0.1;

dt_mix = 3600;
fk_case = 'extrapolated_zmix';
fklog = @(s) [];

%% --------------------------------------------------------------
% 3) Action: apply lambda hybrid
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

% --- 4.1 No modification expected ---
assert(norm(T_new - T_fk) < tol, ...
    'Lambda hybrid modified FK result although no violation occurred');

% --- 4.2 Check energy conservation ---
rel_err = abs((E_new - E_old) - Q_add) / max(abs(Q_add), tol);
assert(rel_err < tol, ...
    'Energy is not conserved in lambda extrapolate downward no violation FK');

disp('✓ TEST PASSED: lambda hybrid / extrapolated_zmix / no violation');

%% --------------------------------------------------------------
% 6) Visual Sanity Check
% ---------------------------------------------------------------
figure; hold on; grid on;
plot(T_old, z, 'k--', 'LineWidth', 1.2);
plot(T_fk,  z, 'b-',  'LineWidth', 1.2);
plot(T_new, z, 'r:',  'LineWidth', 1.5);
xlabel('Temperature [°C]');
ylabel('z [m]');
legend('old profile','FK result','hybrid result','Location','best');
title('Lambda hybrid (no violation) – extrapolated downward');
