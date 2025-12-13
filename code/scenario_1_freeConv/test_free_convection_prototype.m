%% test_free_convection_prototype.m
% Prototype testbench for the free convection module.
% - Creates 34 temperature profiles to trigger all FK branches.
% - Applies one dt_mix update: T_after = T_before + dTdt_free * dt_mix
% - Plots before/after and delta for each profile.

clear; clc; close all;

%% ------------------------ Parameters ----------------------------
% Grid / geometry
Nw   = 60;           % number of water nodes
dz   = 0.25;         % [m]
z    = (0:Nw-1)' * dz;

% Indices defining water sub-volumes & inlets (adapt as needed)
z_w_lower_start_idx = 1;
z_w_lower_end_idx   = round(0.45*Nw);      % lower volume ends ~45%
z_w_upper_start_idx = z_w_lower_end_idx+1;
z_w_upper_end_idx   = Nw;

z_inlet_lower_idx   = round(0.25*Nw);      % inlet somewhere in lower volume
z_inlet_upper_idx   = round(0.75*Nw);      % inlet somewhere in upper volume

% Physical parameters (order-of-magnitude; enough for prototype)
rho_w = 997;         % [kg/m^3]
c_w   = 4180;        % [J/(kg*K)]
A_hws = 0.5;         % [m^2] cross section
flow_abs = 2e-4;     % [m/s] prototype flow magnitude
dt_mix = 900;        % [s] Schäfer fixed FK timestep

% Derived
v_flow = abs(flow_abs);
Vdot   = v_flow * A_hws;      % [m^3/s]
mdot   = rho_w * Vdot;        % [kg/s]

%% ---------------------- Build 34 profiles -----------------------
% We want to trigger all four FK directions:
% A) charging (flow<0)  : upward   when T_upper_in > T_lower_in
% B) charging (flow<0)  : downward when T_upper_in < T_lower_in
% C) discharging (flow>0): upward  when T_lower_in > T_upper_in
% D) discharging (flow>0): downward when T_lower_in < T_upper_in
%
% We create 34 profiles: 8+8+8+8 = 32 plus 2 edge-ish cases.

profiles = cell(34,1);
labels   = strings(34,1);
flows    = zeros(34,1);

k = 1;

% Helper to create smooth-ish profiles
make_profile = @(T0, grad, bumpAmp, bumpCenter, bumpWidth) ...
    (T0 + grad*z + bumpAmp*exp(-((z-bumpCenter)/bumpWidth).^2));

% --- Group A: charging upward (flow<0, T_upper_in > T_lower_in)
flowsA = -flow_abs;
for i = 1:1
    T0 = 30 + 2*i;
    grad = -0.25; % cooler with height (so upper tends to be cooler unless bump)
    bumpAmp = 10 + i;                 % make upper inlet hotter
    bumpCenter = z(z_inlet_upper_idx);
    bumpWidth  = 0.6 + 0.05*i;
    T = make_profile(T0, grad, bumpAmp, bumpCenter, bumpWidth);
    profiles{k} = T;
    labels(k)   = "A charging/upward #" + i;
    flows(k)    = flowsA;
    k = k+1;
end

% --- Group B: charging downward (flow<0, T_upper_in < T_lower_in)
flowsB = -flow_abs;
for i = 1:1
    T0 = 55 - i;
    grad = +0.15; % warmer with height (helps make upper inlet warmer, but we want upper < lower at inlets)
    bumpAmp = 12 + i;                 % make LOWER inlet hotter
    bumpCenter = z(z_inlet_lower_idx);
    bumpWidth  = 0.7;
    T = make_profile(T0, grad, bumpAmp, bumpCenter, bumpWidth);
    profiles{k} = T;
    labels(k)   = "B charging/downward #" + i;
    flows(k)    = flowsB;
    k = k+1;
end

% --- Group C: discharging upward (flow>0, T_lower_in > T_upper_in)
flowsC = +flow_abs;
for i = 1:1
    T0 = 40 + i;
    grad = +0.10;
    bumpAmp = 14 + i;                 % make LOWER inlet hotter than upper inlet
    bumpCenter = z(z_inlet_lower_idx);
    bumpWidth  = 0.5 + 0.05*i;
    T = make_profile(T0, grad, bumpAmp, bumpCenter, bumpWidth);
    profiles{k} = T;
    labels(k)   = "C discharging/upward #" + i;
    flows(k)    = flowsC;
    k = k+1;
end

% --- Group D: discharging downward (flow>0, T_lower_in < T_upper_in)
flowsD = +flow_abs;
for i = 1:1
    T0 = 35 + 2*i;
    grad = -0.10;
    bumpAmp = 16 + i;                 % make UPPER inlet hotter than lower inlet
    bumpCenter = z(z_inlet_upper_idx);
    bumpWidth  = 0.55;
    T = make_profile(T0, grad, bumpAmp, bumpCenter, bumpWidth);
    profiles{k} = T;
    labels(k)   = "D discharging/downward #" + i;
    flows(k)    = flowsD;
    k = k+1;
end

% % --- Two extra “edge-ish” cases (still FK on but minimal effect)
% % 33: almost flat, tiny bump at upper inlet (charging)
% profiles{33} = make_profile(50, 0.0, 0.8, z(z_inlet_upper_idx), 0.8);
% labels(33)   = "Edge charging (weak FK)";
% flows(33)    = -flow_abs;
% 
% % 34: almost flat, tiny bump at lower inlet (discharging)
% profiles{34} = make_profile(50, 0.0, 0.8, z(z_inlet_lower_idx), 0.8);
% labels(34)   = "Edge discharging (weak FK)";
% flows(34)    = +flow_abs;

%% -------------------- Run FK and plot (2 cases per figure) ---------------------------
nCases = 34;
casesPerFig = 2;

for startId = 1:casesPerFig:nCases
    endId = min(startId + casesPerFig - 1, nCases);

    figure('Name', sprintf('Free convection prototype (cases %d-%d)', startId, endId));
    tl = tiledlayout((endId-startId+1), 2, 'TileSpacing','compact', 'Padding','compact');

    for caseId = startId:endId
        T_before = profiles{caseId};
        flow     = flows(caseId);

        % Apply FK module (one dt_mix step)
        T_after = apply_free_convection_step( ...
            T_before, flow, dz, dt_mix, ...
            z_inlet_upper_idx, z_inlet_lower_idx, ...
            z_w_lower_start_idx, z_w_lower_end_idx, ...
            z_w_upper_start_idx, z_w_upper_end_idx, ...
            mdot, c_w, A_hws, rho_w);

        % --- left: before vs after
        nexttile;
        plot(T_before, z, '-', T_after, z, '--', 'LineWidth', 1.5);
        grid on; ylabel('z [m]'); xlabel('T [°C]');
        title(labels(caseId), 'Interpreter','none');
        legend('before','after','Location','best');

        % --- right: delta
        nexttile;
        plot(T_after - T_before, z, '-', 'LineWidth', 1.5);
        grid on; ylabel('z [m]'); xlabel('\DeltaT [K]');
        title('after - before');
        xline(0,'--'); % helps visually
    end

    title(tl, sprintf('Free convection prototype (cases %d-%d)', startId, endId));
end

%% ===================== Local function ===========================
function T_new = apply_free_convection_step( ...
    T, flow, dz, dt_mix, ...
    z_inlet_upper_idx, z_inlet_lower_idx, ...
    z_w_lower_start_idx, z_w_lower_end_idx, ...
    z_w_upper_start_idx, z_w_upper_end_idx, ...
    mdot, c_w, A_hws, rho_w)
% Apply one "fixed" dt_mix FK update step based on the provided snippet.
% Output: updated temperature vector after one FK step.

    % Initialize derivative contribution
    dTdt_free = zeros(size(T));

    % Store inlet temperatures
    T_w_in_upper = T(z_inlet_upper_idx);
    T_w_in_lower = T(z_inlet_lower_idx);

    % -------------------- 8.1 charging (flow < 0) --------------------
    if flow < 0

        z_mix_idx = z_inlet_lower_idx; % start at lower inlet

        if T_w_in_upper > T_w_in_lower
            % upward direction in lower volume
            while (T_w_in_upper > T(z_mix_idx) && z_mix_idx < z_w_lower_end_idx)
                z_mix_idx = z_mix_idx + 1;
            end

            ids_mix = z_inlet_lower_idx:z_mix_idx;
            T_mix = T(ids_mix);
            I_mix = sum((T_w_in_upper - T_mix) * dz); % [K*m]

        else
            % downward direction in lower volume
            while (T_w_in_upper < T(z_mix_idx) && z_mix_idx > z_w_lower_start_idx)
                z_mix_idx = z_mix_idx - 1;
            end

            ids_mix = z_inlet_lower_idx:-1:z_mix_idx;
            T_mix = T(ids_mix);
            I_mix = sum((T_w_in_upper - T_mix) * dz * (-1)); % [K*m]
        end

        Qdot = mdot * c_w * (T_w_in_upper - T_w_in_lower); % [W]

        dz_mix = max((numel(ids_mix)-1) * dz, eps); % [m]
        z_dist = (ids_mix - z_mix_idx).' * dz;       % [m]
        geom_factor = 2 * z_dist / (dz_mix^2);       % [1/m]

        dTdt_free(ids_mix) = (1/dt_mix) * (geom_factor * I_mix + (T_w_in_upper - T_mix)) ...
                             - Qdot * geom_factor / (A_hws * rho_w * c_w);

    % ------------------ 8.2 discharging (flow > 0) -------------------
    elseif flow > 0

        z_mix_idx = z_inlet_upper_idx; % start at upper inlet

        if T_w_in_lower > T_w_in_upper
            % upward direction in upper volume
            while (T_w_in_lower > T(z_mix_idx) && z_mix_idx < z_w_upper_end_idx)
                z_mix_idx = z_mix_idx + 1;
            end

            ids_mix = z_inlet_upper_idx:z_mix_idx;
            T_mix = T(ids_mix);
            I_mix = sum((T_w_in_lower - T_mix) * dz); % [K*m]

        else
            % downward direction in upper volume
            while (T_w_in_lower < T(z_mix_idx) && z_mix_idx > z_w_upper_start_idx)
                z_mix_idx = z_mix_idx - 1;
            end

            % NOTE: fixed to descending index order (otherwise empty range)
            ids_mix = z_inlet_upper_idx:-1:z_mix_idx;
            T_mix = T(ids_mix);
            I_mix = sum((T_w_in_lower - T_mix) * dz * (-1)); % [K*m]
        end

        Qdot = mdot * c_w * (T_w_in_lower - T_w_in_upper); % [W]

        dz_mix = max((numel(ids_mix)-1) * dz, eps); % [m]
        z_dist = (ids_mix - z_mix_idx).' * dz;       % [m]

        % NOTE: fixed (your snippet used z_loc which was undefined)
        geom_factor = 2 * z_dist / (dz_mix^2);       % [1/m]

        dTdt_free(ids_mix) = (1/dt_mix) * (geom_factor * I_mix + (T_w_in_lower - T_mix)) ...
                             - Qdot * geom_factor / (A_hws * rho_w * c_w);
    end

    % Apply one dt_mix step update
    T_new = T + dTdt_free * dt_mix;
end
