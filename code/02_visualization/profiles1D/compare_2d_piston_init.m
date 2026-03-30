clc
clear
close all

%%%------------------------------------------%%%
% 01. Load paths
%%%------------------------------------------%%%
run(fullfile('..','..','..', 'configs', 'paths_local.m'));

DATA_INIT       = fullfile(DATA_BASE, 'init');
DATA_SCEN1_BASE = fullfile(DATA_BASE, 'Modellvergleich1D');

%%%------------------------------------------%%%
% 02. File names
%%%------------------------------------------%%%
file_init_ref  = '20260324_d18_hp16.0_gap0.5_2D_chunk009_v1.mat';
file_init_old  = '20260328_Init_piston_d18_hp18.0_gap0.5_Init2D_chunk009_v3.mat';
file_comsol    = 'T1_1818_900.mat';

%%%------------------------------------------%%%
% 03. Load data
%%%------------------------------------------%%%
S_ref = load(fullfile(DATA_INIT, file_init_ref));
S_old = load(fullfile(DATA_INIT, file_init_old));
S_c   = load(fullfile(DATA_SCEN1_BASE, file_comsol));

%%%------------------------------------------%%%
% 04. Extract piston fields
%%%------------------------------------------%%%
% Adjust here only if one file stores the field differently
T_ref = S_ref.InitOut.state.T_RPf_end;
T_old = S_old.InitOut.state.T_RPf_end;

z_RP  = S_old.InitOut.grid.z_RP;
Nz    = S_old.InitOut.grid.Nz;

Nz_p = Nz(3);
Nr_p = Nz(14);

% piston height coordinate [m]
z_p = 0:0.005:18;

%%%------------------------------------------%%%
% 05. Area-weighted mean piston temperature
%%%------------------------------------------%%%
areas   = z_RP(2,:);
weights = areas / sum(areas);

T_ref_mean = T_ref * weights';
T_old_mean = T_old * weights';

%%%------------------------------------------%%%
% 06. COMSOL piston profile (already mean)
%%%------------------------------------------%%%
T_comsol = S_c.P(:,1);      % first radial column (axis temperature)

dz_C = 0.05;
z_comsol = (0:length(T_comsol)-1)' * dz_C;

dz_M = 0.005;
z_matlab = (0:length(T_ref_mean)-1)' * dz_M;

T_comsol_interp = interp1(z_comsol, T_comsol, z_matlab, 'linear', 'extrap');

%%%------------------------------------------%%%
% 07. Plot
%%%------------------------------------------%%%
figure('Color','w');
hold on
grid on
box on

plot(T_ref_mean, z_p, 'k-',  'LineWidth', 1.8)
plot(T_old_mean, z_p, 'b--', 'LineWidth', 1.8)
plot(T_comsol_interp, z_matlab, 'g-', 'LineWidth', 1.8)

xlabel('Temperature [^\circC]')
ylabel('Piston height [m]')
title('Comparison of mean piston temperature profiles')
legend({'MATLAB reference', ...
        'MATLAB flux radial operator', ...
        'COMSOL'}, ...
        'Location','best')

ylim([0 18])

%%% optional:
% xlim([40 80])


%%

t0 = 0;
t1 = 2890*3600;                         % 120 days in s
t_end = 35064*3600+3600*1700;           % 1532 days in s or 4.2 years

time_short = t0:900:t1;                 % short time array in 15 minutes steps
time_long = t0:900:t_end;               % long time array in 15 minutes steps

T = zeros(3,length(time_long));         % 3-row temperature array (only third row relevant)

% Custom linear temperature profile for 120 days
for k = 1:length(time_short)
    T(3,k) = 11+273.15 + time_long(k)*(286.25-(11+273.15))/(t1);
end
% Custom temperature profile for 4.2 years, as function of sin and exponential terms.
for k = length(time_short):length(time_long)
    T(1,k) = 11+79*(1-exp(-time_long(k)/(35064*3600/3))); % an 1
    T(2,k) = 11+44*(1-exp(-time_long(k)/(35064*3600/3))); % an 1
    T(3,k) = 264.75+((T(1,k)-T(2,k))/2)*sin(time_long(k)*(2*pi)/(3600*8766)+3600*8764/4)+((T(1,k)-T(2,k))/2)+T(2,k);
end