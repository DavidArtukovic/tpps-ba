%% TIMPESTEPS.M
% ---------------------------------------------------------------
% PURPOSE:
%   Generates hydraulic and thermal load profiles for hourly and
%   15-minute resolution simulations of the thermo-hydraulic system.
%   The profile spans over 20 days
%   The resulting profiles are stored in 'SzenarioComsol.mat' and
%   used as boundary conditions for subsequent COMSOL or MATLAB-based
%   simulation workflows.
%
% OUPTUTS:
%   - t_hour : [2 x 480] matrix
%       Row 1 → thermal charging/discharging schedule (hourly)
%       Row 2 → hydraulic piston movement schedule (hourly)
%
%   - t_900 : [2 x 1920] matrix
%       Row 1 → thermal charging/discharging schedule (15-min)
%       Row 2 → hydraulic piston movement schedule (15-min)
%
% STRUCUTRE:
%   1. Initialize arrays for hour-based and 15-min-based time grids
%   2. Construct hydraulic profile:
%        - Potential discharge phases
%        - Holding phases (stillstand)
%        - Charging phases
%   3. Construct thermal profile:
%        - Periodic charging windows
%        - Terminal discharging phase
%   4. Save profiles to a .mat file for later simulation use
%
% NOTES:
%   All profile ramps (charging/discharging) are linear.

%
% ---------------------------------------------------------------

clc
clear

%%%------------------------------------------%%%
% 01 Create 20 days charge/discharge scenario
%%%------------------------------------------%%%

% t_900(1,i) = 1 -> thermal charge
% t_900(1,i) = -1 -> thermal discharge
% t_900(1,i) = 0 -> halting

%%% Untersuchungszeitraum %%
filenameSIM = 'SzenarioComsol.mat';
t_hour = ones(2,480); % Zeitschritte in Stunde
t_900 = ones(2,480*4); % Zeitschritte in 900s (15min)

% Note first 10 days is halting phase

% Potentielles Entladen Stunden
for i = 259:259+15
t_hour(2,i)=t_hour(2,i-1)-1/16;
end

% Potenteiller Stillstand nach Entladen Stunden
t_hour(2,275:298) = 0;

% Potenteilles Laden Viertelstunden
for i = 283:283+15 
t_hour(2,i)=t_hour(2,i-1)+1/16;
end

% Potenteilles Entladen Viertelstunden
for i = 258*4+1:258*4+16*4
t_900(2,i)=t_900(2,i-1)-1/(16*4);
end

% Potenteiller Stillstand nach Entladen Viertelstunden
t_900(2,274*4:298*4) = 0;
 
% Potenteilles Laden Viertelstunden
for i = 282*4+1:283*4+15*4
t_900(2,i)=t_900(2,i-1)+1/(16*4);
end

% thermisches Laden Stunden
t_hour(1,:)=0;
for i = 0:9
t_hour(1,7+24*i:7+11+24*i) = 1;
end
% thermisches Entladen Stunden
t_hour(1,323:end) = -1;

% thermisches Laden Viertelstunden
t_900(1,:)=0;

for i = 0:9
t_900(1,1+6*4+24*4*i:1+6*4+12*4-1+24*4*i) = 1;
end

% thermisches Entladen Viertelstunden
t_900(1,322*4+1:end) = -1;

save(filenameSIM,"t_900","t_hour");