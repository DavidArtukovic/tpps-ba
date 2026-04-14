% -------------------------------------------------------------
% Relative project paths (repo-based, portable)
% -------------------------------------------------------------

% Project root = parent of /configs
PROJECT_BASE = fileparts(mfilename('fullpath'));
PROJECT_BASE = fileparts(PROJECT_BASE);

% Code
CODE_BASE    = fullfile(PROJECT_BASE, 'code');

% Add model paths
addpath(CODE_BASE);
addpath(fullfile(CODE_BASE, '01_model', 'scenario_1_freeConv_extended'));
addpath(fullfile(CODE_BASE, '01_model', 'scenario_1_freeConv_extended', 'utils_fk'));
