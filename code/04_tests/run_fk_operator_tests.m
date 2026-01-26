%% Run FK operator tests
% ---------------------------------------------------------------
% Executes all free-convection operator tests.
% Stops immediately if a test fails.
%
% Intended for:
%   - local development
%   - BA documentation consistency
% ---------------------------------------------------------------

clear; clc;

%% Add model paths
run(fullfile('..', '..','configs','paths_relative.m'));

this_file = mfilename('fullpath');
this_dir  = fileparts(this_file);

fprintf('============================================\n');
fprintf(' Running FK operator test suite\n');
fprintf('============================================\n\n');

% FK operator test directory
test_dir = fullfile(this_dir, 'fk_operator');
assert(isfolder(test_dir), ...
    'Test directory not found: %s', test_dir);

test_files = {
    'test_fk_uniform_shift_upward'
    'test_fk_uniform_shift_downward'
    'test_fk_internal_upward.m'
    'test_fk_internal_downward.m'
    'test_fk_extrapolated_upward.m'
    'test_fk_extrapolated_downward.m'
    'test_fk_lambda_extrapolated_downward_no_violation.m'
    'test_fk_lambda_extrapolated_downward_violation.m'
    'test_fk_lambda_extrapolated_upward_violation.m'
};

n_fail = 0;
%%
for k = 1:numel(test_files)

    test_name = test_files{k};
    fprintf('→ Running %s ...\n', test_name);

    try
        run(fullfile(test_dir, test_name));
        fprintf('  ✓ PASSED\n\n');
    catch ME
        fprintf('  ✗ FAILED\n');
        fprintf('  %s\n\n', ME.message);
        n_fail = n_fail + 1;
    end
end

fprintf('--------------------------------------------\n');

if n_fail == 0
    fprintf(' All FK operator tests PASSED successfully.\n');
else
    fprintf(' Test suite aborted due to failure.\n');
end

fprintf('--------------------------------------------\n');
