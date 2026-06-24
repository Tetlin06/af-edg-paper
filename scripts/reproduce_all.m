% PAPER_REPO_PATH_SETUP
this_script = string(mfilename('fullpath'));
if strlength(this_script) == 0
    repo_root = pwd;
else
    repo_root = fileparts(fileparts(char(this_script)));
end
cd(repo_root);
addpath(genpath(fullfile(repo_root, 'src')));

run(fullfile('scripts', 'reproduce_exact_constraints.m'));
run(fullfile('scripts', 'reproduce_noisy_constraints.m'));
run(fullfile('scripts', 'reproduce_state_cases.m'));
