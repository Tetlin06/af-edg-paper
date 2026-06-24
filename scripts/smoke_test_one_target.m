% PAPER_REPO_PATH_SETUP
this_script = string(mfilename('fullpath'));
if strlength(this_script) == 0
    repo_root = pwd;
else
    repo_root = fileparts(fileparts(char(this_script)));
end
cd(repo_root);
addpath(genpath(fullfile(repo_root, 'src')));

out_dir = fullfile('results', 'generated');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

target_csv = fullfile('configs', 'exact_targets.csv');
paper_check_targets(target_csv);

targets = readtable(target_csv, 'FileType', 'text', 'VariableNamingRule', 'preserve');
cfg = paper_cfg_from_row(target_csv, targets, 1, "AF_3D", 5, 0, NaN);
cfg.opts.maxit = 2;
[result, ~] = paper_run_cfg(cfg);
writetable(paper_result_row(result), fullfile(out_dir, 'smoke_test_one_target.csv'));
