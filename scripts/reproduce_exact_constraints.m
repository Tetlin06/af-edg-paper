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
init_methods = ["Random_rand", "Random_randn", "Floyd", "AF_3D", "AF_rank"];

rows = table();
for r = 1:height(targets)
    for im = 1:numel(init_methods)
        cfg = paper_cfg_from_row(target_csv, targets, r, init_methods(im), 5, 0, NaN);
        cfg.af_rank_jitter = 0;
        [result, ~] = paper_run_cfg(cfg);
        rows = paper_append_row(rows, paper_result_row(result));
    end
end

writetable(rows, fullfile(out_dir, 'exact_constraints_all_inits.csv'));
writetable(rows(rows.init_method == "AF_3D", :), fullfile(out_dir, 'exact_constraints_paper_table.csv'));
