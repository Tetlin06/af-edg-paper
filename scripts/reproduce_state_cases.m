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

target_csv = fullfile('configs', 'state_targets.csv');
paper_check_targets(target_csv);

targets = readtable(target_csv, 'FileType', 'text', 'VariableNamingRule', 'preserve');

rows = table();
for r = 1:height(targets)
    cfg = paper_cfg_from_row(target_csv, targets, r, "AF_3D", 5, 0, NaN);
    [result, ~] = paper_run_cfg(cfg);
    row = paper_result_row(result);
    row.case_name = paper_scalar_text(targets.("Case"), r);
    row.target_state = paper_scalar_text(targets.("Target state"), r);
    row = movevars(row, {'case_name','target_state'}, 'Before', 'afdb_id');
    rows = paper_append_row(rows, row);
end

writetable(rows, fullfile(out_dir, 'state_cases.csv'));
