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

target_csv = fullfile('configs', 'noisy_targets.csv');
paper_check_targets(target_csv);

targets = readtable(target_csv, 'FileType', 'text', 'VariableNamingRule', 'preserve');
eta_values = [0.01, 0.05];
noise_seeds = 1:5;

runs = table();
for r = 1:height(targets)
    for e = 1:numel(eta_values)
        eta = eta_values(e);
        K = eta / 3;
        for s = noise_seeds
            cfg = paper_cfg_from_row(target_csv, targets, r, "AF_rank", 6, K, s);
            cfg.af_rank_jitter = 1e-3;
            [result, ~] = paper_run_cfg(cfg);
            row = paper_result_row(result);
            row.eta = eta;
            row = movevars(row, 'eta', 'After', 'pdb_id');
            runs = paper_append_row(runs, row);
        end
    end
end

writetable(runs, fullfile(out_dir, 'noisy_constraints_runs.csv'));

summary = table();
for r = 1:height(targets)
    afdb_id = paper_scalar_text(targets.("AFDB ID"), r);
    pdb_id = upper(paper_scalar_text(targets.("PDB ID"), r));
    for e = 1:numel(eta_values)
        eta = eta_values(e);
        mask = runs.afdb_id == afdb_id & abs(runs.eta - eta) < 1e-12;
        x = runs(mask, :);
        row = table( ...
            afdb_id, pdb_id, eta, ...
            mean(x.edg_rmsd_ca), std(x.edg_rmsd_ca), ...
            mean(x.edg_gdt_ts), std(x.edg_gdt_ts), ...
            mean(x.edg_gdt_ha), std(x.edg_gdt_ha), ...
            mean(x.edg_lddt), std(x.edg_lddt), ...
            'VariableNames', {'afdb_id','pdb_id','eta', ...
            'rmsd_ca_mean','rmsd_ca_sd','gdt_ts_mean','gdt_ts_sd', ...
            'gdt_ha_mean','gdt_ha_sd','edg_lddt_mean','edg_lddt_sd'});
        summary = paper_append_row(summary, row);
    end
end

writetable(summary, fullfile(out_dir, 'noisy_constraints_summary.csv'));
