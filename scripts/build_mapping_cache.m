% PAPER_REPO_PATH_SETUP
this_script = string(mfilename('fullpath'));
if strlength(this_script) == 0
    repo_root = pwd;
else
    repo_root = fileparts(fileparts(char(this_script)));
end
cd(repo_root);
addpath(genpath(fullfile(repo_root, 'src')));

out_csv = fullfile('configs', 'pdb_uniprot_residue_map.csv');
target_files = [ ...
    string(fullfile('configs', 'exact_targets.csv')), ...
    string(fullfile('configs', 'noisy_targets.csv')), ...
    string(fullfile('configs', 'state_targets.csv'))];

T_all = table();
for k = 1:numel(target_files)
    paper_check_targets(target_files(k));
    tmp_csv = fullfile(tempdir, sprintf('paper_sifts_%d.csv', k));
    T = build_pdb_uniprot_mapping_cache(target_files(k), ...
        'OutCSV', tmp_csv, ...
        'SiftsDir', fullfile('data', 'SIFTS'), ...
        'SkipMissing', false, ...
        'Verbose', false);
    T_all = paper_append_row(T_all, T);
end

if isempty(T_all)
    error('No SIFTS mapping rows were generated.');
end

key = string(T_all.AFDB_ID) + "|" + string(T_all.PDB_ID) + "|" + ...
      string(T_all.Chain_ID) + "|" + string(T_all.PDB_resNum) + "|" + ...
      string(T_all.PDB_iCode) + "|" + string(T_all.UniProt_resNum);
[~, ia] = unique(key, 'stable');
T_all = T_all(ia, :);

writetable(T_all, out_csv);
