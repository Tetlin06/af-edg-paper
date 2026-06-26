% BUILD_MAPPING_CACHE
%
% Build the local PDB-to-UniProt residue mapping cache from the committed
% target sheets and local SIFTS XML files.

close all force;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(fullfile(repoRoot, 'src')));

targetFiles = [
    "configs/exact_targets.csv"
    "configs/noisy_targets.csv"
    "configs/state_targets.csv"
];

T_all = table();

for k = 1:numel(targetFiles)
    f = targetFiles(k);

    if ~isfile(f)
        continue;
    end

    T = readtable(f, 'FileType', 'text', 'VariableNamingRule', 'preserve');

    required = {'AFDB ID', 'PDB ID', 'Chain ID'};
    for r = 1:numel(required)
        if ~any(strcmp(required{r}, T.Properties.VariableNames))
            error('Target file %s is missing required column: %s', f, required{r});
        end
    end

    T = T(:, required);
    T_all = [T_all; T]; %#ok<AGROW>
end

if isempty(T_all)
    error('No target rows found.');
end

[~, ia] = unique(strcat(string(T_all.("AFDB ID")), "|", string(T_all.("PDB ID")), "|", string(T_all.("Chain ID"))), 'stable');
T_all = T_all(ia, :);

tmpCSV = fullfile('configs', '_mapping_targets_tmp.csv');
writetable(T_all, tmpCSV);
cleanupObj = onCleanup(@() delete_if_exists(tmpCSV)); %#ok<NASGU>

build_pdb_uniprot_mapping_cache( ...
    tmpCSV, ...
    'SiftsDir', fullfile('data', 'SIFTS'), ...
    'OutCSV', fullfile('configs', 'pdb_uniprot_residue_map.csv'), ...
    'SkipMissing', false, ...
    'Verbose', false);

function delete_if_exists(path)
    if isfile(path)
        delete(path);
    end
end
