function paper_check_targets(target_csv)
    T = readtable(target_csv, 'FileType', 'text', 'VariableNamingRule', 'preserve');

    required = {"AFDB ID", "PDB ID", "Chain ID"};
    for k = 1:numel(required)
        if ~any(strcmp(required{k}, T.Properties.VariableNames))
            error('Target CSV is missing required column: %s', required{k});
        end
    end

    chains = strtrim(string(T.("Chain ID")));
    bad = chains == "" | contains(upper(chains), "TODO");

    if any(bad)
        rows = find(bad);
        error('Target CSV %s has missing/TODO chain IDs at row(s): %s', ...
            target_csv, strjoin(string(rows.'), ', '));
    end
end
