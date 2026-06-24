function cfg = validate_paper_cfg(cfg)
%VALIDATE_PAPER_CFG Validate the intentionally narrow paper-only config.
    required = ["target_csv", "afdb_id", "pdb_id", "chain_id", ...
                "domain_cut", "cutoff", "K", "noise_seed", ...
                "init_method", "init_seed", "af_rank_jitter", ...
                "opts", "lsopts", "mapping_csv"];

    for k = 1:numel(required)
        if ~isfield(cfg, required(k))
            error('Missing required cfg field: %s', required(k));
        end
    end

    banned = ["subsample", "edge_keep", "edgekeep", "add_edges", ...
              "num_to_add", "numaf", "pae", "plddt", ...
              "cluster", "spectral", "laplacian", "pymol", "export"];

    names = string(fieldnames(cfg));
    low = lower(names);
    for b = banned
        hit = contains(low, lower(b));
        if any(hit)
            error('Out-of-scope field found in cfg: %s', names(find(hit, 1)));
        end
    end

    allowed = ["Random_rand", "Random_randn", "Floyd", "AF_3D", "AF_rank"];
    if ~any(string(cfg.init_method) == allowed)
        error('init_method must be one of: %s', strjoin(allowed, ', '));
    end

    if ~isscalar(cfg.cutoff) || ~isnumeric(cfg.cutoff) || ~(cfg.cutoff == 5 || cfg.cutoff == 6)
        error('Paper-only cutoff must be 5 or 6 Angstrom.');
    end

    if ~isscalar(cfg.K) || ~isnumeric(cfg.K) || ~isfinite(cfg.K) || cfg.K < 0
        error('K must be a nonnegative finite scalar.');
    end

    cfg.af_chain_id = 'A';
    cfg.atom_selection = "safe";
    cfg.min_seq_sep = 0;
    cfg.protect_chemistry = true;
end
