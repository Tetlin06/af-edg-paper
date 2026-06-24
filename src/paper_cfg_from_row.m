function cfg = paper_cfg_from_row(target_csv, T, r, init_method, cutoff, K, noise_seed)
%PAPER_CFG_FROM_ROW Construct a paper-only cfg from one target table row.
    [opts, lsopts] = paper_default_solver_opts();

    cfg = struct();
    cfg.target_csv = char(string(target_csv));
    cfg.mapping_csv = fullfile('configs', 'pdb_uniprot_residue_map.csv');

    cfg.afdb_id = char(paper_scalar_text(T.("AFDB ID"), r));
    cfg.pdb_id = char(upper(paper_scalar_text(T.("PDB ID"), r)));
    cfg.chain_id = char(paper_scalar_text(T.("Chain ID"), r));

    cfg.domain_cut = [];
    cfg.cutoff = cutoff;
    cfg.K = K;
    cfg.noise_seed = noise_seed;

    cfg.init_method = string(init_method);
    cfg.init_seed = 47;
    cfg.af_rank_jitter = 0;

    cfg.opts = opts;
    cfg.lsopts = lsopts;
end
