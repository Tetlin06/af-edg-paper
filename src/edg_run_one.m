function result = edg_run_one(cfg)
%EDG_RUN_ONE Paper-only all-atom AlphaFold-initialized EDG reconstruction.
%
% This public runner intentionally does not support:
%   subsampling, clustering, added true edges, added AlphaFold edges,
%   PAE/pLDDT edge selection, or PyMOL export.

    cfg = validate_paper_cfg(cfg);

    if exist('src', 'dir')
        addpath(genpath('src'));
    end

    truePdbPath = fullfile('data', 'PDB', char(string(cfg.pdb_id) + ".pdb"));
    afPdbPath   = fullfile('data', 'AFDB', char(string(cfg.afdb_id) + ".pdb"));

    if ~isfile(truePdbPath)
        error('Experimental PDB file not found: %s', truePdbPath);
    end
    if ~isfile(afPdbPath)
        error('AlphaFold PDB file not found: %s', afPdbPath);
    end

    out_align = align_true_vs_af_by_uniprot_mapping( ...
        cfg.target_csv, cfg.afdb_id, ...
        'MappingCSV', cfg.mapping_csv, ...
        'BuildCacheIfMissing', true, ...
        'Verbose', false);

    sub = subset_aligned_pairs(out_align, cfg.domain_cut);

    allpairs = build_aligned_all_atom_pairs( ...
        sub, truePdbPath, cfg.chain_id, afPdbPath, cfg.af_chain_id, ...
        'AtomSelection', cfg.atom_selection, ...
        'Verbose', false);

    [DistSq, Weight, graph_info] = build_allatom_dist_weight_matrix( ...
        allpairs.coords_true_all, cfg.cutoff, ...
        'AtomMeta', allpairs.atom_meta, ...
        'MinSeqSep', cfg.min_seq_sep, ...
        'Verbose', false);

    noise_info = struct();
    if cfg.K > 0
        [DistSq, noise_info] = apply_chemistry_aware_noise( ...
            DistSq, Weight, cfg.K, cfg.noise_seed, ...
            "all_atom", allpairs.atom_meta, cfg.protect_chemistry, ...
            'Verbose', false);
    end

    Pinit3 = allpairs.coords_af_all - mean(allpairs.coords_af_all, 1);
    P0 = make_paper_initialization(cfg, Pinit3, DistSq, Weight);

    if cfg.K > 0
        [GCor, IPM_Recon, solver_output] = alternating_completion_noisy( ...
            DistSq, Weight, P0, cfg.opts, cfg.lsopts);
    else
        [GCor, IPM_Recon, solver_output] = alternating_completion( ...
            DistSq, Weight, P0, cfg.opts, cfg.lsopts);
    end

    ca_rows = allpairs.ca_atom_rows(:);
    if isempty(ca_rows)
        error('No matched C-alpha rows were found after all-atom matching.');
    end

    predCA = GCor(ca_rows, :);
    trueCA = allpairs.coords_true_all(ca_rows, :);
    afCA   = allpairs.coords_af_all(ca_rows, :);

    metrics_edg = evaluate_procrustes_metrics(predCA, trueCA);
    metrics_af  = evaluate_procrustes_metrics(afCA, trueCA);
    lddt_edg    = edg_lddt_score(metrics_edg.aligned, trueCA);

    result = struct();
    result.afdb_id = string(cfg.afdb_id);
    result.pdb_id = string(cfg.pdb_id);
    result.chain_id = string(cfg.chain_id);
    result.cutoff = cfg.cutoff;
    result.K = cfg.K;
    result.noise_seed = cfg.noise_seed;
    result.init_method = string(cfg.init_method);
    result.init_seed = cfg.init_seed;
    result.af_rank_jitter = cfg.af_rank_jitter;

    result.n_atoms = size(DistSq, 1);
    result.n_ca = numel(ca_rows);
    result.n_edges = nnz(triu(Weight, 1));
    result.edge_density_upper = result.n_edges / (result.n_atoms * (result.n_atoms - 1) / 2);

    result.af_rmsd_ca = metrics_af.rmsd;
    result.edg_rmsd_ca = metrics_edg.rmsd;
    result.edg_gdt_ts = metrics_edg.gdt_ts;
    result.edg_gdt_ha = metrics_edg.gdt_ha;
    result.edg_lddt = lddt_edg.global_score100;

    result.numit = solver_output.numit;
    result.recon_error = solver_output.ReconError;
    result.graph_info = graph_info;
    result.noise_info = noise_info;
    result.solver_output = solver_output;
    result.IPM_Recon = IPM_Recon;
end
