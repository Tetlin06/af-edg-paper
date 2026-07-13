% RUN_NOISY_SWEEP
% -------------------------------------------------------------------------
% Paper-only noisy AF-rank sweep.
%
% Run from the repository root with:
%   run('scripts/run_noisy_sweep.m')
%
%
%   outRoot/edgekeep_100/cutoff_06/<AFDB_ID>/K_0p00333/
%       noise_seed_003/shared/
%       noise_seed_003/runs/init_AF_rank_jitter1e-3/
%
%   outRoot/edgekeep_100/cutoff_06/<AFDB_ID>/K_0p0167/
%       noise_seed_003/shared/
%       noise_seed_003/runs/init_AF_rank_jitter1e-3/
%
% The remaining seeds use noise_seed_021, noise_seed_450,
% noise_seed_666, and noise_seed_987.

close all force;

% =========================================================================
% FIXED PAPER SETTINGS
% =========================================================================

targetCSV = fullfile('configs', 'noisy_targets.csv');
outRoot   = fullfile('out', 'init_cutoff_sweep');

trueCutoff = 6;

% Nominal noise levels reported in the manuscript. The implementation uses
% K = eta / 3 when perturbing squared distances.
nominalNoiseLevels = [0.01, 0.05];
KValues = nominalNoiseLevels / 3;

% Exact five noise seeds used for every protein and noise level.
noiseSeeds = [3, 21, 450, 666, 987];

% Exact AF-rank initialization used for every noisy run.
afRankSeed   = 47;
afRankJitter = 1e-3;

opts = struct();
opts.r           = 10000;
opts.printenergy = 1;
opts.printerror  = 1;
opts.rank        = 10;
opts.maxit       = 10000;
opts.tol         = 1e-5;
opts.lamda       = opts.r;

lsopts = struct();
lsopts.maxit = 30;
lsopts.xtol  = 1e-8;
lsopts.gtol  = 1e-8;
lsopts.ftol  = 1e-10;
lsopts.alpha = 1e-3;
lsopts.rho   = 1e-4;
lsopts.sigma = 0.1;
lsopts.eta   = 0.8;

% =========================================================================
% RUN
% =========================================================================

manifest = run_noisy_sweep_impl( ...
    targetCSV, ...
    outRoot, ...
    trueCutoff, ...
    nominalNoiseLevels, ...
    KValues, ...
    noiseSeeds, ...
    afRankSeed, ...
    afRankJitter, ...
    opts, ...
    lsopts);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function manifest = run_noisy_sweep_impl( ...
    targetCSV, outRoot, trueCutoff, nominalNoiseLevels, KValues, ...
    noiseSeeds, afRankSeed, afRankJitter, opts, lsopts)

    repoRoot = find_repo_root();
    cd(repoRoot);
    addpath(genpath(fullfile(repoRoot, 'src')));

    targetCSV = char(string(targetCSV));
    outRoot   = char(string(outRoot));
    mappingCSV = fullfile('configs', 'pdb_uniprot_residue_map.csv');

    validate_inputs( ...
        targetCSV, mappingCSV, trueCutoff, nominalNoiseLevels, KValues, ...
        noiseSeeds, afRankSeed, afRankJitter, opts);

    ensure_dir(outRoot);

    % Save one immutable description of the complete paper sweep at outRoot,
    % matching the top-level convention used by run_init_cutoff_sweep.
    cfg = struct();
    cfg.TargetCSV = targetCSV;
    cfg.OutRoot = outRoot;
    cfg.EdgeKeepFraction = 1;
    cfg.EdgeKeepFolder = "edgekeep_100";
    cfg.TrueCutoff = trueCutoff;
    cfg.NominalNoiseLevels = nominalNoiseLevels(:).';
    cfg.KValues = KValues(:).';
    cfg.NoiseSeeds = noiseSeeds(:).';
    cfg.EdgNodeLevel = "all_atom";
    cfg.MappingMethod = "uniprot";
    cfg.AllatomMinSeqSep = 0;
    cfg.AllatomAtomSelection = "safe";
    cfg.ProtectChemistryFromNoise = true;
    cfg.AFRankSeed = afRankSeed;
    cfg.AFRankJitter = afRankJitter;
    cfg.Solver = "alternating_completion_noisy";
    cfg.opts = opts;
    cfg.lsopts = lsopts;
    save(fullfile(outRoot, 'sweep_config.mat'), 'cfg', '-v7.3');

    targets = readtable( ...
        targetCSV, ...
        'FileType', 'text', ...
        'VariableNamingRule', 'preserve');

    validate_target_table(targets, targetCSV);

    afdbIDs = strtrim(string(targets.("AFDB ID")));
    afdbIDs = afdbIDs(afdbIDs ~= "");

    edgeKeepFraction = 1;
    edgeRoot = fullfile(outRoot, format_edgekeep_folder(edgeKeepFraction));
    cutoffRoot = fullfile(edgeRoot, format_cutoff_folder(trueCutoff));
    initVariant = char("init_AF_rank_jitter" + string(format_jitter(afRankJitter)));

    manifestRecords = init_manifest_records();
    manifestPath = fullfile(outRoot, 'manifest.csv');

    totalRuns = numel(afdbIDs) * numel(KValues) * numel(noiseSeeds);
    runIndex = 0;

    for pp = 1:numel(afdbIDs)
        afdbID = afdbIDs(pp);
        proteinRoot = fullfile(cutoffRoot, safe_name(afdbID));

        for kk = 1:numel(KValues)
            K = KValues(kk);
            nominalEta = nominalNoiseLevels(kk);
            kRoot = fullfile(proteinRoot, format_K_folder(K));

            for ss = 1:numel(noiseSeeds)
                noiseSeed = noiseSeeds(ss);
                noiseRoot = fullfile(kRoot, format_noise_seed_folder(noiseSeed));
                sharedDir = fullfile(noiseRoot, 'shared');
                runDir = fullfile(noiseRoot, 'runs', initVariant);

                ensure_dir(sharedDir);
                ensure_dir(runDir);

                runIndex = runIndex + 1;

                fprintf('\n============================================================\n');
                fprintf('Run %d/%d | AFDB=%s | eta=%.4g | K=%.8g | seed=%d\n', ...
                    runIndex, totalRuns, char(afdbID), nominalEta, K, noiseSeed);
                fprintf('============================================================\n');

                [problem, sharedStatus, sharedMsg] = build_shared_problem_safely( ...
                    targets, ...
                    targetCSV, ...
                    mappingCSV, ...
                    afdbID, ...
                    trueCutoff, ...
                    nominalEta, ...
                    K, ...
                    noiseSeed, ...
                    afRankSeed, ...
                    afRankJitter, ...
                    opts, ...
                    sharedDir);

                if sharedStatus ~= "ok"
                    write_text(fullfile(runDir, 'status.txt'), "failed_shared_build");
                    write_text(fullfile(runDir, 'log.txt'), sharedMsg);

                    output = struct();
                    output.status = "failed_shared_build";
                    output.error_message = sharedMsg;
                    output.edge_keep_fraction = edgeKeepFraction;
                    output.edge_keep_folder = "edgekeep_100";
                    output.K = K;
                    output.noise_seed = noiseSeed;
                    output.noise_seed_folder = string(format_noise_seed_folder(noiseSeed));
                    output.init_info = fixed_init_info(initVariant, afRankSeed, afRankJitter);
                    output.numit = NaN;
                    output.ReconError = NaN;
                    save(fullfile(runDir, 'solver_output.mat'), 'output', '-v7.3');

                    rec = make_manifest_record( ...
                        runIndex, totalRuns, edgeKeepFraction, noiseSeed, ...
                        afdbID, "", "", trueCutoff, K, initVariant, ...
                        "failed_shared_build", sharedMsg, ...
                        NaN, NaN, NaN, NaN, NaN, runDir);

                    manifestRecords = append_record(manifestRecords, rec);
                    write_manifest(manifestRecords, manifestPath);
                    continue;
                end

                save_shared_problem(problem, sharedDir);

                [status, msg, elapsedSec] = run_one_af_rank_noisy( ...
                    problem, initVariant, afRankSeed, afRankJitter, ...
                    opts, lsopts, runDir);

                rec = make_manifest_record( ...
                    runIndex, totalRuns, edgeKeepFraction, noiseSeed, ...
                    problem.afdb_id, problem.pdb_id, problem.chain_id, ...
                    trueCutoff, K, initVariant, status, msg, ...
                    problem.n_nodes, problem.n_edges, problem.coverage, ...
                    problem.connectivity.n_components, elapsedSec, runDir);

                manifestRecords = append_record(manifestRecords, rec);
                write_manifest(manifestRecords, manifestPath);
            end
        end
    end

    manifest = struct2table(manifestRecords);
    write_manifest(manifestRecords, manifestPath);
end

% -------------------------------------------------------------------------
function validate_inputs( ...
    targetCSV, mappingCSV, trueCutoff, nominalNoiseLevels, KValues, ...
    noiseSeeds, afRankSeed, afRankJitter, opts)

    if ~isfile(targetCSV)
        error('Target CSV not found: %s', targetCSV);
    end

    if ~isfile(mappingCSV)
        error(['Mapping cache not found: %s\n', ...
               'Run scripts/build_mapping_cache.m first.'], mappingCSV);
    end

    if ~isscalar(trueCutoff) || ~isnumeric(trueCutoff) || ...
            ~isfinite(trueCutoff) || trueCutoff <= 0
        error('trueCutoff must be a positive finite numeric scalar.');
    end

    if ~isnumeric(nominalNoiseLevels) || isempty(nominalNoiseLevels) || ...
            any(~isfinite(nominalNoiseLevels)) || any(nominalNoiseLevels <= 0)
        error('nominalNoiseLevels must contain positive finite numeric values.');
    end

    if ~isnumeric(KValues) || isempty(KValues) || ...
            any(~isfinite(KValues)) || any(KValues <= 0)
        error('KValues must contain positive finite numeric values.');
    end

    if numel(nominalNoiseLevels) ~= numel(KValues)
        error('nominalNoiseLevels and KValues must have equal length.');
    end

    if any(abs(KValues(:) - nominalNoiseLevels(:) / 3) > 1e-12)
        error('Each K value must equal its nominal noise level divided by 3.');
    end

    if ~isnumeric(noiseSeeds) || isempty(noiseSeeds) || any(~isfinite(noiseSeeds))
        error('noiseSeeds must contain finite numeric seeds.');
    end

    if ~isscalar(afRankSeed) || ~isnumeric(afRankSeed) || ~isfinite(afRankSeed)
        error('afRankSeed must be a finite numeric scalar.');
    end

    if ~isscalar(afRankJitter) || ~isnumeric(afRankJitter) || ...
            ~isfinite(afRankJitter) || afRankJitter < 0
        error('afRankJitter must be a nonnegative finite numeric scalar.');
    end

    if ~isfield(opts, 'rank') || opts.rank < 3
        error('opts.rank must be at least 3.');
    end
end

% -------------------------------------------------------------------------
function [problem, status, msg] = build_shared_problem_safely( ...
    targets, targetCSV, mappingCSV, afdbID, trueCutoff, nominalEta, K, ...
    noiseSeed, afRankSeed, afRankJitter, opts, sharedDir)

    logPath = fullfile(sharedDir, 'log.txt');
    statusPath = fullfile(sharedDir, 'status.txt');

    if isfile(logPath)
        delete(logPath);
    end

    captured = "";

    try
        captured = evalc(['problem = build_shared_problem(', ...
            'targets, targetCSV, mappingCSV, afdbID, trueCutoff, nominalEta, K, ', ...
            'noiseSeed, afRankSeed, afRankJitter, opts);']);

        status = "ok";
        msg = "";
        write_text(statusPath, status);
        write_text(logPath, captured);

    catch ME
        problem = [];
        status = "failed_shared_build";
        msg = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
        write_text(statusPath, status + newline + string(ME.message));
        write_text(logPath, string(captured) + newline + msg);
    end
end

% -------------------------------------------------------------------------
function problem = build_shared_problem( ...
    targets, targetCSV, mappingCSV, afdbID, trueCutoff, nominalEta, K, ...
    noiseSeed, afRankSeed, afRankJitter, opts)

    row = targets(strcmpi(strtrim(string(targets.("AFDB ID"))), string(afdbID)), :);

    if height(row) ~= 1
        error('Expected exactly one row for AFDB ID %s. Found %d.', ...
            afdbID, height(row));
    end

    pdbID = scalar_text(row.("PDB ID"), 1);
    chainID = scalar_text(row.("Chain ID"), 1);

    truePdbPath = fullfile('data', 'PDB', string(pdbID) + ".pdb");
    afPdbPath   = fullfile('data', 'AFDB', string(afdbID) + ".pdb");

    if ~isfile(truePdbPath)
        error('True PDB file not found: %s', truePdbPath);
    end

    if ~isfile(afPdbPath)
        error('AlphaFold PDB file not found: %s', afPdbPath);
    end

    % Official PDB/SIFTS-to-UniProt mapping, followed by the same full matched
    % subset construction used by run_init_cutoff_sweep.
    outAlign = align_true_vs_af_by_uniprot_mapping( ...
        targetCSV, ...
        afdbID, ...
        'MappingCSV', mappingCSV, ...
        'BuildCacheIfMissing', false, ...
        'Verbose', false);

    sub = subset_aligned_pairs(outAlign, []);

    % Residue-level AF coordinates are retained only to construct the same
    % CA-evaluation metadata saved by run_init_cutoff_sweep.
    [~, coordsAfResidue] = load_ca_plddt(afPdbPath, 'A');
    coordsAfMatchedCA = coordsAfResidue(sub.idxA, :);

    allpairs = build_aligned_all_atom_pairs( ...
        sub, ...
        truePdbPath, ...
        char(chainID), ...
        afPdbPath, ...
        'A', ...
        'AtomSelection', "safe", ...
        'Verbose', false);

    [DistSq, Weight, allatomInfo] = build_allatom_dist_weight_matrix( ...
        allpairs.coords_true_all, ...
        trueCutoff, ...
        'AtomMeta', allpairs.atom_meta, ...
        'MinSeqSep', 0, ...
        'Verbose', false);

    % The complete cutoff graph is used exactly as built. There is no edge
    % subsampling and there are no added true or AlphaFold-derived edges.
    Weight = double(Weight ~= 0);
    Weight = double((Weight + Weight') > 0);
    Weight(1:size(Weight, 1)+1:end) = 0;

    [DistSq, chemNoiseInfo] = apply_chemistry_aware_noise( ...
        DistSq, ...
        Weight, ...
        K, ...
        noiseSeed, ...
        "all_atom", ...
        allpairs.atom_meta, ...
        true, ...
        'Verbose', false);

    caRows = allpairs.ca_atom_rows(:);
    backboneRows = allpairs.backbone_rows(:);

    if isempty(caRows)
        error('No CA evaluation rows were found.');
    end

    connectivity = graph_connectivity_info(Weight);

    Pinit3 = allpairs.coords_af_all - mean(allpairs.coords_af_all, 1);
    pointInitial = af_embed_highdim( ...
        Pinit3, opts.rank, afRankSeed, afRankJitter);

    [subCaEval, coordsAfCaEval] = build_sub_ca_eval( ...
        sub, allpairs, caRows, coordsAfMatchedCA);

    problem = struct();
    problem.afdb_id = string(afdbID);
    problem.pdb_id = string(pdbID);
    problem.chain_id = string(chainID);

    problem.edge_keep_fraction = 1;
    problem.edge_keep_folder = "edgekeep_100";
    problem.true_cutoff = trueCutoff;
    problem.nominal_eta = nominalEta;
    problem.K = K;
    problem.noise_seed = noiseSeed;
    problem.noise_seed_folder = string(format_noise_seed_folder(noiseSeed));
    problem.edg_node_level = "all_atom";

    problem.truePdbPath = string(truePdbPath);
    problem.afPdbPath = string(afPdbPath);

    problem.sub = sub;
    problem.sub_ca_eval = subCaEval;
    problem.coords_af_ca_eval = coordsAfCaEval;
    problem.allpairs = allpairs;
    problem.allatom_info = allatomInfo;

    problem.DistSq = DistSq;
    problem.Weight = Weight;
    problem.coords_true_nodes = allpairs.coords_true_all;
    problem.coords_af_nodes = allpairs.coords_af_all;
    problem.Pinit3 = Pinit3;
    problem.pointInitial = pointInitial;
    problem.ca_eval_rows = caRows;
    problem.backbone_eval_rows = backboneRows;

    problem.af_rank_seed = afRankSeed;
    problem.af_rank_jitter = afRankJitter;
    problem.chem_noise_info = chemNoiseInfo;
    problem.connectivity = connectivity;

    problem.n_nodes = size(Weight, 1);
    problem.n_edges = nnz(Weight) / 2;
    problem.coverage = nnz(Weight) / numel(Weight);
end

% -------------------------------------------------------------------------
function save_shared_problem(problem, sharedDir)
    ensure_dir(sharedDir);

    true_cloud_all = problem.coords_true_nodes;
    af_cloud_all = problem.coords_af_nodes;
    true_cloud_CA = problem.coords_true_nodes(problem.ca_eval_rows, :);
    af_cloud_CA = problem.coords_af_nodes(problem.ca_eval_rows, :);

    save(fullfile(sharedDir, 'clouds.mat'), ...
        'true_cloud_all', ...
        'af_cloud_all', ...
        'true_cloud_CA', ...
        'af_cloud_CA', ...
        '-v7.3');

    % Preserve the run_init_cutoff_sweep initials.mat convention, while
    % storing only the one initialization that this stripped runner uses.
    initials = struct();
    initials.(afrank_fieldname(problem.af_rank_jitter)) = problem.pointInitial;

    save(fullfile(sharedDir, 'initials.mat'), 'initials', '-v7.3');

    W_used = sparse(problem.Weight);
    Dsq_used = sparse_weighted_dist_sq(problem.DistSq, problem.Weight);

    constraint_info = struct();
    constraint_info.edge_keep_fraction = problem.edge_keep_fraction;
    constraint_info.edge_keep_folder = problem.edge_keep_folder;
    constraint_info.true_cutoff = problem.true_cutoff;
    constraint_info.nominal_eta = problem.nominal_eta;
    constraint_info.K = problem.K;
    constraint_info.noise_seed = problem.noise_seed;
    constraint_info.noise_seed_folder = problem.noise_seed_folder;
    constraint_info.n_nodes = problem.n_nodes;
    constraint_info.n_edges = problem.n_edges;
    constraint_info.coverage = problem.coverage;
    constraint_info.n_components = problem.connectivity.n_components;
    constraint_info.component_sizes = problem.connectivity.component_sizes;
    constraint_info.chem_noise_info = problem.chem_noise_info;

    save(fullfile(sharedDir, 'constraints.mat'), ...
        'W_used', ...
        'Dsq_used', ...
        'constraint_info', ...
        '-v7.3');

    ca_rows = problem.ca_eval_rows;
    backbone_rows = problem.backbone_eval_rows;
    atom_meta = problem.allpairs.atom_meta;
    sub_ca_eval = problem.sub_ca_eval;
    coords_af_ca_eval = problem.coords_af_ca_eval;

    metadata = struct();
    metadata.afdb_id = problem.afdb_id;
    metadata.pdb_id = problem.pdb_id;
    metadata.chain_id = problem.chain_id;
    metadata.edge_keep_fraction = problem.edge_keep_fraction;
    metadata.edge_keep_folder = problem.edge_keep_folder;
    metadata.true_cutoff = problem.true_cutoff;
    metadata.nominal_eta = problem.nominal_eta;
    metadata.K = problem.K;
    metadata.noise_seed = problem.noise_seed;
    metadata.noise_seed_folder = problem.noise_seed_folder;
    metadata.edg_node_level = problem.edg_node_level;
    metadata.truePdbPath = problem.truePdbPath;
    metadata.afPdbPath = problem.afPdbPath;
    metadata.af_rank_seed = problem.af_rank_seed;
    metadata.af_rank_jitter = problem.af_rank_jitter;
    metadata.connectivity = problem.connectivity;
    metadata.allatom_info = problem.allatom_info;
    metadata.chem_noise_info = problem.chem_noise_info;

    save(fullfile(sharedDir, 'metadata.mat'), ...
        'ca_rows', ...
        'backbone_rows', ...
        'atom_meta', ...
        'sub_ca_eval', ...
        'coords_af_ca_eval', ...
        'metadata', ...
        '-v7.3');
end

% -------------------------------------------------------------------------
function [status, msg, elapsedSec] = run_one_af_rank_noisy( ...
    problem, initVariant, afRankSeed, afRankJitter, opts, lsopts, runDir)

    ensure_dir(runDir);

    solvedPath = fullfile(runDir, 'solved_cloud.mat');
    outputPath = fullfile(runDir, 'solver_output.mat');
    statusPath = fullfile(runDir, 'status.txt');
    logPath = fullfile(runDir, 'log.txt');

    if isfile(logPath)
        delete(logPath);
    end

    tStart = tic;
    captured = "";

    try
        pointInitial = problem.pointInitial;

        oldFigVisible = get(0, 'DefaultFigureVisible');
        set(0, 'DefaultFigureVisible', 'off');
        cleanupFig = onCleanup(@() set(0, 'DefaultFigureVisible', oldFigVisible)); 

        captured = evalc(['[GCor, ~, output] = alternating_completion_noisy(', ...
            'problem.DistSq, problem.Weight, pointInitial, opts, lsopts);']);

        close all force;

        % No handedness/chirality post-processing is present in this runner.
        solved_cloud_all = GCor;
        solved_cloud_CA = GCor(problem.ca_eval_rows, :);

        save(solvedPath, ...
            'solved_cloud_all', ...
            'solved_cloud_CA', ...
            '-v7.3');

        output.status = "ok";
        output.edge_keep_fraction = problem.edge_keep_fraction;
        output.edge_keep_folder = problem.edge_keep_folder;
        output.nominal_eta = problem.nominal_eta;
        output.K = problem.K;
        output.noise_seed = problem.noise_seed;
        output.noise_seed_folder = problem.noise_seed_folder;
        output.init_info = fixed_init_info(initVariant, afRankSeed, afRankJitter);

        save(outputPath, 'output', '-v7.3');

        status = "ok";
        msg = "";
        elapsedSec = toc(tStart);
        write_text(statusPath, status);
        write_text(logPath, captured);

    catch ME
        status = "failed";
        msg = string(ME.message);
        elapsedSec = toc(tStart);
        report = string(getReport(ME, 'extended', 'hyperlinks', 'off'));

        output = struct();
        output.status = status;
        output.edge_keep_fraction = problem.edge_keep_fraction;
        output.edge_keep_folder = problem.edge_keep_folder;
        output.nominal_eta = problem.nominal_eta;
        output.K = problem.K;
        output.noise_seed = problem.noise_seed;
        output.noise_seed_folder = problem.noise_seed_folder;
        output.init_info = fixed_init_info(initVariant, afRankSeed, afRankJitter);
        output.error_message = msg;
        output.error_report = report;
        output.numit = NaN;
        output.ReconError = NaN;

        save(outputPath, 'output', '-v7.3');
        write_text(statusPath, status + newline + msg);
        write_text(logPath, captured + newline + report);
    end
end

% -------------------------------------------------------------------------
function initInfo = fixed_init_info(initVariant, afRankSeed, afRankJitter)
    initInfo = struct();
    initInfo.variant = string(initVariant);
    initInfo.method = "AF_rank";
    initInfo.random_type = "";
    initInfo.seed = afRankSeed;
    initInfo.jitter = afRankJitter;
end

% -------------------------------------------------------------------------
function [subCaEval, coordsAfCaEval] = build_sub_ca_eval( ...
    sub, allpairs, caEvalRows, coordsAfMatchedCA)

    subCaEval = sub;

    caResiduePairIndex = double( ...
        allpairs.atom_meta.residuePairIndex(caEvalRows));

    if any(~isfinite(caResiduePairIndex)) || ...
            any(caResiduePairIndex < 1) || ...
            any(caResiduePairIndex > numel(sub.idxT))
        error('Invalid residuePairIndex values for CA evaluation rows.');
    end

    caResiduePairIndex = caResiduePairIndex(:);
    subFields = fieldnames(subCaEval);
    nSubRows = numel(sub.idxT);

    for ff = 1:numel(subFields)
        fieldName = subFields{ff};
        value = subCaEval.(fieldName);

        if isnumeric(value) || islogical(value) || isstring(value) || iscell(value)
            if ismatrix(value) && size(value, 1) == nSubRows
                subCaEval.(fieldName) = value(caResiduePairIndex, :);
            elseif isvector(value) && numel(value) == nSubRows
                subCaEval.(fieldName) = value(caResiduePairIndex);
            end
        end
    end

    staleBreakFields = { ...
        'breakAfter_true', 'breakAfter_af', 'breakAfter_any', ...
        'breakAfter_true_matched', 'breakAfter_af_matched', ...
        'breakAfter_any_matched'};

    for bb = 1:numel(staleBreakFields)
        if isfield(subCaEval, staleBreakFields{bb})
            subCaEval = rmfield(subCaEval, staleBreakFields{bb});
        end
    end

    coordsAfCaEval = coordsAfMatchedCA(caResiduePairIndex, :);
end

% -------------------------------------------------------------------------
function info = graph_connectivity_info(Weight)
    W = double(Weight ~= 0);
    W = double((W + W') > 0);
    W(1:size(W, 1)+1:end) = 0;

    G = graph(W);
    bins = conncomp(G);
    componentSizes = accumarray(bins(:), 1);

    info = struct();
    info.n_components = numel(componentSizes);
    info.component_sizes = componentSizes(:);
    info.largest_component_size = max(componentSizes);
    info.is_connected = info.n_components == 1;
end

% -------------------------------------------------------------------------
function DsqUsed = sparse_weighted_dist_sq(DistSq, Weight)
    mask = Weight ~= 0;
    [I, J] = find(mask);
    values = DistSq(sub2ind(size(DistSq), I, J));
    DsqUsed = sparse(I, J, values, size(DistSq, 1), size(DistSq, 2));
end

% -------------------------------------------------------------------------
function validate_target_table(targets, targetCSV)
    required = {'AFDB ID', 'PDB ID', 'Chain ID'};

    for kk = 1:numel(required)
        if ~any(strcmp(required{kk}, targets.Properties.VariableNames))
            error('Target CSV %s is missing required column: %s', ...
                targetCSV, required{kk});
        end
    end
end

% -------------------------------------------------------------------------
function records = init_manifest_records()
    records = struct( ...
        'run_index', {}, ...
        'total_runs', {}, ...
        'edge_keep_fraction', {}, ...
        'edge_keep_folder', {}, ...
        'cutoff', {}, ...
        'K', {}, ...
        'noise_seed', {}, ...
        'noise_seed_folder', {}, ...
        'afdb_id', {}, ...
        'pdb_id', {}, ...
        'chain_id', {}, ...
        'edg_node_level', {}, ...
        'init_variant', {}, ...
        'init_method', {}, ...
        'jitter', {}, ...
        'status', {}, ...
        'message', {}, ...
        'n_nodes', {}, ...
        'n_edges', {}, ...
        'coverage', {}, ...
        'n_components', {}, ...
        'elapsed_sec', {}, ...
        'run_dir', {} ...
    );
end

% -------------------------------------------------------------------------
function rec = make_manifest_record( ...
    runIndex, totalRuns, edgeKeepFraction, noiseSeed, ...
    afdbID, pdbID, chainID, cutoff, K, initVariant, status, msg, ...
    nNodes, nEdges, coverage, nComponents, elapsedSec, runDir)

    rec = struct();
    rec.run_index = runIndex;
    rec.total_runs = totalRuns;
    rec.edge_keep_fraction = edgeKeepFraction;
    rec.edge_keep_folder = string(format_edgekeep_folder(edgeKeepFraction));
    rec.cutoff = cutoff;
    rec.K = K;
    rec.noise_seed = noiseSeed;
    rec.noise_seed_folder = string(format_noise_seed_folder(noiseSeed));
    rec.afdb_id = string(afdbID);
    rec.pdb_id = string(pdbID);
    rec.chain_id = string(chainID);
    rec.edg_node_level = "all_atom";
    rec.init_variant = string(initVariant);
    rec.init_method = "AF_rank";
    rec.jitter = 1e-3;
    rec.status = string(status);
    rec.message = truncate_msg(string(msg), 500);
    rec.n_nodes = nNodes;
    rec.n_edges = nEdges;
    rec.coverage = coverage;
    rec.n_components = nComponents;
    rec.elapsed_sec = elapsedSec;
    rec.run_dir = string(runDir);
end

% -------------------------------------------------------------------------
function records = append_record(records, rec)
    if isempty(records)
        records = rec;
    else
        records(end+1) = rec; 
    end
end

% -------------------------------------------------------------------------
function write_manifest(records, manifestPath)
    if isempty(records)
        return;
    end

    T = struct2table(records);
    writetable(T, manifestPath);
end

% -------------------------------------------------------------------------
function repoRoot = find_repo_root()
    thisFile = mfilename('fullpath');

    if isempty(thisFile)
        repoRoot = pwd;
        return;
    end

    here = fileparts(thisFile);

    if exist(fullfile(here, 'src'), 'dir')
        repoRoot = here;
        return;
    end

    oneUp = fileparts(here);

    if exist(fullfile(oneUp, 'src'), 'dir')
        repoRoot = oneUp;
        return;
    end

    repoRoot = pwd;
end

% -------------------------------------------------------------------------
function textValue = scalar_text(column, index)
    value = column(index);

    if iscell(value)
        textValue = string(value{1});
    else
        textValue = string(value);
    end

    textValue = strtrim(textValue);
end

% -------------------------------------------------------------------------
function folder = format_edgekeep_folder(keepFraction)
    percent = 100 * keepFraction;

    if abs(percent - round(percent)) < 1e-10
        folder = sprintf('edgekeep_%03d', round(percent));
    else
        textValue = sprintf('%.3g', percent);
        textValue = strrep(textValue, '.', 'p');
        textValue = strrep(textValue, '-', 'm');
        folder = ['edgekeep_', textValue];
    end
end

% -------------------------------------------------------------------------
function folder = format_cutoff_folder(cutoff)
    if abs(cutoff - round(cutoff)) < 1e-12
        folder = sprintf('cutoff_%02d', round(cutoff));
    else
        folder = sprintf('cutoff_%s', ...
            strrep(sprintf('%.3g', cutoff), '.', 'p'));
    end
end

% -------------------------------------------------------------------------
function folder = format_K_folder(K)
    if K == 0
        folder = 'K_0';
    else
        textValue = sprintf('%.3g', K);
        textValue = strrep(textValue, '.', 'p');
        textValue = strrep(textValue, '-', 'm');
        folder = ['K_', textValue];
    end
end

% -------------------------------------------------------------------------
function folder = format_noise_seed_folder(seed)
    if abs(seed - round(seed)) < 1e-12
        folder = sprintf('noise_seed_%03d', round(seed));
    else
        textValue = sprintf('%.6g', seed);
        textValue = strrep(textValue, '.', 'p');
        textValue = strrep(textValue, '-', 'm');
        folder = ['noise_seed_', textValue];
    end
end

% -------------------------------------------------------------------------
function fieldName = afrank_fieldname(jitter)
    tag = format_jitter(jitter);
    tag = strrep(tag, '-', 'm');
    tag = strrep(tag, '+', 'p');
    tag = strrep(tag, '.', 'p');
    fieldName = char("afrank_jitter" + string(tag));
end

% -------------------------------------------------------------------------
function textValue = format_jitter(value)
    if value == 0
        textValue = '0';
    else
        textValue = strrep(sprintf('%.0e', value), 'e-0', 'e-');
        textValue = strrep(textValue, 'e+0', 'e');
        textValue = strrep(textValue, 'e+', 'e');
    end
end

% -------------------------------------------------------------------------
function name = safe_name(value)
    name = regexprep(char(string(value)), '[^\w.-]', '_');
end

% -------------------------------------------------------------------------
function output = truncate_msg(message, maxLength)
    message = string(message);

    if strlength(message) > maxLength
        output = extractBefore(message, maxLength) + "...";
    else
        output = message;
    end
end

% -------------------------------------------------------------------------
function ensure_dir(directory)
    if ~exist(directory, 'dir')
        mkdir(directory);
    end
end

% -------------------------------------------------------------------------
function write_text(path, textValue)
    fileID = fopen(path, 'w');

    if fileID == -1
        error('Could not write text file: %s', path);
    end

    cleanupObject = onCleanup(@() fclose(fileID)); 

    payload = char(string(textValue));
    fwrite(fileID, payload, 'char');

    if isempty(payload) || payload(end) ~= char(10)
        fwrite(fileID, char(10), 'char');
    end
end
