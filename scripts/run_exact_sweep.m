% RUN_EXACT_SWEEP
%
% Run with:
%   run('scripts/run_exact_sweep.m')
%
% Output layout:
%   outRoot/cutoff_05/K_0/<AFDB_ID>/shared/
%   outRoot/cutoff_05/K_0/<AFDB_ID>/runs/init_rand_seed47/
%   outRoot/cutoff_05/K_0/<AFDB_ID>/runs/init_randn_seed47/
%   outRoot/cutoff_05/K_0/<AFDB_ID>/runs/init_floyd/
%   outRoot/cutoff_05/K_0/<AFDB_ID>/runs/init_AF_3D/

close all force;

% =========================================================================
% FIXED PAPER SETTINGS
% =========================================================================

targetCSV  = fullfile('configs', 'exact_targets.csv');
outRoot    = fullfile('out', 'init_cutoff_sweep');
trueCutoff = 5;
randomSeed = 47;

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

manifest = run_exact_sweep_impl(targetCSV, outRoot, trueCutoff, randomSeed, opts, lsopts);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function manifest = run_exact_sweep_impl(targetCSV, outRoot, trueCutoff, randomSeed, opts, lsopts)

    repoRoot = find_repo_root();
    cd(repoRoot);
    addpath(genpath(fullfile(repoRoot, 'src')));

    targetCSV = char(string(targetCSV));
    outRoot   = char(string(outRoot));
    mappingCSV = fullfile('configs', 'pdb_uniprot_residue_map.csv');

    if ~isfile(targetCSV)
        error('Target CSV not found: %s', targetCSV);
    end

    if ~isfile(mappingCSV)
        error(['Mapping cache not found: %s\n', ...
               'Run scripts/build_mapping_cache.m first.'], mappingCSV);
    end

    if ~isscalar(trueCutoff) || ~isnumeric(trueCutoff) || ~isfinite(trueCutoff) || trueCutoff <= 0
        error('trueCutoff must be a positive numeric scalar.');
    end

    if ~isscalar(randomSeed) || ~isnumeric(randomSeed) || ~isfinite(randomSeed)
        error('randomSeed must be a finite numeric scalar.');
    end

    if ~isfield(opts, 'rank') || opts.rank < 3
        error('opts.rank must be at least 3.');
    end

    ensure_dir(outRoot);

    run_config = struct();
    run_config.targetCSV = targetCSV;
    run_config.outRoot = outRoot;
    run_config.trueCutoff = trueCutoff;
    run_config.randomSeed = randomSeed;
    run_config.opts = opts;
    run_config.lsopts = lsopts;
    save(fullfile(outRoot, 'sweep_config.mat'), 'run_config', '-v7.3');

    targets = readtable(targetCSV, 'FileType', 'text', 'VariableNamingRule', 'preserve');
    validate_target_table(targets, targetCSV);

    afdbIDs = strtrim(string(targets.("AFDB ID")));
    afdbIDs = afdbIDs(afdbIDs ~= "");

    initCases = struct( ...
        'variant', {'init_rand_seed47', 'init_randn_seed47', 'init_floyd', 'init_AF_3D'}, ...
        'kind', {'rand', 'randn', 'floyd', 'af3d'} ...
    );

    manifestRecords = init_manifest_records();
    manifestPath = fullfile(outRoot, 'manifest.csv');

    totalRuns = numel(afdbIDs) * numel(initCases);
    runIndex = 0;

    % Fixed exact-branch hierarchy: edgekeep -> cutoff -> protein -> K.
    cutoffRoot = fullfile(outRoot, 'edgekeep_100', 'cutoff_05');

    for pp = 1:numel(afdbIDs)
        afdbID = afdbIDs(pp);
        proteinRoot = fullfile(cutoffRoot, safe_name(afdbID), 'K_0');
        sharedDir = fullfile(proteinRoot, 'shared');
        ensure_dir(sharedDir);

        [problem, sharedStatus, sharedMsg] = build_shared_problem_safely( ...
            targets, targetCSV, mappingCSV, afdbID, trueCutoff, randomSeed, opts, sharedDir);

        if sharedStatus ~= "ok"
            for ii = 1:numel(initCases)
                runIndex = runIndex + 1;
                initCase = initCases(ii);
                runDir = fullfile(proteinRoot, 'runs', initCase.variant);
                ensure_dir(runDir);

                write_text(fullfile(runDir, 'status.txt'), "failed_shared_build");
                write_text(fullfile(runDir, 'log.txt'), sharedMsg);

                output = struct();
                output.status = "failed_shared_build";
                output.error_message = sharedMsg;
                output.initialization = string(initCase.kind);
                output.K = 0;
                save(fullfile(runDir, 'solver_output.mat'), 'output', '-v7.3');

                rec = make_manifest_record(runIndex, totalRuns, afdbID, "", "", trueCutoff, ...
                    initCase.variant, initCase.kind, "failed_shared_build", sharedMsg, ...
                    NaN, NaN, NaN, NaN, runDir);

                manifestRecords = append_record(manifestRecords, rec);
            end

            write_manifest(manifestRecords, manifestPath);
            continue;
        end

        save_shared_problem(problem, sharedDir);

        for ii = 1:numel(initCases)
            runIndex = runIndex + 1;
            initCase = initCases(ii);
            runDir = fullfile(proteinRoot, 'runs', initCase.variant);
            ensure_dir(runDir);

            [status, msg, elapsedSec] = run_one_initialization(problem, initCase, opts, lsopts, runDir);

            rec = make_manifest_record(runIndex, totalRuns, ...
                problem.afdb_id, problem.pdb_id, problem.chain_id, trueCutoff, ...
                initCase.variant, initCase.kind, status, msg, ...
                problem.n_nodes, problem.n_edges, problem.coverage, elapsedSec, runDir);

            manifestRecords = append_record(manifestRecords, rec);
            write_manifest(manifestRecords, manifestPath);
        end
    end

    manifest = struct2table(manifestRecords);
    write_manifest(manifestRecords, manifestPath);
end

% -------------------------------------------------------------------------
function [problem, status, msg] = build_shared_problem_safely( ...
    targets, targetCSV, mappingCSV, afdbID, trueCutoff, randomSeed, opts, sharedDir)

    logPath = fullfile(sharedDir, 'log.txt');
    statusPath = fullfile(sharedDir, 'status.txt');

    if isfile(logPath)
        delete(logPath);
    end

    captured = "";

    try
        captured = evalc(['problem = build_shared_problem(', ...
            'targets, targetCSV, mappingCSV, afdbID, trueCutoff, randomSeed, opts);']);

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
function problem = build_shared_problem(targets, targetCSV, mappingCSV, afdbID, trueCutoff, randomSeed, opts)

    row = targets(strcmpi(strtrim(string(targets.("AFDB ID"))), string(afdbID)), :);

    if height(row) ~= 1
        error('Expected exactly one row for AFDB ID %s. Found %d.', afdbID, height(row));
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

    alignment = align_true_vs_af_by_uniprot_mapping( ...
        targetCSV, afdbID, ...
        'MappingCSV', mappingCSV, ...
        'BuildCacheIfMissing', false, ...
        'Verbose', false);

    % Use the complete mapped residue set, matching DomainCut=[] in the
    % original full sweep.
    matchedPairs = subset_aligned_pairs(alignment, []);

    allpairs = build_aligned_all_atom_pairs( ...
        matchedPairs, ...
        truePdbPath, ...
        char(chainID), ...
        afPdbPath, ...
        'A', ...
        'AtomSelection', "safe", ...
        'Verbose', false);

    [DistSq, Weight, allatomInfo] = build_allatom_dist_weight_matrix( ...
        allpairs.coords_true_all, trueCutoff, ...
        'AtomMeta', allpairs.atom_meta, ...
        'MinSeqSep', 0, ...
        'Verbose', false);

    Weight = double(Weight ~= 0);
    Weight = double((Weight + Weight') > 0);
    Weight(1:size(Weight, 1)+1:end) = 0;

    caRows = allpairs.ca_atom_rows(:);
    backboneRows = allpairs.backbone_rows(:);

    if isempty(caRows)
        error('No CA rows found for output.');
    end

    connectivity = graph_connectivity_info(Weight);
    Pinit3 = allpairs.coords_af_all - mean(allpairs.coords_af_all, 1);

    rng(randomSeed, 'twister');
    randPoint = rand(size(DistSq, 1), opts.rank);

    rng(randomSeed, 'twister');
    randnPoint = randn(size(DistSq, 1), opts.rank);

    floydPoint = [];
    floydApplicable = false;
    floydReason = "";

    if connectivity.n_components > 1
        floydReason = "not_applicable_disconnected";
    else
        [~, floydPointCandidate, floydRows] = floyd_warshall_init(DistSq, Weight, 3);

        if numel(floydRows) ~= size(DistSq, 1)
            floydReason = "not_applicable_internal_inf_trim";
        else
            floydPoint = floydPointCandidate;
            floydApplicable = true;
        end
    end

    initials = struct();
    initials.rand = randPoint;
    initials.randn = randnPoint;
    initials.floyd = floydPoint;
    initials.floyd_applicable = floydApplicable;
    initials.floyd_reason = floydReason;
    initials.af3d = Pinit3;

    problem = struct();
    problem.afdb_id = string(afdbID);
    problem.pdb_id = string(pdbID);
    problem.chain_id = string(chainID);
    problem.true_cutoff = trueCutoff;
    problem.K = 0;
    problem.truePdbPath = string(truePdbPath);
    problem.afPdbPath = string(afPdbPath);
    problem.alignment = alignment;
    problem.allpairs = allpairs;
    problem.allatom_info = allatomInfo;
    problem.DistSq = DistSq;
    problem.Weight = Weight;
    problem.coords_true_all = allpairs.coords_true_all;
    problem.coords_af_all = allpairs.coords_af_all;
    problem.coords_true_CA = allpairs.coords_true_all(caRows, :);
    problem.coords_af_CA = allpairs.coords_af_all(caRows, :);
    problem.ca_rows = caRows;
    problem.backbone_rows = backboneRows;
    problem.atom_meta = allpairs.atom_meta;
    problem.connectivity = connectivity;
    problem.initials = initials;
    problem.n_nodes = size(Weight, 1);
    problem.n_edges = nnz(triu(Weight, 1));
    problem.coverage = nnz(Weight) / numel(Weight);
end

% -------------------------------------------------------------------------
function save_shared_problem(problem, sharedDir)
    ensure_dir(sharedDir);

    true_cloud_all = problem.coords_true_all;
    af_cloud_all = problem.coords_af_all;
    true_cloud_CA = problem.coords_true_CA;
    af_cloud_CA = problem.coords_af_CA;

    save(fullfile(sharedDir, 'clouds.mat'), ...
        'true_cloud_all', 'af_cloud_all', 'true_cloud_CA', 'af_cloud_CA', '-v7.3');

    initials = problem.initials;
    save(fullfile(sharedDir, 'initials.mat'), 'initials', '-v7.3');

    W_used = sparse(problem.Weight);
    Dsq_used = sparse_weighted_dist_sq(problem.DistSq, problem.Weight);

    constraint_info = struct();
    constraint_info.true_cutoff = problem.true_cutoff;
    constraint_info.K = problem.K;
    constraint_info.n_nodes = problem.n_nodes;
    constraint_info.n_edges = problem.n_edges;
    constraint_info.coverage = problem.coverage;
    constraint_info.n_components = problem.connectivity.n_components;
    constraint_info.component_sizes = problem.connectivity.component_sizes;

    save(fullfile(sharedDir, 'constraints.mat'), ...
        'W_used', 'Dsq_used', 'constraint_info', '-v7.3');

    ca_rows = problem.ca_rows;
    backbone_rows = problem.backbone_rows;
    atom_meta = problem.atom_meta;
    alignment = problem.alignment;

    metadata = struct();
    metadata.afdb_id = problem.afdb_id;
    metadata.pdb_id = problem.pdb_id;
    metadata.chain_id = problem.chain_id;
    metadata.true_cutoff = problem.true_cutoff;
    metadata.K = problem.K;
    metadata.truePdbPath = problem.truePdbPath;
    metadata.afPdbPath = problem.afPdbPath;
    metadata.connectivity = problem.connectivity;
    metadata.allatom_info = problem.allatom_info;

    save(fullfile(sharedDir, 'metadata.mat'), ...
        'ca_rows', 'backbone_rows', 'atom_meta', 'alignment', 'metadata', '-v7.3');
end

% -------------------------------------------------------------------------
function [status, msg, elapsedSec] = run_one_initialization(problem, initCase, opts, lsopts, runDir)
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
        [pointInitial, skipStatus] = initial_point(problem, initCase.kind);

        if skipStatus ~= ""
            status = skipStatus;
            msg = skipStatus;

            output = struct();
            output.status = status;
            output.initialization = string(initCase.kind);
            output.reason = skipStatus;
            output.K = problem.K;
            output.numit = NaN;
            output.ReconError = NaN;

            save(outputPath, 'output', '-v7.3');
            write_text(statusPath, status);
            write_text(logPath, "");

            elapsedSec = toc(tStart);
            return;
        end

        oldFigVisible = get(0, 'DefaultFigureVisible');
        set(0, 'DefaultFigureVisible', 'off');
        cleanupFig = onCleanup(@() set(0, 'DefaultFigureVisible', oldFigVisible)); 

        captured = evalc(['[GCor, IPM_Recon, output] = alternating_completion(', ...
            'problem.DistSq, problem.Weight, pointInitial, opts, lsopts);']);

        close all force;

        % Save the raw solver coordinates directly. No post-solver
        % reflection, handedness, chirality, or other transform is present.
        solved_cloud_all = GCor;
        solved_cloud_CA = GCor(problem.ca_rows, :);
        save(solvedPath, 'solved_cloud_all', 'solved_cloud_CA', '-v7.3');

        output.status = "ok";
        output.initialization = string(initCase.kind);
        output.K = problem.K;
        save(outputPath, 'output', 'IPM_Recon', '-v7.3');

        status = "ok";
        msg = "";
        elapsedSec = toc(tStart);

        write_text(statusPath, status);
        write_text(logPath, string(captured));

    catch ME
        status = "failed";
        msg = string(ME.message);
        elapsedSec = toc(tStart);
        report = getReport(ME, 'extended', 'hyperlinks', 'off');

        output = struct();
        output.status = status;
        output.initialization = string(initCase.kind);
        output.K = problem.K;
        output.error_message = msg;
        output.error_report = string(report);
        output.numit = NaN;
        output.ReconError = NaN;

        save(outputPath, 'output', '-v7.3');
        write_text(statusPath, status + newline + msg);
        write_text(logPath, string(captured) + newline + string(report));
    end
end

% -------------------------------------------------------------------------
function [P0, skipStatus] = initial_point(problem, kind)
    skipStatus = "";

    switch lower(string(kind))
        case "rand"
            P0 = problem.initials.rand;
        case "randn"
            P0 = problem.initials.randn;
        case "floyd"
            if ~problem.initials.floyd_applicable
                P0 = [];
                skipStatus = string(problem.initials.floyd_reason);
                return;
            end
            P0 = problem.initials.floyd;
        case "af3d"
            P0 = problem.initials.af3d;
        otherwise
            error('Unknown initialization: %s.', string(kind));
    end
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
function Dsq_used = sparse_weighted_dist_sq(DistSq, Weight)
    mask = Weight ~= 0;
    [I, J] = find(mask);
    vals = DistSq(sub2ind(size(DistSq), I, J));
    Dsq_used = sparse(I, J, vals, size(DistSq, 1), size(DistSq, 2));
end

% -------------------------------------------------------------------------
function validate_target_table(targets, targetCSV)
    required = {'AFDB ID', 'PDB ID', 'Chain ID'};
    for k = 1:numel(required)
        if ~any(strcmp(required{k}, targets.Properties.VariableNames))
            error('Target CSV %s is missing required column: %s', targetCSV, required{k});
        end
    end
end

% -------------------------------------------------------------------------
function records = init_manifest_records()
    records = struct( ...
        'run_index', {}, 'total_runs', {}, 'afdb_id', {}, 'pdb_id', {}, ...
        'chain_id', {}, 'cutoff', {}, 'K', {}, 'init_variant', {}, ...
        'init_method', {}, 'status', {}, 'message', {}, 'n_nodes', {}, ...
        'n_edges', {}, 'coverage', {}, 'elapsed_sec', {}, 'run_dir', {});
end

% -------------------------------------------------------------------------
function rec = make_manifest_record(runIndex, totalRuns, afdbID, pdbID, chainID, cutoff, ...
    initVariant, initMethod, status, msg, nNodes, nEdges, coverage, elapsedSec, runDir)

    rec = struct();
    rec.run_index = runIndex;
    rec.total_runs = totalRuns;
    rec.afdb_id = string(afdbID);
    rec.pdb_id = string(pdbID);
    rec.chain_id = string(chainID);
    rec.cutoff = cutoff;
    rec.K = 0;
    rec.init_variant = string(initVariant);
    rec.init_method = string(initMethod);
    rec.status = string(status);
    rec.message = truncate_msg(string(msg), 500);
    rec.n_nodes = nNodes;
    rec.n_edges = nEdges;
    rec.coverage = coverage;
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
function txt = scalar_text(col, idx)
    v = col(idx);
    if iscell(v)
        txt = string(v{1});
    else
        txt = string(v);
    end
    txt = strtrim(txt);
end

% -------------------------------------------------------------------------
function s = safe_name(x)
    s = regexprep(char(string(x)), '[^\w.-]', '_');
end

% -------------------------------------------------------------------------
function msg = truncate_msg(msg, maxLen)
    msg = string(msg);
    if strlength(msg) > maxLen
        msg = extractBefore(msg, maxLen) + "...";
    end
end

% -------------------------------------------------------------------------
function ensure_dir(d)
    if ~exist(d, 'dir')
        mkdir(d);
    end
end

% -------------------------------------------------------------------------
function write_text(path, txt)
    fid = fopen(path, 'w');
    if fid < 0
        error('Could not open file for writing: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid)); 
    fprintf(fid, '%s', char(string(txt)));
end

