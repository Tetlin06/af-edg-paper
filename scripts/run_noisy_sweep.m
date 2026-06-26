% RUN_NOISY_SWEEP
%
% Run with:
%   run('scripts/run_noisy_sweep.m')
%
% Output layout:
%   outRoot/cutoff_06/noise_0p003333/noise_seed_006/<AFDB_ID>/shared/
%   outRoot/cutoff_06/noise_0p003333/noise_seed_006/<AFDB_ID>/runs/init_AF_rank_jitter0p001/

close all force;

% =========================================================================
% USER SETTINGS
% =========================================================================

targetCSV = fullfile('configs', 'noisy_targets.csv');
outRoot   = fullfile('out', 'noisy_sweep');

trueCutoff = 6;
noise      = 0.01 / 3;
noiseSeed  = 6;
afJitter   = 1e-3;

opts = struct();
opts.r           = 10000;
opts.printenergy = 1;
opts.printerror  = 1;
opts.rank        = 10;
opts.maxit       = 2000;
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
    targetCSV, outRoot, trueCutoff, noise, noiseSeed, afJitter, opts, lsopts);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function manifest = run_noisy_sweep_impl(targetCSV, outRoot, trueCutoff, noise, noiseSeed, afJitter, opts, lsopts)

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

    if ~isscalar(noise) || ~isnumeric(noise) || ~isfinite(noise) || noise <= 0
        error('noise must be a positive numeric scalar.');
    end

    if ~isscalar(noiseSeed) || ~isnumeric(noiseSeed) || ~isfinite(noiseSeed)
        error('noiseSeed must be a finite numeric scalar.');
    end

    if ~isscalar(afJitter) || ~isnumeric(afJitter) || ~isfinite(afJitter) || afJitter < 0
        error('afJitter must be a nonnegative numeric scalar.');
    end

    if ~isfield(opts, 'rank') || opts.rank < 3
        error('opts.rank must be at least 3.');
    end

    ensure_dir(outRoot);

    run_config = struct();
    run_config.targetCSV = targetCSV;
    run_config.outRoot = outRoot;
    run_config.trueCutoff = trueCutoff;
    run_config.noise = noise;
    run_config.noiseSeed = noiseSeed;
    run_config.afJitter = afJitter;
    run_config.opts = opts;
    run_config.lsopts = lsopts;
    save(fullfile(outRoot, 'sweep_config.mat'), 'run_config', '-v7.3');

    targets = readtable(targetCSV, 'FileType', 'text', 'VariableNamingRule', 'preserve');
    validate_target_table(targets, targetCSV);

    afdbIDs = strtrim(string(targets.("AFDB ID")));
    afdbIDs = afdbIDs(afdbIDs ~= "");

    initVariant = "init_AF_rank_jitter" + format_float_token(afJitter);

    manifestRecords = init_manifest_records();
    manifestPath = fullfile(outRoot, 'manifest.csv');

    totalRuns = numel(afdbIDs);
    runIndex = 0;

    baseRoot = fullfile( ...
        outRoot, ...
        format_cutoff_folder(trueCutoff), ...
        format_noise_folder(noise), ...
        format_noise_seed_folder(noiseSeed));

    for pp = 1:numel(afdbIDs)
        runIndex = runIndex + 1;
        afdbID = afdbIDs(pp);
        proteinRoot = fullfile(baseRoot, safe_name(afdbID));
        sharedDir = fullfile(proteinRoot, 'shared');
        runDir = fullfile(proteinRoot, 'runs', char(initVariant));

        ensure_dir(sharedDir);
        ensure_dir(runDir);

        [problem, sharedStatus, sharedMsg] = build_shared_problem_safely( ...
            targets, targetCSV, mappingCSV, afdbID, trueCutoff, noise, noiseSeed, afJitter, opts, sharedDir);

        if sharedStatus ~= "ok"
            write_text(fullfile(runDir, 'status.txt'), "failed_shared_build");
            write_text(fullfile(runDir, 'log.txt'), sharedMsg);

            output = struct();
            output.status = "failed_shared_build";
            output.error_message = sharedMsg;
            output.initialization = "AF_rank";
            output.noise = noise;
            output.noise_seed = noiseSeed;
            output.af_jitter = afJitter;
            save(fullfile(runDir, 'solver_output.mat'), 'output', '-v7.3');

            rec = make_manifest_record(runIndex, totalRuns, afdbID, "", "", trueCutoff, ...
                noise, noiseSeed, afJitter, initVariant, "failed_shared_build", sharedMsg, ...
                NaN, NaN, NaN, NaN, runDir);

            manifestRecords = append_record(manifestRecords, rec);
            write_manifest(manifestRecords, manifestPath);
            continue;
        end

        save_shared_problem(problem, sharedDir);
        [status, msg, elapsedSec] = run_one_af_rank(problem, opts, lsopts, runDir);

        rec = make_manifest_record(runIndex, totalRuns, ...
            problem.afdb_id, problem.pdb_id, problem.chain_id, trueCutoff, ...
            noise, noiseSeed, afJitter, initVariant, status, msg, ...
            problem.n_nodes, problem.n_edges, problem.coverage, elapsedSec, runDir);

        manifestRecords = append_record(manifestRecords, rec);
        write_manifest(manifestRecords, manifestPath);
    end

    manifest = struct2table(manifestRecords);
    write_manifest(manifestRecords, manifestPath);
end

% -------------------------------------------------------------------------
function [problem, status, msg] = build_shared_problem_safely( ...
    targets, targetCSV, mappingCSV, afdbID, trueCutoff, noise, noiseSeed, afJitter, opts, sharedDir)

    logPath = fullfile(sharedDir, 'log.txt');
    statusPath = fullfile(sharedDir, 'status.txt');

    if isfile(logPath)
        delete(logPath);
    end

    captured = "";

    try
        captured = evalc(['problem = build_shared_problem(', ...
            'targets, targetCSV, mappingCSV, afdbID, trueCutoff, noise, noiseSeed, afJitter, opts);']);

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
function problem = build_shared_problem(targets, targetCSV, mappingCSV, afdbID, trueCutoff, noise, noiseSeed, afJitter, opts)

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

    allpairs = build_aligned_all_atom_pairs( ...
        alignment, truePdbPath, char(chainID), afPdbPath, 'A');

    [DistSq, Weight, allatomInfo] = build_allatom_dist_weight_matrix( ...
        allpairs.coords_true_all, trueCutoff, ...
        'AtomMeta', allpairs.atom_meta, ...
        'MinSeqSep', 0, ...
        'Verbose', false);

    Weight = double(Weight ~= 0);
    Weight = double((Weight + Weight') > 0);
    Weight(1:size(Weight, 1)+1:end) = 0;

    [DistSq, noiseInfo] = apply_chemistry_aware_noise( ...
        DistSq, Weight, noise, noiseSeed, "all_atom", allpairs.atom_meta, true, 'Verbose', false);

    caRows = allpairs.ca_atom_rows(:);
    backboneRows = allpairs.backbone_rows(:);

    if isempty(caRows)
        error('No CA rows found for output.');
    end

    Pinit3 = allpairs.coords_af_all - mean(allpairs.coords_af_all, 1);
    [pointInitial, afQ] = af_embed_highdim(Pinit3, opts.rank, 47, afJitter);

    problem = struct();
    problem.afdb_id = string(afdbID);
    problem.pdb_id = string(pdbID);
    problem.chain_id = string(chainID);
    problem.true_cutoff = trueCutoff;
    problem.noise = noise;
    problem.noise_seed = noiseSeed;
    problem.af_jitter = afJitter;
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
    problem.pointInitial = pointInitial;
    problem.afQ = afQ;
    problem.Pinit3 = Pinit3;
    problem.noise_info = noiseInfo;
    problem.connectivity = graph_connectivity_info(Weight);
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

    initial = problem.pointInitial;
    Pinit3 = problem.Pinit3;
    afQ = problem.afQ;
    af_jitter = problem.af_jitter;

    save(fullfile(sharedDir, 'initials.mat'), ...
        'initial', 'Pinit3', 'afQ', 'af_jitter', '-v7.3');

    W_used = sparse(problem.Weight);
    Dsq_used = sparse_weighted_dist_sq(problem.DistSq, problem.Weight);

    constraint_info = struct();
    constraint_info.true_cutoff = problem.true_cutoff;
    constraint_info.noise = problem.noise;
    constraint_info.noise_seed = problem.noise_seed;
    constraint_info.n_nodes = problem.n_nodes;
    constraint_info.n_edges = problem.n_edges;
    constraint_info.coverage = problem.coverage;
    constraint_info.n_components = problem.connectivity.n_components;
    constraint_info.component_sizes = problem.connectivity.component_sizes;
    constraint_info.noise_info = problem.noise_info;

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
    metadata.noise = problem.noise;
    metadata.noise_seed = problem.noise_seed;
    metadata.af_jitter = problem.af_jitter;
    metadata.truePdbPath = problem.truePdbPath;
    metadata.afPdbPath = problem.afPdbPath;
    metadata.connectivity = problem.connectivity;
    metadata.allatom_info = problem.allatom_info;

    save(fullfile(sharedDir, 'metadata.mat'), ...
        'ca_rows', 'backbone_rows', 'atom_meta', 'alignment', 'metadata', '-v7.3');
end

% -------------------------------------------------------------------------
function [status, msg, elapsedSec] = run_one_af_rank(problem, opts, lsopts, runDir)
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

        captured = evalc(['[GCor, IPM_Recon, output] = alternating_completion_noisy(', ...
            'problem.DistSq, problem.Weight, pointInitial, opts, lsopts);']);

        close all force;

        solved_cloud_all = GCor;
        solved_cloud_CA = GCor(problem.ca_rows, :);
        save(solvedPath, 'solved_cloud_all', 'solved_cloud_CA', '-v7.3');

        output.status = "ok";
        output.initialization = "AF_rank";
        output.noise = problem.noise;
        output.noise_seed = problem.noise_seed;
        output.af_jitter = problem.af_jitter;
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
        output.initialization = "AF_rank";
        output.noise = problem.noise;
        output.noise_seed = problem.noise_seed;
        output.af_jitter = problem.af_jitter;
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
        'chain_id', {}, 'cutoff', {}, 'noise', {}, 'noise_seed', {}, ...
        'af_jitter', {}, 'init_variant', {}, 'init_method', {}, 'status', {}, ...
        'message', {}, 'n_nodes', {}, 'n_edges', {}, 'coverage', {}, ...
        'elapsed_sec', {}, 'run_dir', {});
end

% -------------------------------------------------------------------------
function rec = make_manifest_record(runIndex, totalRuns, afdbID, pdbID, chainID, cutoff, ...
    noise, noiseSeed, afJitter, initVariant, status, msg, nNodes, nEdges, coverage, elapsedSec, runDir)

    rec = struct();
    rec.run_index = runIndex;
    rec.total_runs = totalRuns;
    rec.afdb_id = string(afdbID);
    rec.pdb_id = string(pdbID);
    rec.chain_id = string(chainID);
    rec.cutoff = cutoff;
    rec.noise = noise;
    rec.noise_seed = noiseSeed;
    rec.af_jitter = afJitter;
    rec.init_variant = string(initVariant);
    rec.init_method = "AF_rank";
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
function folder = format_cutoff_folder(cutoff)
    if abs(cutoff - round(cutoff)) < 1e-12
        folder = sprintf('cutoff_%02d', round(cutoff));
    else
        folder = sprintf('cutoff_%s', strrep(sprintf('%.3g', cutoff), '.', 'p'));
    end
end

% -------------------------------------------------------------------------
function folder = format_noise_folder(noise)
    folder = ['noise_', char(format_float_token(noise))];
end

% -------------------------------------------------------------------------
function folder = format_noise_seed_folder(seed)
    if abs(seed - round(seed)) < 1e-12
        folder = sprintf('noise_seed_%03d', round(seed));
    else
        folder = ['noise_seed_', char(format_float_token(seed))];
    end
end

% -------------------------------------------------------------------------
function token = format_float_token(x)
    if abs(x - round(x)) < 1e-12
        token = string(round(x));
    else
        token = string(sprintf('%.6g', x));
    end
    token = strrep(token, '.', 'p');
    token = strrep(token, '-', 'm');
    token = strrep(token, '+', '');
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
