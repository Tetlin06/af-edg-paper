function out = align_true_vs_af_by_uniprot_mapping(varargin)
% ALIGN_TRUE_VS_AF_BY_UNIPROT_MAPPING
% -------------------------------------------------------------------------
% Build TRUE PDB <-> AlphaFold residue correspondence using official
% PDB -> UniProt mapping, rather than Needleman-Wunsch sequence alignment.
%
% Drop-in replacement for align_true_vs_af_by_full_sequence in main.m.
%
% MODE 1: targets.csv lookup
%   out = align_true_vs_af_by_uniprot_mapping('sheets/targets.csv','Q5VSL9-2')
%
% MODE 2: explicit paths
%   out = align_true_vs_af_by_uniprot_mapping(truePdbPath,trueChain,afPdbPath,afChain, ...
%                                            'PDBID','7K36','AFDBID','Q5VSL9-2')
%
% Required cache file:
%   sheets/pdb_uniprot_residue_map.csv
%
% You can generate it with:
%   build_pdb_uniprot_mapping_cache('sheets/targets.csv')

    if nargin < 2
        error('Not enough inputs.');
    end

    firstArg = string(varargin{1});

    if endsWith(firstArg, ".csv")
        % MODE 1: lookup from targets.csv
        targetCSV = string(varargin{1});
        afdbID = string(varargin{2});
        extraArgs = varargin(3:end);

        if ~isfile(targetCSV)
            error('targets.csv not found: %s', targetCSV);
        end

        Ttargets = readtable(targetCSV, 'FileType','text', 'VariableNamingRule','preserve');
        afdbCol = string(Ttargets.("AFDB ID"));
        row = Ttargets(strcmpi(strtrim(afdbCol), afdbID), :);

        if height(row) ~= 1
            error('Expected exactly one row in %s for AFDB ID %s. Found %d.', targetCSV, afdbID, height(row));
        end

        pdbID = upper(scalar_text(row.("PDB ID"), 1));
        trueChain = char(scalar_text(row.("Chain ID"), 1));

        truePdbPath = fullfile('data','PDB', pdbID + ".pdb");
        afPdbPath   = fullfile('data','AFDB', afdbID + ".pdb");
        afChain     = 'A';

        varargin = extraArgs;

    else
        % MODE 2: explicit paths
        if nargin < 4
            error('Explicit mode requires truePdbPath, trueChain, afPdbPath, afChain.');
        end

        targetCSV = "";
        truePdbPath = string(varargin{1});
        trueChain   = char(string(varargin{2}));
        afPdbPath   = string(varargin{3});
        afChain     = char(string(varargin{4}));
        afdbID      = "";
        pdbID       = "";

        varargin = varargin(5:end);
    end

    % ---------------- options ----------------
    ip = inputParser;
    ip.addParameter('MappingCSV', fullfile('sheets','pdb_uniprot_residue_map.csv'), @(s)ischar(s)||isstring(s));
    ip.addParameter('BuildCacheIfMissing', true, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('SiftsDir', fullfile('data','SIFTS'), @(s)ischar(s)||isstring(s));
    ip.addParameter('AFDBID', afdbID, @(s)ischar(s)||isstring(s));
    ip.addParameter('PDBID', pdbID, @(s)ischar(s)||isstring(s));
    ip.addParameter('UniProtID', "", @(s)ischar(s)||isstring(s));
    ip.addParameter('AllowReflection', false, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('MakeFigure', false, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('SaveFigureTo', '', @(s)ischar(s)||isstring(s));
    ip.addParameter('Verbose', true, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('MinAlignedPairs', 30, @(x)isnumeric(x)&&isscalar(x));
    ip.addParameter('RequireObserved', true, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('PerfectMatchesOnly', false, @(x)islogical(x)||isnumeric(x));
    ip.parse(varargin{:});
    opt = ip.Results;

    mappingCSV = string(opt.MappingCSV);
    afdbID = string(opt.AFDBID);
    pdbID = upper(string(opt.PDBID));
    uniWanted = string(opt.UniProtID);

    if strlength(uniWanted) == 0 && strlength(afdbID) > 0
        uniWanted = base_uniprot(afdbID);
    end

    if strlength(pdbID) == 0
        [~, stem, ~] = fileparts(truePdbPath);
        pdbID = upper(string(stem));
    end

    if ~isfile(truePdbPath)
        error('True PDB not found: %s', string(truePdbPath));
    end
    if ~isfile(afPdbPath)
        error('AF PDB not found: %s', string(afPdbPath));
    end

    if ~isfile(mappingCSV)
        if strlength(targetCSV) > 0 && logical(opt.BuildCacheIfMissing)
            build_pdb_uniprot_mapping_cache(targetCSV, ...
                'SiftsDir', opt.SiftsDir, ...
                'OutCSV', mappingCSV, ...
                'Verbose', logical(opt.Verbose));
        else
            error('Mapping CSV not found: %s', mappingCSV);
        end
    end

    % ---------------- load CA coordinates through shared loader ----------------
    [seqT_ca, coordsT_ca, residT_ca, resseqT_ca, unknownT_ca, statsT_ca] = ...
        load_ca_sequence_coords(truePdbPath, trueChain, ...
            'AltLocPolicy', 'blank_or_A', ...
            'Verbose', false);

    [seqA_ca, coordsA_ca, residA_ca, resseqA_ca, unknownA_ca, statsA_ca] = ...
        load_ca_sequence_coords(afPdbPath, afChain, ...
            'AltLocPolicy', 'blank_or_A', ...
            'Verbose', false);

    if strlength(seqT_ca) == 0 || strlength(seqA_ca) == 0
        error('Empty CA sequence. True CA=%d | AF CA=%d. Check chain IDs.', ...
            strlength(seqT_ca), strlength(seqA_ca));
    end

    % ---------------- load/filter official mapping ----------------
    M = readtable(mappingCSV, 'VariableNamingRule','preserve');

    required = {'PDB_ID','Chain_ID','PDB_resNum','PDB_iCode','UniProt_ID','UniProt_resNum','Observed'};
    for k = 1:numel(required)
        if ~any(strcmp(required{k}, M.Properties.VariableNames))
            error('Mapping CSV is missing required column: %s', required{k});
        end
    end

    keep = strcmpi(string(M.PDB_ID), pdbID) & strcmpi(string(M.Chain_ID), string(trueChain));

    if any(strcmp('AFDB_ID', M.Properties.VariableNames)) && strlength(afdbID) > 0
        keep_afdb = strcmpi(string(M.AFDB_ID), afdbID);
        if any(keep & keep_afdb)
            keep = keep & keep_afdb;
        end
    end

    if strlength(uniWanted) > 0
        keep = keep & strcmpi(base_uniprot(string(M.UniProt_ID)), base_uniprot(uniWanted));
    end

    obsVec = observed_to_logical(M.Observed);

    if logical(opt.RequireObserved)
        keep = keep & obsVec;
    end

    M = M(keep, :);

    if height(M) == 0
        error('No official mapping rows found for PDB=%s chain=%s UniProt=%s in %s.', ...
            pdbID, string(trueChain), uniWanted, mappingCSV);
    end

    M = sortrows(M, {'UniProt_resNum','PDB_resNum'});

    % ---------------- true PDB CA row lookup by author residue id ----------------
    trueKeyToRow = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(residT_ca)
        trueKeyToRow(char(residT_ca(i))) = i;
    end

    afResSeqToRow = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(resseqA_ca)
        key = char(string(resseqA_ca(i)));
        if ~isKey(afResSeqToRow, key)
            afResSeqToRow(key) = i;
        end
    end

    nRows = height(M);
    idxT = nan(nRows,1);
    idxA = nan(nRows,1);
    trueFullPos = nan(nRows,1);
    afFullPos = nan(nRows,1);

    for r = 1:nRows
        pdbKey = make_pdb_residue_key(M.PDB_resNum(r), M.PDB_iCode(r));
        if isKey(trueKeyToRow, char(pdbKey))
            idxT(r) = trueKeyToRow(char(pdbKey));
        end

        u = to_double_scalar(M.UniProt_resNum(r));
        if isnan(u)
            continue;
        end

        u = round(u);

        trueFullPos(r) = u;
        afFullPos(r) = u;

        % For AFDB monomer files, residue row usually equals UniProt position.
        % Fallback uses AF residue numbering if the row position is not valid.
        if u >= 1 && u <= numel(resseqA_ca)
            idxA(r) = u;
        else
            key = char(string(u));
            if isKey(afResSeqToRow, key)
                idxA(r) = afResSeqToRow(key);
            end
        end
    end

    ok = ~isnan(idxT) & ~isnan(idxA) & isfinite(trueFullPos) & isfinite(afFullPos);

    if ~all(ok) && logical(opt.Verbose)
        fprintf('[UNIPROT-MAP] Dropping %d mapped rows without both TRUE and AF CA rows.\n', nnz(~ok));
    end

    idxT = round(idxT(ok));
    idxA = round(idxA(ok));
    trueFullPos = round(trueFullPos(ok));
    afFullPos = round(afFullPos(ok));
    M_used = M(ok,:);

    % Remove accidental duplicate pairs, keeping first.
    pairMat = [idxT(:), idxA(:)];
    [~, ia] = unique(pairMat, 'rows', 'stable');
    idxT = idxT(ia);
    idxA = idxA(ia);
    trueFullPos = trueFullPos(ia);
    afFullPos = afFullPos(ia);
    M_used = M_used(ia,:);

    % Optional: drop official mapped residue pairs whose actual residue
    % letters differ between TRUE PDB and AlphaFold/UniProt.
    %
    % These are sequence conflicts, not necessarily alignment failures.
    nSequenceConflictsDropped = 0;
    sequenceConflictRowsDropped = M_used([],:);

    if logical(opt.PerfectMatchesOnly)
        cT_tmp = char(seqT_ca);
        cA_tmp = char(seqA_ca);

        keepPerfect = cT_tmp(idxT(:)) == cA_tmp(idxA(:));

        nSequenceConflictsDropped = nnz(~keepPerfect);
        sequenceConflictRowsDropped = M_used(~keepPerfect, :);

        if logical(opt.Verbose)
            fprintf('[UNIPROT-MAP] PerfectMatchesOnly: dropped %d / %d sequence-conflict mapped residues.\n', ...
                nSequenceConflictsDropped, numel(keepPerfect));

            if nSequenceConflictsDropped > 0
                bad = find(~keepPerfect(:));
                for rr = bad(:).'
                    fprintf('  conflict: UniProt %d | TRUE idxT %d = %c | AF idxA %d = %c\n', ...
                        M_used.UniProt_resNum(rr), ...
                        idxT(rr), cT_tmp(idxT(rr)), ...
                        idxA(rr), cA_tmp(idxA(rr)));
                end
            end
        end

        idxT = idxT(keepPerfect);
        idxA = idxA(keepPerfect);
        trueFullPos = trueFullPos(keepPerfect);
        afFullPos = afFullPos(keepPerfect);
        M_used = M_used(keepPerfect, :);
    end

    nPairs = numel(idxT);
    if nPairs < opt.MinAlignedPairs
        error('Too few official mapped CA pairs (%d). Check SIFTS mapping, PDB chain, or AFDB ID.', nPairs);
    end

    TRU = coordsT_ca(idxT, :);
    AF  = coordsA_ca(idxA, :);

    cT = char(seqT_ca);
    cA = char(seqA_ca);
    matchesCommon = sum(cT(idxT) == cA(idxA));
    identityCommon = matchesCommon / nPairs;

    % ---------------- full-position maps ----------------
    maxFull = max([trueFullPos(:); afFullPos(:); numel(resseqA_ca)]);
    maxFull = max(1, ceil(maxFull));

    trueFullToCa = zeros(maxFull,1);
    afFullToCa = zeros(maxFull,1);
    trueCaToFull = zeros(numel(resseqT_ca),1);
    afCaToFull = zeros(numel(resseqA_ca),1);

    for k = 1:nPairs
        fpT = trueFullPos(k);
        fpA = afFullPos(k);

        if fpT >= 1 && fpT <= maxFull
            trueFullToCa(fpT) = idxT(k);
        end

        if fpA >= 1 && fpA <= maxFull
            afFullToCa(fpA) = idxA(k);
        end

        trueCaToFull(idxT(k)) = fpT;
        afCaToFull(idxA(k)) = fpA;
    end

    % ---------------- trace breaks ----------------
    breakAfter_true = false(nPairs,1);
    breakAfter_af = false(nPairs,1);

    if nPairs > 1
        breakAfter_true(1:end-1) = ...
            (diff(idxT(:)) ~= 1) | ...
            (diff(double(resseqT_ca(idxT(:)))) ~= 1) | ...
            (diff(double(trueFullPos(:))) ~= 1);

        breakAfter_af(1:end-1) = ...
            (diff(idxA(:)) ~= 1) | ...
            (diff(double(resseqA_ca(idxA(:)))) ~= 1) | ...
            (diff(double(afFullPos(:))) ~= 1);
    end

    breakAfter_any = breakAfter_true | breakAfter_af;

    % ---------------- Procrustes baseline AF vs TRUE ----------------
    [d, Z, tform] = procrustes(TRU, AF, ...
        'Scaling', false, ...
        'Reflection', logical(opt.AllowReflection));

    rmsd = sqrt(mean(sum((Z - TRU).^2, 2)));

    % ---------------- reporting ----------------
    if logical(opt.Verbose)
        fprintf('\n[UNIPROT-MAP ALIGN]\n');
        fprintf('  Mapping source: %s\n', mappingCSV);
        fprintf('  PDB=%s chain=%s | AFDB=%s AF chain=%s | UniProt=%s\n', ...
            pdbID, string(trueChain), afdbID, string(afChain), uniWanted);
        fprintf('  CA loader: load_ca_sequence_coords -> load_ca_records -> load_pdb_atom_table\n');
        fprintf('  True CA residues modeled: %d (unknown X: %d)\n', strlength(seqT_ca), unknownT_ca);
        fprintf('  AF   CA residues:         %d (unknown X: %d)\n', strlength(seqA_ca), unknownA_ca);
        fprintf('  Official mapping rows used after CA intersection: %d\n', nPairs);
        fprintf('  Common CA identity: %.2f%%\n', 100*identityCommon);
        fprintf('  UniProt observed range used: %d..%d\n', min(trueFullPos), max(trueFullPos));
        fprintf('  Procrustes AF->TRUE: d=%.6f | RMSD=%.3f A | reflection=%d\n', d, rmsd, logical(opt.AllowReflection));
        fprintf('  Trace breaks true: %d | AF: %d | either: %d\n', ...
            nnz(breakAfter_true), nnz(breakAfter_af), nnz(breakAfter_any));

        if max(M_used.UniProt_resNum) > numel(resseqA_ca)
            fprintf(2, '  [WARN] Some mapped UniProt positions exceed AF CA count. Check isoform numbering.\n');
        end

        if identityCommon < 0.80
            fprintf(2, '  [WARN] Common CA identity below 80%%. Check target isoform and chain.\n');
        end
    end

    % ---------------- optional plot ----------------
    if logical(opt.MakeFigure)
        figure('Name','Official PDB-UniProt mapped common CA: True vs AF');

        subplot(1,2,1);
        plot3(TRU(:,1),TRU(:,2),TRU(:,3),'.'); hold on;
        plot3(AF(:,1), AF(:,2), AF(:,3), '.');
        axis equal; grid on; view(3);
        title('Before rigid alignment');
        legend('True common CA','AF common CA','Location','best');

        subplot(1,2,2);
        plot3(TRU(:,1),TRU(:,2),TRU(:,3),'.'); hold on;
        plot3(Z(:,1), Z(:,2), Z(:,3), '.');
        axis equal; grid on; view(3);
        title(sprintf('After alignment (RMSD %.3f A)', rmsd));
        legend('True common CA','AF aligned common CA','Location','best');
        rotate3d on;

        if strlength(string(opt.SaveFigureTo)) > 0
            saveas(gcf, string(opt.SaveFigureTo));
        end
    end

    % ---------------- output struct, backward-compatible ----------------
    out = struct();
    out.truePdbPath = string(truePdbPath);
    out.afPdbPath = string(afPdbPath);
    out.trueChain = string(trueChain);
    out.afChain = string(afChain);
    out.mapping_method = "uniprot_sifts";
    out.mapping_csv = mappingCSV;
    out.pdb_id = pdbID;
    out.afdb_id = afdbID;
    out.uniprot_id = uniWanted;
    out.mapping_rows_used = M_used;
    out.perfect_matches_only = logical(opt.PerfectMatchesOnly);
    out.n_sequence_conflicts_dropped = nSequenceConflictsDropped;
    out.sequence_conflict_rows_dropped = sequenceConflictRowsDropped;

    out.seq_true = seqT_ca;
    out.seq_af = seqA_ca;
    out.seq_true_full = "";
    out.seq_af_full = seqA_ca;
    out.true_seqres_resnames = strings(0,1);

    out.seq_true_ca = seqT_ca;
    out.seq_af_ca = seqA_ca;

    out.resid_true = residT_ca;
    out.resid_af = residA_ca;
    out.resseqs_true = resseqT_ca;
    out.resseqs_af = resseqA_ca;

    out.alignment_true = "official_pdb_uniprot_mapping";
    out.alignment_af = "official_pdb_uniprot_mapping";
    out.full_alignment_true = out.alignment_true;
    out.full_alignment_af = out.alignment_af;
    out.full_map_true = trueFullPos(:);
    out.full_map_af = afFullPos(:);

    out.true_full_to_ca_idx = trueFullToCa(:);
    out.true_ca_to_full_pos = trueCaToFull(:);
    out.af_full_to_ca_idx = afFullToCa(:);
    out.af_ca_to_full_pos = afCaToFull(:);

    out.idx_true_matched = idxT(:);
    out.idx_af_matched = idxA(:);

    out.true_fullpos_matched = trueFullPos(:);
    out.af_fullpos_matched = afFullPos(:);

    out.resid_true_matched = residT_ca(idxT);
    out.resid_af_matched = residA_ca(idxA);
    out.resseqs_true_matched = resseqT_ca(idxT);
    out.resseqs_af_matched = resseqA_ca(idxA);

    out.coords_true_matched = TRU;
    out.coords_af_matched = AF;
    out.coords_af_aligned = Z;

    out.identity = identityCommon;
    out.identity_common_ca = identityCommon;
    out.identity_full = identityCommon;
    out.n_aligned_pairs = nPairs;
    out.n_common_ca_pairs = nPairs;
    out.n_full_non_gap_pairs = nPairs;
    out.nw_score = NaN;

    out.breakAfter_true_matched = breakAfter_true;
    out.breakAfter_af_matched = breakAfter_af;
    out.breakAfter_any_matched = breakAfter_any;

    out.procrustes_d = d;
    out.tform = tform;
    out.rmsd = rmsd;

    out.unknown_true = unknownT_ca;
    out.unknown_af = unknownA_ca;
    out.extract_stats_true_ca = statsT_ca;
    out.extract_stats_af = statsA_ca;
end

% =====================================================================
function x = scalar_text(col, r)
    if iscell(col)
        x = string(col{r});
    else
        x = string(col(r));
    end
end

% =====================================================================
function obs = observed_to_logical(col)
    if islogical(col)
        obs = col;
    elseif isnumeric(col)
        obs = col ~= 0;
    else
        s = lower(strtrim(string(col)));
        obs = (s == "true") | (s == "1") | (s == "yes");
    end
end

% =====================================================================
function key = make_pdb_residue_key(resNum, iCode)
    rn = str2double(string(resNum));
    if isnan(rn)
        key = "";
        return;
    end

    ic = string(iCode);
    ic = strtrim(ic);

    if ismissing(ic) || strlength(ic) == 0 || strcmpi(ic, "NaN") || strcmpi(ic, "<missing>")
        key = string(sprintf('%d', rn));
    else
        key = string(sprintf('%d%s', rn, char(ic)));
    end
end

% =====================================================================
function b = base_uniprot(acc)
    b = string(acc);
    b = regexprep(b, '^AF-', '');
    b = regexprep(b, '-F\d+$', '');
    b = regexprep(b, '-.*$', '');
end

% =====================================================================
function x = to_double_scalar(v)
    if isempty(v)
        x = NaN;
        return;
    end

    if isnumeric(v)
        x = double(v(1));
        return;
    end

    s = string(v);
    if isempty(s) || ismissing(s(1))
        x = NaN;
        return;
    end

    x = str2double(s(1));
end