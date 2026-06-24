function [AT, stats] = load_all_atom_records(pdbPathOrID, chainID, varargin)
% LOAD_ALL_ATOM_RECORDS
% -------------------------------------------------------------------------
% All-atom wrapper around load_pdb_atom_table.m.
%
% This is the all-atom analogue of load_ca_records.m.
%
% Purpose:
%   Build a clean atom-level table for future all-atom EDG experiments.
%
% Important:
%   - This does NOT replace load_ca_records.m yet.
%   - This does NOT modify main.m yet.
%   - This does NOT solve atom matching yet.
%   - This only creates a clean atom table for one PDB chain.
%
% Default behavior:
%   - Keeps ATOM records only.
%   - Excludes HETATM by default.
%   - Removes hydrogens by default.
%   - Keeps blank altLoc or altLoc A.
%   - Deduplicates alternate atoms by residue + atom name.
%
% Outputs:
%   AT    : table with one row per atom
%   stats : summary struct
%
% Important AT columns:
%   atomIndex       : row index in this cleaned all-atom table
%   residueIndex    : residue/block index inside this cleaned table
%   atomKey         : unique atom key, e.g. "A:70:CA"
%   residueAtomKey  : residue-local atom key, e.g. "70A:CA"
%   coords          : [x y z]
%   isCA            : true for CA atom
%   isBackbone      : true for N, CA, C, O, OXT
%   isHeavyAtom     : true for non-hydrogen atoms
%
% Example:
%   [AT, stats] = load_all_atom_records('data/AFDB/P52799.pdb','A','Verbose',true);

    if nargin < 2
        chainID = [];
    end

    ip = inputParser;
    ip.FunctionName = mfilename;

    ip.addParameter('ModelIndex', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1 && floor(x) == x);
    ip.addParameter('IncludeHetero', false, @(x)islogical(x) || isnumeric(x));
    ip.addParameter('AltLocPolicy', "blank_or_A", @(x)ischar(x) || isstring(x));
    ip.addParameter('RemoveHydrogens', true, @(x)islogical(x) || isnumeric(x));
    ip.addParameter('Deduplicate', true, @(x)islogical(x) || isnumeric(x));
    ip.addParameter('Verbose', false, @(x)islogical(x) || isnumeric(x));

    ip.parse(varargin{:});
    opt = ip.Results;

    includeHetero = logical(opt.IncludeHetero);
    removeHydrogens = logical(opt.RemoveHydrogens);
    doDedup = logical(opt.Deduplicate);
    verbose = logical(opt.Verbose);

    % ---------------------------------------------------------------------
    % 1) Load full atom table from the shared parser
    % ---------------------------------------------------------------------
    T = load_pdb_atom_table( ...
        pdbPathOrID, ...
        'ChainID', chainID, ...
        'ModelIndex', opt.ModelIndex, ...
        'IncludeHetero', includeHetero, ...
        'AltLocPolicy', opt.AltLocPolicy, ...
        'Verbose', verbose);

    stats = struct();
    stats.nAtomRowsLoaded = height(T);

    % ---------------------------------------------------------------------
    % 2) Keep requested record types and usable coordinates
    % ---------------------------------------------------------------------
    T.recordType = upper(strtrim(string(T.recordType)));
    T.atomName   = upper(strtrim(string(T.atomName)));
    T.resName    = upper(strtrim(string(T.resName)));
    T.chainID    = strtrim(string(T.chainID));
    T.altLoc     = strtrim(string(T.altLoc));
    T.iCode      = strtrim(string(T.iCode));
    T.resid      = strtrim(string(T.resid));
    T.chainResid = strtrim(string(T.chainResid));
    T.element    = upper(strtrim(string(T.element)));

    if includeHetero
        keepRecord = (T.recordType == "ATOM") | (T.recordType == "HETATM");
    else
        keepRecord = (T.recordType == "ATOM");
    end

    hasCoord = isfinite(T.x) & isfinite(T.y) & isfinite(T.z);
    hasResSeq = isfinite(T.resSeq);
    hasAtomName = T.atomName ~= "";

    AT = T(keepRecord & hasCoord & hasResSeq & hasAtomName, :);

    stats.nRowsAfterRecordCoordFilter = height(AT);
    stats.nDropped_recordCoordAtomName = height(T) - height(AT);

    if isempty(AT)
        error('No usable atom records found after record/coordinate filtering.');
    end

    % ---------------------------------------------------------------------
    % 3) Remove hydrogens by default
    % ---------------------------------------------------------------------
    if removeHydrogens
        isH = infer_hydrogen_mask(AT);
        stats.nHydrogensDropped = nnz(isH);
        AT = AT(~isH, :);
    else
        stats.nHydrogensDropped = 0;
    end

    if isempty(AT)
        error('No atoms remain after hydrogen filtering.');
    end

    % ---------------------------------------------------------------------
    % 4) Build atom keys before deduplication
    % ---------------------------------------------------------------------
    AT.atomKey = AT.chainResid + ":" + AT.atomName;
    AT.residueAtomKey = AT.resid + ":" + AT.atomName;

    % ---------------------------------------------------------------------
    % 5) Deduplicate alternate atoms if needed
    % ---------------------------------------------------------------------
    if doDedup
        [AT, nDuplicateRowsDropped] = deduplicate_atom_rows(AT);
    else
        nDuplicateRowsDropped = 0;
    end

    stats.nDuplicateRowsDropped = nDuplicateRowsDropped;

    % ---------------------------------------------------------------------
    % 6) Sort in residue/source order
    % ---------------------------------------------------------------------
    AT = sortrows(AT, {'resSeq', 'iCode', 'sourceIndex'});

    % ---------------------------------------------------------------------
    % 7) Add atom/residue indices and useful flags
    % ---------------------------------------------------------------------
    AT.atomIndex = (1:height(AT)).';
    AT.coords = [AT.x, AT.y, AT.z];

    uniqueResidues = unique(AT.chainResid, 'stable');

    residueIndex = zeros(height(AT), 1);
    for r = 1:numel(uniqueResidues)
        residueIndex(AT.chainResid == uniqueResidues(r)) = r;
    end

    AT.residueIndex = residueIndex;

    AT.isCA = upper(AT.atomName) == "CA";
    AT.isBackbone = ismember(upper(AT.atomName), ["N", "CA", "C", "O", "OXT"]);
    AT.isHeavyAtom = ~infer_hydrogen_mask(AT);

    % Move convenience columns near the front if possible
    frontVars = {'atomIndex','residueIndex','atomKey','residueAtomKey'};
    AT = movevars(AT, frontVars, 'Before', 1);

    % ---------------------------------------------------------------------
    % 8) Stats
    % ---------------------------------------------------------------------
    stats.nAtoms = height(AT);
    stats.nResidues = numel(uniqueResidues);
    stats.nCA = nnz(AT.isCA);
    stats.nBackbone = nnz(AT.isBackbone);
    stats.nHeavyAtom = nnz(AT.isHeavyAtom);
    stats.includeHetero = includeHetero;
    stats.removeHydrogens = removeHydrogens;
    stats.altLocPolicy = string(opt.AltLocPolicy);

    if isempty(chainID)
        stats.chainID = "";
    else
        stats.chainID = string(chainID);
    end

    if verbose
        fprintf('\n[load_all_atom_records]\n');
        fprintf('  Atoms after cleaning: %d\n', stats.nAtoms);
        fprintf('  Residues represented: %d\n', stats.nResidues);
        fprintf('  CA atoms: %d\n', stats.nCA);
        fprintf('  Backbone atoms: %d\n', stats.nBackbone);
        fprintf('  Hydrogens dropped: %d\n', stats.nHydrogensDropped);
        fprintf('  Duplicate atom rows dropped: %d\n', stats.nDuplicateRowsDropped);
        fprintf('  Include HETATM: %d\n', stats.includeHetero);
        fprintf('  Residue range: %s..%s\n', AT.resid(1), AT.resid(end));
    end
end

% =========================================================================
function [ATout, nDropped] = deduplicate_atom_rows(AT)
% Keep one atom row per atomKey.
%
% Preference:
%   blank altLoc > altLoc A > higher occupancy > earlier sourceIndex

    keys = AT.atomKey;
    uniqueKeys = unique(keys, 'stable');

    keep = false(height(AT), 1);

    for k = 1:numel(uniqueKeys)
        idx = find(keys == uniqueKeys(k));

        if numel(idx) == 1
            keep(idx) = true;
            continue;
        end

        bestLocal = choose_best_atom_row(AT(idx, :));
        keep(idx(bestLocal)) = true;
    end

    ATout = AT(keep, :);
    nDropped = height(AT) - height(ATout);
end

% =========================================================================
function bestLocal = choose_best_atom_row(ATsub)
% Return row index within ATsub.

    n = height(ATsub);
    score = zeros(n, 1);

    for i = 1:n
        alt = upper(strtrim(string(ATsub.altLoc(i))));

        if alt == ""
            score(i) = score(i) + 1000;
        elseif alt == "A"
            score(i) = score(i) + 500;
        end

        occ = ATsub.occupancy(i);
        if isfinite(occ)
            score(i) = score(i) + occ;
        end

        src = ATsub.sourceIndex(i);
        if isfinite(src)
            score(i) = score(i) - 1e-9 * src;
        end
    end

    [~, bestLocal] = max(score);
end

% =========================================================================
function isH = infer_hydrogen_mask(T)
% Infer hydrogen atoms.
%
% Prefer element column when available.
% If element is missing, infer from atomName by taking the first alphabetic
% character. For protein atom names, CA/CB/CD/etc. correctly infer as carbon.

    n = height(T);
    isH = false(n, 1);

    for i = 1:n
        elem = upper(strtrim(string(T.element(i))));

        if elem == ""
            elem = infer_element_from_atom_name(T.atomName(i));
        end

        isH(i) = (elem == "H");
    end
end

% =========================================================================
function elem = infer_element_from_atom_name(atomName)
    s = upper(strtrim(char(string(atomName))));

    if isempty(s)
        elem = "";
        return;
    end

    letters = regexp(s, '[A-Z]', 'match');

    if isempty(letters)
        elem = "";
    else
        elem = string(letters{1});
    end
end