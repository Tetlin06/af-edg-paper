function [CA, stats] = load_ca_records(pdbPathOrID, chainID, varargin)
% LOAD_CA_RECORDS
% -------------------------------------------------------------------------
% CA-specific wrapper around load_pdb_atom_table.m.
%
% This function is the shared source for all C-alpha residue-level loading.
% It should be used by:
%   - load_ca_plddt.m
%   - load_ca_sequence_coords.m
%   - alignment functions that need CA sequence/coords
%
% It does NOT read the PDB directly. It calls load_pdb_atom_table.m.
%
% Outputs:
%   CA    : table with one row per residue CA atom
%   stats : struct with useful counts
%
% Important:
%   - Keeps only ATOM records.
%   - Keeps only atomName == "CA".
%   - Filters to the selected chain.
%   - Default altLoc behavior keeps blank or A.
%   - If duplicate CA records remain, keeps one per residue, preferring:
%       blank altLoc > altLoc A > higher occupancy > earlier sourceIndex.
%
% CA table columns inherited from load_pdb_atom_table include:
%   atomName, altLoc, resName, chainID, resSeq, iCode, resid, chainResid,
%   x, y, z, occupancy, tempFactor, element
%
% Additional columns added here:
%   aa1    : one-letter residue code
%   coords : [x y z] coordinate row
%
% Example:
%   [CA, stats] = load_ca_records('data/AFDB/P52799.pdb', 'A');
%   height(CA)
%   CA.coords(1:5,:)

    if nargin < 2
        chainID = [];
    end

    ip = inputParser;
    ip.FunctionName = mfilename;

    ip.addParameter('ModelIndex', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1 && floor(x) == x);
    ip.addParameter('AltLocPolicy', "blank_or_A", @(x)ischar(x) || isstring(x));
    ip.addParameter('Deduplicate', true, @(x)islogical(x) || isnumeric(x));
    ip.addParameter('Verbose', false, @(x)islogical(x) || isnumeric(x));

    ip.parse(varargin{:});
    opt = ip.Results;

    verbose = logical(opt.Verbose);
    doDedup = logical(opt.Deduplicate);

    % ---------------------------------------------------------------------
    % 1) Load full atom table from the shared parser
    % ---------------------------------------------------------------------
    T = load_pdb_atom_table( ...
        pdbPathOrID, ...
        'ChainID', chainID, ...
        'ModelIndex', opt.ModelIndex, ...
        'IncludeHetero', false, ...
        'AltLocPolicy', opt.AltLocPolicy, ...
        'Verbose', verbose);

    stats = struct();
    stats.nAtomRowsLoaded = height(T);

    % ---------------------------------------------------------------------
    % 2) Keep only ATOM C-alpha records with usable coordinates/residue nums
    % ---------------------------------------------------------------------
    isCA = upper(T.atomName) == "CA";
    isAtom = T.recordType == "ATOM";

    hasCoord = isfinite(T.x) & isfinite(T.y) & isfinite(T.z);
    hasResSeq = isfinite(T.resSeq);

    CA = T(isAtom & isCA & hasCoord & hasResSeq, :);

    stats.nCA_before_dedup = height(CA);
    stats.nDropped_nonCA_or_unusable = height(T) - height(CA);

    if isempty(CA)
        error('No usable CA atoms found for requested PDB/chain.');
    end

    % ---------------------------------------------------------------------
    % 3) Deduplicate one CA per residue if needed
    % ---------------------------------------------------------------------
    if doDedup
        [CA, nDuplicateRowsDropped] = deduplicate_ca_rows(CA);
    else
        nDuplicateRowsDropped = 0;
    end

    stats.nDuplicateRowsDropped = nDuplicateRowsDropped;
    stats.nCA_after_dedup = height(CA);

    % ---------------------------------------------------------------------
    % 4) Sort in biological residue order
    % ---------------------------------------------------------------------
    CA = sortrows(CA, {'resSeq', 'iCode', 'sourceIndex'});

    % ---------------------------------------------------------------------
    % 5) Add residue sequence and coordinate convenience columns
    % ---------------------------------------------------------------------
    aa1 = strings(height(CA), 1);
    unknownCount = 0;

    for i = 1:height(CA)
        aa = aa3_to_1(CA.resName(i));
        aa1(i) = string(aa);

        if aa == 'X'
            unknownCount = unknownCount + 1;
        end
    end

    CA.aa1 = aa1;
    CA.coords = [CA.x, CA.y, CA.z];

    stats.nResidues = height(CA);
    stats.unknownCount = unknownCount;

    if isempty(chainID)
        stats.chainID = "";
    else
        stats.chainID = string(chainID);
    end

    if verbose
        fprintf('\n[load_ca_records]\n');
        fprintf('  CA before dedup: %d\n', stats.nCA_before_dedup);
        fprintf('  Duplicate CA rows dropped: %d\n', stats.nDuplicateRowsDropped);
        fprintf('  CA after dedup: %d\n', stats.nCA_after_dedup);
        fprintf('  Unknown residue codes X: %d\n', stats.unknownCount);
        fprintf('  Residue range: %s..%s\n', CA.resid(1), CA.resid(end));
    end
end

% =========================================================================
function [CAout, nDropped] = deduplicate_ca_rows(CA)
% Keep one CA row per chainResid.
%
% Preference:
%   blank altLoc > altLoc A > higher occupancy > earlier sourceIndex

    keys = CA.chainResid;
    uniqueKeys = unique(keys, 'stable');

    keep = false(height(CA), 1);

    for k = 1:numel(uniqueKeys)
        idx = find(keys == uniqueKeys(k));

        if numel(idx) == 1
            keep(idx) = true;
            continue;
        end

        bestLocal = choose_best_ca_row(CA(idx, :));
        keep(idx(bestLocal)) = true;
    end

    CAout = CA(keep, :);
    nDropped = height(CA) - height(CAout);
end

% =========================================================================
function bestLocal = choose_best_ca_row(CAsub)
% Return row index within CAsub.

    n = height(CAsub);

    score = zeros(n, 1);

    for i = 1:n
        alt = upper(strtrim(string(CAsub.altLoc(i))));

        if alt == ""
            score(i) = score(i) + 1000;
        elseif alt == "A"
            score(i) = score(i) + 500;
        end

        occ = CAsub.occupancy(i);
        if isfinite(occ)
            score(i) = score(i) + occ;
        end

        src = CAsub.sourceIndex(i);
        if isfinite(src)
            score(i) = score(i) - 1e-9 * src;
        end
    end

    [~, bestLocal] = max(score);
end

% =========================================================================
function aa = aa3_to_1(resName)
% Convert 3-letter residue name to 1-letter residue code.
% Unknown or unsupported residues become X.
%
% Includes common modified residue mappings:
%   MSE -> M
%   SEC -> U
%   PYL -> O
%   ASX -> B
%   GLX -> Z

    r = upper(strtrim(char(string(resName))));

    switch r
        case 'ALA'
            aa = 'A';
        case 'ARG'
            aa = 'R';
        case 'ASN'
            aa = 'N';
        case 'ASP'
            aa = 'D';
        case 'CYS'
            aa = 'C';
        case 'GLN'
            aa = 'Q';
        case 'GLU'
            aa = 'E';
        case 'GLY'
            aa = 'G';
        case 'HIS'
            aa = 'H';
        case 'ILE'
            aa = 'I';
        case 'LEU'
            aa = 'L';
        case 'LYS'
            aa = 'K';
        case {'MET', 'MSE'}
            aa = 'M';
        case 'PHE'
            aa = 'F';
        case 'PRO'
            aa = 'P';
        case 'SER'
            aa = 'S';
        case 'THR'
            aa = 'T';
        case 'TRP'
            aa = 'W';
        case 'TYR'
            aa = 'Y';
        case 'VAL'
            aa = 'V';
        case 'SEC'
            aa = 'U';
        case 'PYL'
            aa = 'O';
        case 'ASX'
            aa = 'B';
        case 'GLX'
            aa = 'Z';
        otherwise
            aa = 'X';
    end
end