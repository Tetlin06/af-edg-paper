function [DistSq_out, info] = apply_chemistry_aware_noise( ...
    DistSq_in, Weight, K, noiseSeed, edgNodeLevel, atomMeta, protectChemistry, varargin)
% APPLY_CHEMISTRY_AWARE_NOISE
% Adds multiplicative noise only to non-chemistry graph edges.
%
% Chemistry/backbone-like edges are kept clean:
%
%   same residue:
%       N-CA, CA-C, C-O, CA-CB
%
%   adjacent residues:
%       C_i-N_{i+1}, CA_i-CA_{i+1}
%       CA_i-N_{i+1}, C_i-CA_{i+1}, O_i-N_{i+1}, O_i-CA_{i+1}
%
% In CA mode, only adjacent CA_i-CA_{i+1} edges are protected.
%
% Usage:
%   [DistSq, info] = apply_chemistry_aware_noise( ...
%       DistSq, Weight, K, noiseSeed, edgNodeLevel, atomMeta, protectChemistry);
%
% Optional:
%   'Verbose', true/false
%
% Notes:
%   - Silent by default.
%   - Only edges with Weight ~= 0 are considered.
%   - If K = 0, the matrix is returned unchanged.
%   - If protectChemistry = false, all graph edges receive noise.

    % ---------------------------------------------------------------------
    % Defaults
    % ---------------------------------------------------------------------
    if nargin < 7 || isempty(protectChemistry)
        protectChemistry = true;
    end

    if nargin < 6
        atomMeta = [];
    end

    if nargin < 5 || isempty(edgNodeLevel)
        edgNodeLevel = "ca";
    end

    if nargin < 4
        noiseSeed = [];
    end

    ip = inputParser;
    ip.addParameter('Verbose', false, @(x)islogical(x) || isnumeric(x));
    ip.parse(varargin{:});
    verbose = logical(ip.Results.Verbose);

    % ---------------------------------------------------------------------
    % Basic checks
    % ---------------------------------------------------------------------
    if ~ismatrix(DistSq_in) || size(DistSq_in,1) ~= size(DistSq_in,2)
        error('DistSq_in must be square.');
    end

    N = size(DistSq_in, 1);

    if any(size(Weight) ~= [N N])
        error('Weight must have the same size as DistSq_in.');
    end

    if ~isscalar(K) || ~isnumeric(K) || ~isfinite(K) || K < 0
        error('K must be a nonnegative finite scalar.');
    end

    if any(~isfinite(DistSq_in(:)))
        error('DistSq_in contains NaN or Inf.');
    end

    % ---------------------------------------------------------------------
    % Clean / normalize graph mask
    % ---------------------------------------------------------------------
    W = double(Weight ~= 0);
    W = double((W + W') > 0);
    W(1:N+1:end) = 0;

    edgeMask = logical(W);

    % ---------------------------------------------------------------------
    % Initialize output and info
    % ---------------------------------------------------------------------
    DistSq_out = DistSq_in;

    info = struct();
    info.K = K;
    info.noiseSeed = noiseSeed;
    info.edgNodeLevel = string(edgNodeLevel);
    info.protectChemistry = logical(protectChemistry);
    info.status = "clean";
    info.nEdges = nnz(triu(edgeMask, 1));
    info.nChemistryEdges = 0;
    info.nNoisyEdges = 0;
    info.fracChemistryEdges = 0;
    info.fracNoisyEdges = 0;

    if K == 0
        if verbose
            fprintf('[NOISE] Clean matrix | K=0\n');
        end
        return;
    end

    % ---------------------------------------------------------------------
    % Build noise matrix
    % ---------------------------------------------------------------------
    if ~isempty(noiseSeed)
        rng(noiseSeed, 'twister');
    end

    Noise = randn(N, N);
    Noise = 0.5 * (Noise + Noise');
    Noise(1:N+1:end) = 0;

    % ---------------------------------------------------------------------
    % Build chemistry-protection mask
    % ---------------------------------------------------------------------
    if logical(protectChemistry)
        chemMask = make_chemistry_mask(W, edgNodeLevel, atomMeta);
    else
        chemMask = false(N, N);
    end

    chemMask = chemMask & edgeMask;
    chemMask = chemMask | chemMask.';
    chemMask(1:N+1:end) = false;

    noiseMask = edgeMask & ~chemMask;

    % ---------------------------------------------------------------------
    % Apply noise only to non-chemistry graph edges
    % ---------------------------------------------------------------------
    DistSq_out(noiseMask) = DistSq_in(noiseMask) + ...
        K * Noise(noiseMask) .* DistSq_in(noiseMask);

    % Restore exact chemistry/backbone edges.
    DistSq_out(chemMask) = DistSq_in(chemMask);

    % Keep matrix symmetric and hollow.
    DistSq_out = 0.5 * (DistSq_out + DistSq_out');
    DistSq_out(1:N+1:end) = 0;

    % ---------------------------------------------------------------------
    % Fill diagnostics
    % ---------------------------------------------------------------------
    info.status = "noisy_nonchemistry_only";
    info.nChemistryEdges = nnz(triu(chemMask, 1));
    info.nNoisyEdges = nnz(triu(noiseMask, 1));
    info.fracChemistryEdges = info.nChemistryEdges / max(info.nEdges, 1);
    info.fracNoisyEdges = info.nNoisyEdges / max(info.nEdges, 1);

    if verbose
        fprintf('[NOISE] Chemistry-aware multiplicative noise | K=%.4g\n', K);

        if ~isempty(noiseSeed)
            fprintf('        noise seed: %d\n', noiseSeed);
        end

        fprintf('        graph edges:                %d\n', info.nEdges);
        fprintf('        protected chemistry edges:  %d\n', info.nChemistryEdges);
        fprintf('        noisy non-chemistry edges:  %d\n', info.nNoisyEdges);
    end
end

% =========================================================================
function chemMask = make_chemistry_mask(Weight, edgNodeLevel, atomMeta)

    N = size(Weight, 1);
    chemMask = false(N, N);

    mode = lower(string(edgNodeLevel));

    % ---------------------------------------------------------------------
    % CA-only fallback
    % ---------------------------------------------------------------------
    if mode == "ca"
        for i = 1:(N - 1)
            chemMask(i, i+1) = true;
            chemMask(i+1, i) = true;
        end

        chemMask = chemMask & logical(Weight ~= 0);
        chemMask = chemMask | chemMask.';
        chemMask(1:N+1:end) = false;
        return;
    end

    % ---------------------------------------------------------------------
    % If not all-atom, quietly protect nothing.
    % Paper-only runner protects nothing for unknown node modes.
    % ---------------------------------------------------------------------
    if mode ~= "all_atom"
        return;
    end

    if isempty(atomMeta) || ~istable(atomMeta)
        return;
    end

    required = ["residuePairIndex", "atomName"];
    for k = 1:numel(required)
        if ~ismember(required(k), string(atomMeta.Properties.VariableNames))
            error('atomMeta is missing required column: %s', required(k));
        end
    end

    if height(atomMeta) ~= N
        error('atomMeta height must match matrix size.');
    end

    atomName = upper(strtrim(string(atomMeta.atomName)));
    resIdx = double(atomMeta.residuePairIndex(:));

    if ismember('trueResSeq', atomMeta.Properties.VariableNames)
        trueResSeq = double(atomMeta.trueResSeq(:));
    else
        trueResSeq = resIdx;
    end

    if ismember('alignedPairPos', atomMeta.Properties.VariableNames)
        alignedPairPos = double(atomMeta.alignedPairPos(:));
    else
        alignedPairPos = resIdx;
    end

    % ---------------------------------------------------------------------
    % Same-residue covalent/local backbone geometry
    % ---------------------------------------------------------------------
    sameResiduePairs = [
        "N",  "CA"
        "CA", "C"
        "C",  "O"
        "CA", "CB"
    ];

    for p = 1:size(sameResiduePairs, 1)
        chemMask = add_same_residue_pair( ...
            chemMask, atomName, resIdx, ...
            sameResiduePairs(p,1), sameResiduePairs(p,2));
    end

    % ---------------------------------------------------------------------
    % Adjacent-residue peptide/backbone geometry
    % ---------------------------------------------------------------------
    adjacentPairs = [
        "C",  "N"     % peptide bond
        "CA", "CA"    % CA_i to CA_{i+1}
        "CA", "N"     % peptide-plane helper
        "C",  "CA"    % peptide-plane helper
        "O",  "N"     % peptide-plane helper
        "O",  "CA"    % peptide-plane helper
    ];

    uniqueResidues = unique(resIdx, 'stable');

    for rr = 1:(numel(uniqueResidues) - 1)
        r1 = uniqueResidues(rr);
        r2 = uniqueResidues(rr + 1);

        rows1 = find(resIdx == r1);
        rows2 = find(resIdx == r2);

        if isempty(rows1) || isempty(rows2)
            continue;
        end

        % Avoid fake peptide links across gaps/domain cuts.
        seqOK = trueResSeq(rows2(1)) == trueResSeq(rows1(1)) + 1;
        alignOK = alignedPairPos(rows2(1)) == alignedPairPos(rows1(1)) + 1;

        if ~(seqOK && alignOK)
            continue;
        end

        for p = 1:size(adjacentPairs, 1)
            chemMask = add_adjacent_residue_pair( ...
                chemMask, atomName, resIdx, r1, r2, ...
                adjacentPairs(p,1), adjacentPairs(p,2));
        end
    end

    chemMask = chemMask & logical(Weight ~= 0);
    chemMask = chemMask | chemMask.';
    chemMask(1:N+1:end) = false;
end

% =========================================================================
function M = add_same_residue_pair(M, atomName, resIdx, atomA, atomB)

    residues = unique(resIdx(:).', 'stable');

    for r = residues
        ia = find(resIdx == r & atomName == atomA, 1, 'first');
        ib = find(resIdx == r & atomName == atomB, 1, 'first');

        if ~isempty(ia) && ~isempty(ib)
            M(ia, ib) = true;
            M(ib, ia) = true;
        end
    end
end

% =========================================================================
function M = add_adjacent_residue_pair(M, atomName, resIdx, r1, r2, atomA, atomB)

    ia = find(resIdx == r1 & atomName == atomA, 1, 'first');
    ib = find(resIdx == r2 & atomName == atomB, 1, 'first');

    if ~isempty(ia) && ~isempty(ib)
        M(ia, ib) = true;
        M(ib, ia) = true;
    end
end