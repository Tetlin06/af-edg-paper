function sub = subset_aligned_pairs(out_align, sel)
% SUBSET_ALIGNED_PAIRS
% Select a subset of matched alignment pairs after sequence alignment.
%
% Supported forms for sel:
%   []                 -> keep ALL matched pairs
%   scalar N           -> keep 1:N
%   [a b]              -> keep a:b
%   [a1 b1; a2 b2]     -> keep UNION of multiple ranges
%   logical mask       -> keep where mask==true, length = nPairs
%   index vector       -> explicit indices into matched pairs
%
% Notes:
% - Ranges are in aligned-pair index space, NOT raw residue numbering.
% - This version preserves full-sequence mapping fields when present.
% - It also computes trace-break fields for matched-trace bookkeeping.
%   A trace should break whenever the selected rows are not continuous in
%   TRUE, in AF, or because the user selected a non-contiguous subset.

    if ~isfield(out_align, 'n_aligned_pairs')
        error('out_align must contain field n_aligned_pairs.');
    end

    nPairs = out_align.n_aligned_pairs;

    % ---- normalize sel ----
    if nargin < 2 || isempty(sel)
        pos = 1:nPairs;

    elseif isscalar(sel)
        pos = 1:min(sel, nPairs);

    elseif islogical(sel)
        if numel(sel) ~= nPairs
            error('Logical mask must have length %d.', nPairs);
        end
        pos = find(sel);

    elseif isnumeric(sel) && size(sel,2) == 2 && size(sel,1) >= 1
        % Single range [a b] OR multiple ranges [a1 b1; a2 b2; ...]
        pos = [];

        for r = 1:size(sel,1)
            a = sel(r,1);
            b = sel(r,2);

            if ~isscalar(a) || ~isscalar(b) || ~isnumeric(a) || ~isnumeric(b)
                error('Each range row must be numeric [a b].');
            end

            if isinf(a)
                error('Lower bound cannot be Inf.');
            end
            if isinf(b)
                b = nPairs;
            end

            a = max(1, ceil(a));
            b = min(nPairs, floor(b));

            if a > b
                error('Invalid range [%g %g] with nPairs=%d.', sel(r,1), sel(r,2), nPairs);
            end

            pos = [pos, a:b]; 
        end

        pos = unique(pos, 'stable');

    else
        % Explicit index vector
        pos = sel(:)';
    end

    if isempty(pos)
        error('Subset is empty after selection.');
    end
    if any(pos < 1) || any(pos > nPairs)
        error('Subset indices must lie in 1..%d.', nPairs);
    end

    pos = unique(pos(:), 'stable');
    N = numel(pos);

    % ---- required fields ----
    requiredFields = { ...
        'idx_true_matched', ...
        'idx_af_matched', ...
        'coords_true_matched', ...
        'coords_af_matched', ...
        'resseqs_true_matched', ...
        'resseqs_af_matched', ...
        'resid_true_matched', ...
        'resid_af_matched'};

    for k = 1:numel(requiredFields)
        if ~isfield(out_align, requiredFields{k})
            error('out_align is missing required field: %s', requiredFields{k});
        end
    end

    % ---- build output ----
    sub = struct();
    sub.pos = pos(:);

    sub.idxT = out_align.idx_true_matched(pos);
    sub.idxA = out_align.idx_af_matched(pos);

    sub.coords_true = out_align.coords_true_matched(pos, :);
    sub.coords_af   = out_align.coords_af_matched(pos, :);

    sub.resseq_true = out_align.resseqs_true_matched(pos);
    sub.resseq_af   = out_align.resseqs_af_matched(pos);

    sub.resid_true  = out_align.resid_true_matched(pos);
    sub.resid_af    = out_align.resid_af_matched(pos);

    % ---- preserve full-sequence alignment fields when available ----
    if isfield(out_align, 'true_fullpos_matched')
        sub.true_fullpos = out_align.true_fullpos_matched(pos);
        sub.true_fullpos_matched = sub.true_fullpos;
    end

    if isfield(out_align, 'af_fullpos_matched')
        sub.af_fullpos = out_align.af_fullpos_matched(pos);
        sub.af_fullpos_matched = sub.af_fullpos;
    end

    % Optional aliases / metadata useful for debugging.
    if isfield(out_align, 'identity')
        sub.identity_source = out_align.identity;
    end
    if isfield(out_align, 'identity_common_ca')
        sub.identity_common_ca_source = out_align.identity_common_ca;
    end
    if isfield(out_align, 'identity_full')
        sub.identity_full_source = out_align.identity_full;
    end
    if isfield(out_align, 'truePdbPath')
        sub.truePdbPath = out_align.truePdbPath;
    end
    if isfield(out_align, 'afPdbPath')
        sub.afPdbPath = out_align.afPdbPath;
    end
    if isfield(out_align, 'trueChain')
        sub.trueChain = out_align.trueChain;
    end
    if isfield(out_align, 'afChain')
        sub.afChain = out_align.afChain;
    end

    % ---- recompute trace breaks for this selected subset ----
    [breakAfter_true, breakAfter_af] = infer_subset_breaks(out_align, sub, pos, nPairs);
    breakAfter_any = breakAfter_true | breakAfter_af;

    sub.breakAfter_true = breakAfter_true;
    sub.breakAfter_af   = breakAfter_af;
    sub.breakAfter_any  = breakAfter_any;

    % Backward/forward-compatible aliases.
    sub.breakAfter_true_matched = breakAfter_true;
    sub.breakAfter_af_matched   = breakAfter_af;
    sub.breakAfter_any_matched  = breakAfter_any;

    % Paper-only public runner is intentionally quiet.
end

% =========================================================================
function [breakAfter_true, breakAfter_af] = infer_subset_breaks(out_align, sub, pos, nPairs)
% Build break vectors for the selected subset.
%
% breakAfter_true(i)=true means insert TER after row i for true-side
% continuity. breakAfter_af is the same idea for AF-side continuity.
%
% A break is inserted if:
%   - selected aligned-pair positions are not consecutive, OR
%   - full-sequence positions jump, OR
%   - CA-row indices jump, OR
%   - author residue numbers jump, OR
%   - source alignment already marked a break.

    N = numel(pos);
    breakAfter_true = false(N,1);
    breakAfter_af   = false(N,1);

    if N <= 1
        return;
    end

    pos = pos(:);

    % If the user selected non-contiguous aligned-pair rows, do not draw
    % a continuous trace across the omitted region.
    selectedGap = diff(pos) ~= 1;

    % ---- TRUE continuity ----
    bT = selectedGap;

    if isfield(sub, 'true_fullpos') && numel(sub.true_fullpos) == N
        bT = bT | (diff(double(sub.true_fullpos(:))) ~= 1);
    end

    if isfield(sub, 'idxT') && numel(sub.idxT) == N
        bT = bT | (diff(double(sub.idxT(:))) ~= 1);
    end

    if isfield(sub, 'resseq_true') && numel(sub.resseq_true) == N
        bT = bT | (diff(double(sub.resseq_true(:))) ~= 1);
    end

    if isfield(out_align, 'breakAfter_true_matched')
        src = logical(out_align.breakAfter_true_matched(:));
        bT = bT | source_breaks_for_selected_positions(src, pos, N, nPairs);
    end

    breakAfter_true(1:end-1) = bT;

    % ---- AF continuity ----
    bA = selectedGap;

    if isfield(sub, 'af_fullpos') && numel(sub.af_fullpos) == N
        bA = bA | (diff(double(sub.af_fullpos(:))) ~= 1);
    end

    if isfield(sub, 'idxA') && numel(sub.idxA) == N
        bA = bA | (diff(double(sub.idxA(:))) ~= 1);
    end

    if isfield(sub, 'resseq_af') && numel(sub.resseq_af) == N
        bA = bA | (diff(double(sub.resseq_af(:))) ~= 1);
    end

    if isfield(out_align, 'breakAfter_af_matched')
        src = logical(out_align.breakAfter_af_matched(:));
        bA = bA | source_breaks_for_selected_positions(src, pos, N, nPairs);
    end

    breakAfter_af(1:end-1) = bA;
end

% =========================================================================
function b = source_breaks_for_selected_positions(src, pos, N, nPairs)
% Copy source break flags only when the selected next row is the original
% next row. If the subset skips rows, selectedGap already handles that.

    b = false(N-1,1);

    if isempty(src)
        return;
    end

    maxUsable = min(numel(src), nPairs);
    p = pos(1:end-1);
    q = pos(2:end);

    valid = (p >= 1) & (p <= maxUsable) & (q == p + 1);
    if any(valid)
        b(valid) = src(p(valid));
    end
end