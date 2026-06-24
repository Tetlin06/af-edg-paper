function out = edg_lddt_score(predCoords, trueCoords, cutoff, thresholds)
% EDG_LDDT_SCORE
% Local CA-distance preservation score, scaled 0..100.
%
% Usage:
%   out = edg_lddt_score(predCA, trueCA);
%
% Default:
%   local neighbors: true CA distance <= 15 A
%   error thresholds: [0.5 1 2 4] A
%
% Meaning:
%   100 = local distances are preserved almost perfectly
%     0 = local distances are badly wrong

    if nargin < 3 || isempty(cutoff)
        cutoff = 15;
    end

    if nargin < 4 || isempty(thresholds)
        thresholds = [0.5 1 2 4];
    end

    if size(predCoords,2) ~= 3 || size(trueCoords,2) ~= 3
        error('Inputs must be N x 3 coordinate matrices.');
    end

    if size(predCoords,1) ~= size(trueCoords,1)
        error('predCoords and trueCoords must have the same number of rows.');
    end

    if any(~isfinite(predCoords(:))) || any(~isfinite(trueCoords(:)))
        error('Coordinates contain NaN or Inf.');
    end

    N = size(trueCoords, 1);

    Dtrue = local_pairwise_dist(trueCoords);
    Dpred = local_pairwise_dist(predCoords);

    localMask = (Dtrue > 0) & (Dtrue <= cutoff);

    perResidue = NaN(N, 1);
    neighborCount = zeros(N, 1);

    for i = 1:N
        nbr = localMask(i, :);
        neighborCount(i) = sum(nbr);

        if neighborCount(i) == 0
            continue;
        end

        errors = abs(Dpred(i, nbr) - Dtrue(i, nbr));   % 1 x neighbors

        % Each pair gets partial credit:
        % error <= 0.5, 1, 2, 4 A gives 4/4, 3/4, 2/4, 1/4, or 0/4 credit.
        pass = errors(:) <= thresholds(:).';           % neighbors x thresholds
        perResidue(i) = 100 * mean(pass(:));
    end

    out = struct();
    out.global_score100 = mean(perResidue(isfinite(perResidue)));
    out.per_residue_score100 = perResidue;
    out.neighbor_count = neighborCount;
    out.cutoff = cutoff;
    out.thresholds = thresholds;
end

function D = local_pairwise_dist(X)
    X2 = sum(X.^2, 2);
    D2 = X2 + X2' - 2 * (X * X');
    D2(D2 < 0) = 0;
    D = sqrt(D2);
end