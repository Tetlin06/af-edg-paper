function [DistSq_all, Weight_all, info] = build_allatom_dist_weight_matrix(coords_true_all, cutoff, varargin)
% BUILD_ALLATOM_DIST_WEIGHT_MATRIX
% -------------------------------------------------------------------------
% Build an all-atom squared-distance matrix and binary contact/weight matrix.
%
% Inputs:
%   coords_true_all : [N x 3] true matched all-atom coordinates
%   cutoff          : distance cutoff in Angstroms, e.g. 6.0
%
% Outputs:
%   DistSq_all : [N x N] squared distances
%   Weight_all : [N x N] binary contact mask
%                Weight_all(i,j)=1 if distance <= cutoff, else 0
%                diagonal is always 0
%   info       : diagnostic struct
%
% Example:
%   [DistSq_all, Weight_all, info] = ...
%       build_allatom_dist_weight_matrix(allpairs.coords_true_all, 6.0);

    ip = inputParser;
    ip.FunctionName = mfilename;

    ip.addParameter('MinSeqSep', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
    ip.addParameter('AtomMeta', [], @(x) isempty(x) || istable(x));
    ip.addParameter('Verbose', true, @(x)islogical(x)||isnumeric(x));

    ip.parse(varargin{:});
    opt = ip.Results;

    verbose = logical(opt.Verbose);

    if size(coords_true_all,2) ~= 3
        error('coords_true_all must be N x 3.');
    end

    if ~isscalar(cutoff) || ~isnumeric(cutoff) || cutoff <= 0
        error('cutoff must be a positive scalar in Angstroms.');
    end

    N = size(coords_true_all,1);

    if N < 2
        error('Need at least 2 atoms.');
    end

    if any(~isfinite(coords_true_all(:)))
        error('coords_true_all contains NaN or Inf.');
    end

    % Pairwise squared distances
    X = coords_true_all;
    X2 = sum(X.^2, 2);
    DistSq_all = X2 + X2' - 2 * (X * X');
    DistSq_all(DistSq_all < 0) = 0;
    DistSq_all = 0.5 * (DistSq_all + DistSq_all');

    Dist_all = sqrt(DistSq_all);

    Weight_all = double(Dist_all <= cutoff);
    Weight_all(1:N+1:end) = 0;

    % Optional sequence separation filter.
    % This is atom-level, but uses residuePairIndex if AtomMeta is given.
    if opt.MinSeqSep > 0
        if isempty(opt.AtomMeta) || ~ismember('residuePairIndex', opt.AtomMeta.Properties.VariableNames)
            error('MinSeqSep requires AtomMeta with residuePairIndex column.');
        end

        resIdx = double(opt.AtomMeta.residuePairIndex(:));

        if numel(resIdx) ~= N
            error('AtomMeta height must match coords_true_all rows.');
        end

        sep = abs(resIdx - resIdx');
        Weight_all(sep < opt.MinSeqSep) = 0;
        Weight_all(1:N+1:end) = 0;
    end

    % Diagnostics
    numEdges = nnz(triu(Weight_all,1));
    possibleEdges = N*(N-1)/2;

    info = struct();
    info.N = N;
    info.cutoff = cutoff;
    info.MinSeqSep = opt.MinSeqSep;
    info.numEdges = numEdges;
    info.possibleEdges = possibleEdges;
    info.coverage_upper_triangle = numEdges / possibleEdges;
    info.coverage_full_matrix = nnz(Weight_all) / numel(Weight_all);

    positiveDistances = Dist_all(triu(Weight_all,1) > 0);
    if isempty(positiveDistances)
        info.minContactDist = NaN;
        info.maxContactDist = NaN;
        info.meanContactDist = NaN;
    else
        info.minContactDist = min(positiveDistances);
        info.maxContactDist = max(positiveDistances);
        info.meanContactDist = mean(positiveDistances);
    end

    if verbose
        fprintf('\n[ALL-ATOM DIST/WEIGHT MATRIX]\n');
        fprintf('  Atoms: %d\n', N);
        fprintf('  Cutoff: %.2f A\n', cutoff);
        fprintf('  MinSeqSep: %d\n', opt.MinSeqSep);
        fprintf('  Edges upper triangle: %d / %d\n', numEdges, possibleEdges);
        fprintf('  Coverage upper triangle: %.4f%%\n', 100 * info.coverage_upper_triangle);
        fprintf('  Coverage full matrix: %.4f%%\n', 100 * info.coverage_full_matrix);

        if ~isempty(positiveDistances)
            fprintf('  Contact distance range: %.3f..%.3f A | mean %.3f A\n', ...
                info.minContactDist, info.maxContactDist, info.meanContactDist);
        end
    end
end