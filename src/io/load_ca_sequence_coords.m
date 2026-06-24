function [seq, coords, resid, resseqs, unknownCount, stats, CA] = load_ca_sequence_coords(pdbPathOrID, chainID, varargin)
% LOAD_CA_SEQUENCE_COORDS
% -------------------------------------------------------------------------
% Shared CA sequence/coordinate loader for alignment functions.
%
% This is a thin wrapper around load_ca_records.m.
%
% Purpose:
%   Alignment functions need the CA-level amino acid sequence, coordinates,
%   residue IDs, residue numbers, and unknown residue count.
%
% This avoids having separate internal CA parsers inside:
%   - align_true_vs_af_by_uniprot_mapping.m
%   - align_true_vs_af_by_full_sequence.m
%
% Outputs:
%   seq          : 1-letter amino acid sequence as a string
%   coords       : [N x 3] CA coordinates
%   resid        : [N x 1] residue IDs including insertion code, e.g. "70A"
%   resseqs      : [N x 1] numeric residue numbers
%   unknownCount : number of residues mapped to X
%   stats        : struct from load_ca_records, with added sequence fields
%   CA           : full CA table, optional extra output
%
% Example:
%   [seq, coords, resid, resseqs, unknownCount, stats] = ...
%       load_ca_sequence_coords('data/AFDB/P52799.pdb', 'A');

    if nargin < 2 || isempty(chainID)
        chainID = 'A';
    end

    [CA, stats] = load_ca_records(pdbPathOrID, chainID, varargin{:});

    if isempty(CA)
        error('No CA records returned by load_ca_records.');
    end

    % Build 1-letter sequence from CA.aa1
    if ~ismember('aa1', CA.Properties.VariableNames)
        error('CA table is missing aa1 column. Check load_ca_records.m.');
    end

    aa = string(CA.aa1(:));
    seq = string(strjoin(cellstr(aa), ''));

    % Coordinates
    if ismember('coords', CA.Properties.VariableNames)
        coords = CA.coords;
    else
        coords = [CA.x, CA.y, CA.z];
    end

    % Residue identifiers
    resid = CA.resid;
    resseqs = CA.resSeq;

    % Unknown residues
    unknownCount = sum(aa == "X");

    % Add sequence-specific stats
    stats.seqLength = strlength(seq);
    stats.nCA = height(CA);
    stats.nResidues = height(CA);
    stats.unknownCount = unknownCount;
    stats.firstResid = resid(1);
    stats.lastResid = resid(end);

    if isfield(stats, 'chainID')
        % keep existing
    else
        stats.chainID = string(chainID);
    end
end