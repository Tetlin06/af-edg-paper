function out = build_aligned_all_atom_pairs(sub, truePdbPath, trueChain, afPdbPath, afChain, varargin)
% BUILD_ALIGNED_ALL_ATOM_PAIRS
% -------------------------------------------------------------------------
% Expand an existing residue-level TRUE <-> AF alignment into atom-level
% matched coordinate pairs.
%
% This function assumes residue correspondence has ALREADY been solved by:
%   align_true_vs_af_by_uniprot_mapping.m
%       -> subset_aligned_pairs.m
%
% It does NOT do sequence alignment or UniProt mapping itself.
%
% Input:
%   sub         : output from subset_aligned_pairs(out_align, DomainCut)
%   truePdbPath : path to true experimental PDB
%   trueChain   : true PDB chain ID
%   afPdbPath   : path to AF PDB
%   afChain     : AF chain ID, usually 'A'
%
% Output:
%   out.coords_true_all : [M x 3] matched TRUE atom coordinates
%   out.coords_af_all   : [M x 3] matched AF atom coordinates
%   out.atom_meta       : table describing each matched atom row
%   out.ca_atom_rows    : rows of coords_* corresponding to CA atoms
%   out.backbone_rows   : rows corresponding to N, CA, C, O, OXT
%
% Default matching rule:
%   AtomSelection = "safe"
%
%   For each matched residue pair:
%     - Always keep common backbone atoms.
%     - Keep side-chain atoms only when TRUE and AF represent the same
%       amino acid class.
%
% Why:
%   Matching all side-chain atom names across different residue identities
%   is not biologically meaningful. Backbone atoms are still comparable.
%
% Name-value options:
%   'AtomSelection':
%       "safe"                 default
%       "backbone_only"         only N, CA, C, O, OXT if common
%       "common_atoms"          all common atom names, even if residue differs
%       "strict_same_residue"   skip residue if amino acid class differs
%
%   'BackboneAtoms'    : default ["N","CA","C","O","OXT"]
%   'RemoveHydrogens'  : true default
%   'IncludeHetero'    : false default
%   'AltLocPolicy'     : "blank_or_A" default
%   'KeepLoadedTables' : false default
%   'Verbose'          : true default

    % ---------------- options ----------------
    ip = inputParser;
    ip.FunctionName = mfilename;

    ip.addParameter('AtomSelection', "safe", @(x)ischar(x)||isstring(x));
    ip.addParameter('BackboneAtoms', ["N","CA","C","O","OXT"], @(x)ischar(x)||isstring(x)||iscellstr(x));
    ip.addParameter('RemoveHydrogens', true, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('IncludeHetero', false, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('AltLocPolicy', "blank_or_A", @(x)ischar(x)||isstring(x));
    ip.addParameter('KeepLoadedTables', false, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('Verbose', true, @(x)islogical(x)||isnumeric(x));

    ip.parse(varargin{:});
    opt = ip.Results;

    atomSelection = lower(string(opt.AtomSelection));
    validModes = ["safe","backbone_only","common_atoms","strict_same_residue"];
    if ~any(atomSelection == validModes)
        error('AtomSelection must be "safe", "backbone_only", "common_atoms", or "strict_same_residue".');
    end

    backboneAtoms = upper(string(opt.BackboneAtoms));
    backboneAtoms = backboneAtoms(:).';

    verbose = logical(opt.Verbose);

    % ---------------- required sub fields ----------------
    requiredFields = {'resid_true','resid_af','resseq_true','resseq_af'};
    for k = 1:numel(requiredFields)
        if ~isfield(sub, requiredFields{k})
            error('sub is missing required field: %s', requiredFields{k});
        end
    end

    nResidues = numel(sub.resid_true);

    if numel(sub.resid_af) ~= nResidues
        error('sub.resid_true and sub.resid_af must have the same length.');
    end

    if isfield(sub, 'pos')
        alignedPairPos = sub.pos(:);
    else
        alignedPairPos = (1:nResidues).';
    end

    % ---------------- load all atoms through shared loader ----------------
    [ATtrue, statsTrue] = load_all_atom_records(truePdbPath, trueChain, ...
        'IncludeHetero', logical(opt.IncludeHetero), ...
        'RemoveHydrogens', logical(opt.RemoveHydrogens), ...
        'AltLocPolicy', opt.AltLocPolicy, ...
        'Verbose', false);

    [ATaf, statsAF] = load_all_atom_records(afPdbPath, afChain, ...
        'IncludeHetero', logical(opt.IncludeHetero), ...
        'RemoveHydrogens', logical(opt.RemoveHydrogens), ...
        'AltLocPolicy', opt.AltLocPolicy, ...
        'Verbose', false);

    % ---------------- storage for atom-level rows ----------------
    atomRow_col = zeros(0,1);
    residuePairIndex_col = zeros(0,1);
    alignedPairPos_col = zeros(0,1);

    trueResid_col = strings(0,1);
    afResid_col = strings(0,1);
    trueChainResid_col = strings(0,1);
    afChainResid_col = strings(0,1);

    trueResSeq_col = zeros(0,1);
    afResSeq_col = zeros(0,1);

    trueResName_col = strings(0,1);
    afResName_col = strings(0,1);
    trueAA_col = strings(0,1);
    afAA_col = strings(0,1);
    sameAA_col = false(0,1);

    atomName_col = strings(0,1);
    atomClass_col = strings(0,1);

    trueAtomIndex_col = zeros(0,1);
    afAtomIndex_col = zeros(0,1);
    trueAtomKey_col = strings(0,1);
    afAtomKey_col = strings(0,1);

    xT_col = zeros(0,1);
    yT_col = zeros(0,1);
    zT_col = zeros(0,1);

    xA_col = zeros(0,1);
    yA_col = zeros(0,1);
    zA_col = zeros(0,1);

    isCA_col = false(0,1);
    isBackbone_col = false(0,1);

    % ---------------- residue summary storage ----------------
    sum_residuePairIndex = zeros(nResidues,1);
    sum_alignedPairPos = zeros(nResidues,1);
    sum_trueResid = strings(nResidues,1);
    sum_afResid = strings(nResidues,1);
    sum_trueResName = strings(nResidues,1);
    sum_afResName = strings(nResidues,1);
    sum_trueAA = strings(nResidues,1);
    sum_afAA = strings(nResidues,1);
    sum_sameAA = false(nResidues,1);
    sum_nTrueAtoms = zeros(nResidues,1);
    sum_nAFAtoms = zeros(nResidues,1);
    sum_nCommonAtoms = zeros(nResidues,1);
    sum_nSelectedAtoms = zeros(nResidues,1);
    sum_status = strings(nResidues,1);

    % ---------------- main expansion loop ----------------
    atomCounter = 0;

    trueChainNorm = normalize_chain(trueChain);
    afChainNorm = normalize_chain(afChain);

    for r = 1:nResidues
        trueResid = string(sub.resid_true(r));
        afResid   = string(sub.resid_af(r));

        trueChainResid = trueChainNorm + ":" + trueResid;
        afChainResid   = afChainNorm + ":" + afResid;

        TtrueRes = ATtrue(ATtrue.chainResid == trueChainResid, :);
        TafRes   = ATaf(ATaf.chainResid == afChainResid, :);

        sum_residuePairIndex(r) = r;
        sum_alignedPairPos(r) = alignedPairPos(r);
        sum_trueResid(r) = trueResid;
        sum_afResid(r) = afResid;

        if isempty(TtrueRes)
            sum_status(r) = "missing_true_residue_atoms";
            continue;
        end

        if isempty(TafRes)
            sum_status(r) = "missing_af_residue_atoms";
            continue;
        end

        trueResName = TtrueRes.resName(1);
        afResName   = TafRes.resName(1);

        trueAA = string(aa3_to_1(trueResName));
        afAA   = string(aa3_to_1(afResName));
        sameAA = trueAA == afAA && trueAA ~= "X" && afAA ~= "X";

        sum_trueResName(r) = trueResName;
        sum_afResName(r) = afResName;
        sum_trueAA(r) = trueAA;
        sum_afAA(r) = afAA;
        sum_sameAA(r) = sameAA;

        trueNames = unique(upper(TtrueRes.atomName), 'stable');
        afNames   = unique(upper(TafRes.atomName), 'stable');

        commonNames = trueNames(ismember(trueNames, afNames));
        commonNames = commonNames(:).';

        isCommonBackbone = ismember(commonNames, backboneAtoms);

        switch atomSelection
            case "safe"
                if sameAA
                    selectedNames = commonNames;
                else
                    selectedNames = commonNames(isCommonBackbone);
                end

            case "backbone_only"
                selectedNames = commonNames(isCommonBackbone);

            case "common_atoms"
                selectedNames = commonNames;

            case "strict_same_residue"
                if sameAA
                    selectedNames = commonNames;
                else
                    selectedNames = strings(1,0);
                end
        end

        sum_nTrueAtoms(r) = height(TtrueRes);
        sum_nAFAtoms(r) = height(TafRes);
        sum_nCommonAtoms(r) = numel(commonNames);
        sum_nSelectedAtoms(r) = numel(selectedNames);

        if isempty(selectedNames)
            if sameAA
                sum_status(r) = "no_common_atoms";
            else
                sum_status(r) = "residue_mismatch_no_atoms_kept";
            end
            continue;
        end

        if sameAA
            sum_status(r) = "ok_same_aa";
        else
            sum_status(r) = "ok_backbone_only_due_to_mismatch";
        end

        for a = 1:numel(selectedNames)
            atomName = selectedNames(a);

            it = find(upper(TtrueRes.atomName) == atomName, 1, 'first');
            ia = find(upper(TafRes.atomName) == atomName, 1, 'first');

            if isempty(it) || isempty(ia)
                continue;
            end

            trow = TtrueRes(it, :);
            arow = TafRes(ia, :);

            atomCounter = atomCounter + 1;

            atomRow_col(end+1,1) = atomCounter; 
            residuePairIndex_col(end+1,1) = r; 
            alignedPairPos_col(end+1,1) = alignedPairPos(r); 

            trueResid_col(end+1,1) = trueResid; 
            afResid_col(end+1,1) = afResid; 
            trueChainResid_col(end+1,1) = trueChainResid; 
            afChainResid_col(end+1,1) = afChainResid; 

            trueResSeq_col(end+1,1) = sub.resseq_true(r); 
            afResSeq_col(end+1,1) = sub.resseq_af(r); 

            trueResName_col(end+1,1) = trueResName; 
            afResName_col(end+1,1) = afResName; 
            trueAA_col(end+1,1) = trueAA; 
            afAA_col(end+1,1) = afAA; 
            sameAA_col(end+1,1) = sameAA; 

            atomName_col(end+1,1) = atomName; 

            if ismember(atomName, backboneAtoms)
                atomClass_col(end+1,1) = "backbone"; 
                isBackbone_col(end+1,1) = true; 
            else
                atomClass_col(end+1,1) = "sidechain"; 
                isBackbone_col(end+1,1) = false; 
            end

            trueAtomIndex_col(end+1,1) = trow.atomIndex(1); 
            afAtomIndex_col(end+1,1) = arow.atomIndex(1); 
            trueAtomKey_col(end+1,1) = trow.atomKey(1); 
            afAtomKey_col(end+1,1) = arow.atomKey(1); 

            xT_col(end+1,1) = trow.x(1); 
            yT_col(end+1,1) = trow.y(1); 
            zT_col(end+1,1) = trow.z(1); 

            xA_col(end+1,1) = arow.x(1); 
            yA_col(end+1,1) = arow.y(1); 
            zA_col(end+1,1) = arow.z(1); 

            isCA_col(end+1,1) = (atomName == "CA"); 
        end
    end

    if atomCounter == 0
        error('No matched atom pairs were produced. Check chains, residue IDs, and AtomSelection.');
    end

    % ---------------- package atom metadata ----------------
    atom_meta = table( ...
        atomRow_col, ...
        residuePairIndex_col, ...
        alignedPairPos_col, ...
        trueResid_col, ...
        afResid_col, ...
        trueChainResid_col, ...
        afChainResid_col, ...
        trueResSeq_col, ...
        afResSeq_col, ...
        trueResName_col, ...
        afResName_col, ...
        trueAA_col, ...
        afAA_col, ...
        sameAA_col, ...
        atomName_col, ...
        atomClass_col, ...
        trueAtomIndex_col, ...
        afAtomIndex_col, ...
        trueAtomKey_col, ...
        afAtomKey_col, ...
        xT_col, yT_col, zT_col, ...
        xA_col, yA_col, zA_col, ...
        isCA_col, ...
        isBackbone_col, ...
        'VariableNames', { ...
            'atomRow', ...
            'residuePairIndex', ...
            'alignedPairPos', ...
            'trueResid', ...
            'afResid', ...
            'trueChainResid', ...
            'afChainResid', ...
            'trueResSeq', ...
            'afResSeq', ...
            'trueResName', ...
            'afResName', ...
            'trueAA', ...
            'afAA', ...
            'sameAA', ...
            'atomName', ...
            'atomClass', ...
            'trueAtomIndex', ...
            'afAtomIndex', ...
            'trueAtomKey', ...
            'afAtomKey', ...
            'xTrue', 'yTrue', 'zTrue', ...
            'xAF', 'yAF', 'zAF', ...
            'isCA', ...
            'isBackbone' ...
        } ...
    );

    coords_true_all = [atom_meta.xTrue, atom_meta.yTrue, atom_meta.zTrue];
    coords_af_all   = [atom_meta.xAF,   atom_meta.yAF,   atom_meta.zAF];

    ca_atom_rows = find(atom_meta.isCA);
    backbone_rows = find(atom_meta.isBackbone);
    sidechain_rows = find(~atom_meta.isBackbone);

    % ---------------- package residue summary ----------------
    residue_summary = table( ...
        sum_residuePairIndex, ...
        sum_alignedPairPos, ...
        sum_trueResid, ...
        sum_afResid, ...
        sum_trueResName, ...
        sum_afResName, ...
        sum_trueAA, ...
        sum_afAA, ...
        sum_sameAA, ...
        sum_nTrueAtoms, ...
        sum_nAFAtoms, ...
        sum_nCommonAtoms, ...
        sum_nSelectedAtoms, ...
        sum_status, ...
        'VariableNames', { ...
            'residuePairIndex', ...
            'alignedPairPos', ...
            'trueResid', ...
            'afResid', ...
            'trueResName', ...
            'afResName', ...
            'trueAA', ...
            'afAA', ...
            'sameAA', ...
            'nTrueAtoms', ...
            'nAFAtoms', ...
            'nCommonAtoms', ...
            'nSelectedAtoms', ...
            'status' ...
        } ...
    );

    % ---------------- output struct ----------------
    out = struct();
    out.truePdbPath = string(truePdbPath);
    out.afPdbPath = string(afPdbPath);
    out.trueChain = string(trueChain);
    out.afChain = string(afChain);

    out.atomSelection = atomSelection;
    out.backboneAtoms = backboneAtoms;

    out.coords_true_all = coords_true_all;
    out.coords_af_all = coords_af_all;
    out.atom_meta = atom_meta;
    out.residue_summary = residue_summary;

    out.ca_atom_rows = ca_atom_rows;
    out.backbone_rows = backbone_rows;
    out.sidechain_rows = sidechain_rows;

    out.true_atom_indices = atom_meta.trueAtomIndex;
    out.af_atom_indices = atom_meta.afAtomIndex;

    out.n_atoms = height(atom_meta);
    out.n_residue_pairs_requested = nResidues;
    out.n_residue_pairs_with_atoms = nnz(residue_summary.nSelectedAtoms > 0);
    out.n_ca_atoms = numel(ca_atom_rows);
    out.n_backbone_atoms = numel(backbone_rows);
    out.n_sidechain_atoms = numel(sidechain_rows);

    out.n_residue_mismatches = nnz(~residue_summary.sameAA);
    out.n_residue_mismatches_with_atoms_kept = nnz(~residue_summary.sameAA & residue_summary.nSelectedAtoms > 0);
    out.n_residue_pairs_without_atoms = nnz(residue_summary.nSelectedAtoms == 0);

    out.stats_true_all = statsTrue;
    out.stats_af_all = statsAF;

    if logical(opt.KeepLoadedTables)
        out.ATtrue = ATtrue;
        out.ATaf = ATaf;
    end

    % ---------------- reporting ----------------
    if verbose
        fprintf('\n[ALL-ATOM PAIRS]\n');
        fprintf('  Atom selection mode: %s\n', atomSelection);
        fprintf('  Residue pairs requested: %d\n', out.n_residue_pairs_requested);
        fprintf('  Residue pairs with >=1 atom pair: %d\n', out.n_residue_pairs_with_atoms);
        fprintf('  Matched atom pairs: %d\n', out.n_atoms);
        fprintf('  CA atom pairs: %d\n', out.n_ca_atoms);
        fprintf('  Backbone atom pairs: %d\n', out.n_backbone_atoms);
        fprintf('  Sidechain atom pairs: %d\n', out.n_sidechain_atoms);
        fprintf('  Residue amino-acid mismatches: %d\n', out.n_residue_mismatches);
        fprintf('  Mismatched residues with atoms kept: %d\n', out.n_residue_mismatches_with_atoms_kept);

        if out.n_ca_atoms ~= nResidues
            fprintf(2, '  [WARN] CA atom pairs (%d) != residue pairs requested (%d).\n', ...
                out.n_ca_atoms, nResidues);
        end
    end
end

% =========================================================================
function chainID = normalize_chain(chainID)
    chainID = string(chainID);
    chainID = strtrim(chainID);

    if strlength(chainID) == 0
        chainID = "";
    else
        chainID = extractBetween(chainID, 1, 1);
    end
end

% =========================================================================
function aa = aa3_to_1(resName)
% Local copy so this all-atom matcher does not depend on private helpers.

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
        case 'SEP'
            aa = 'S';
        case 'TPO'
            aa = 'T';
        case 'PTR'
            aa = 'Y';
        case {'HSD','HSE','HSP'}
            aa = 'H';
        case {'CSO','CSD','CSS','CME','CSE'}
            aa = 'C';
        otherwise
            aa = 'X';
    end
end