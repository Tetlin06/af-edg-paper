function out = build_aligned_all_atom_pairs(alignment, truePdbPath, trueChain, afPdbPath, afChain)
% BUILD_ALIGNED_ALL_ATOM_PAIRS
% Expand the full residue-level TRUE <-> AlphaFold alignment into matched
% all-atom coordinate pairs.
%
% Residue matching is assumed to have already been produced by
% align_true_vs_af_by_uniprot_mapping. This function uses the full matched
% alignment and applies the safe atom-matching rule:
%   - always keep common backbone atoms
%   - keep side-chain atoms only when residue identities agree

    required = { ...
        'n_aligned_pairs', ...
        'resid_true_matched', ...
        'resid_af_matched', ...
        'resseqs_true_matched', ...
        'resseqs_af_matched'};

    for k = 1:numel(required)
        if ~isfield(alignment, required{k})
            error('alignment is missing required field: %s', required{k});
        end
    end

    nResidues = alignment.n_aligned_pairs;
    if nResidues <= 0
        error('Alignment contains no matched residue pairs.');
    end

    residTrue = string(alignment.resid_true_matched(:));
    residAF   = string(alignment.resid_af_matched(:));
    resseqTrue = alignment.resseqs_true_matched(:);
    resseqAF   = alignment.resseqs_af_matched(:);

    if numel(residTrue) ~= nResidues || numel(residAF) ~= nResidues
        error('Alignment residue fields do not match n_aligned_pairs.');
    end

    backboneAtoms = ["N", "CA", "C", "O", "OXT"];

    [ATtrue, statsTrue] = load_all_atom_records(truePdbPath, trueChain, 'Verbose', false);
    [ATaf, statsAF] = load_all_atom_records(afPdbPath, afChain, 'Verbose', false);

    trueChainNorm = normalize_chain(trueChain);
    afChainNorm = normalize_chain(afChain);

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

    atomCounter = 0;

    for r = 1:nResidues
        tr = residTrue(r);
        ar = residAF(r);

        trueChainResid = trueChainNorm + ":" + tr;
        afChainResid   = afChainNorm + ":" + ar;

        TtrueRes = ATtrue(ATtrue.chainResid == trueChainResid, :);
        TafRes   = ATaf(ATaf.chainResid == afChainResid, :);

        sum_residuePairIndex(r) = r;
        sum_alignedPairPos(r) = r;
        sum_trueResid(r) = tr;
        sum_afResid(r) = ar;

        if isempty(TtrueRes)
            sum_status(r) = "missing_true_residue_atoms";
            continue;
        end

        if isempty(TafRes)
            sum_status(r) = "missing_af_residue_atoms";
            continue;
        end

        trueResName = string(TtrueRes.resName(1));
        afResName   = string(TafRes.resName(1));

        trueAA = string(aa3_to_1(trueResName));
        afAA   = string(aa3_to_1(afResName));
        sameAA = trueAA == afAA && trueAA ~= "X" && afAA ~= "X";

        sum_trueResName(r) = trueResName;
        sum_afResName(r) = afResName;
        sum_trueAA(r) = trueAA;
        sum_afAA(r) = afAA;
        sum_sameAA(r) = sameAA;

        trueNames = unique(upper(string(TtrueRes.atomName)), 'stable');
        afNames   = unique(upper(string(TafRes.atomName)), 'stable');
        commonNames = trueNames(ismember(trueNames, afNames));
        commonNames = commonNames(:).';

        if sameAA
            selectedNames = commonNames;
        else
            selectedNames = commonNames(ismember(commonNames, backboneAtoms));
        end

        sum_nTrueAtoms(r) = height(TtrueRes);
        sum_nAFAtoms(r) = height(TafRes);
        sum_nCommonAtoms(r) = numel(commonNames);
        sum_nSelectedAtoms(r) = numel(selectedNames);

        if isempty(selectedNames)
            if sameAA
                sum_status(r) = "no_common_atoms";
            else
                sum_status(r) = "residue_mismatch_no_backbone_atoms";
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

            it = find(upper(string(TtrueRes.atomName)) == atomName, 1, 'first');
            ia = find(upper(string(TafRes.atomName)) == atomName, 1, 'first');

            if isempty(it) || isempty(ia)
                continue;
            end

            trow = TtrueRes(it, :);
            arow = TafRes(ia, :);

            atomCounter = atomCounter + 1;

            atomRow_col(end+1,1) = atomCounter; %#ok<AGROW>
            residuePairIndex_col(end+1,1) = r; %#ok<AGROW>
            alignedPairPos_col(end+1,1) = r; %#ok<AGROW>

            trueResid_col(end+1,1) = tr; %#ok<AGROW>
            afResid_col(end+1,1) = ar; %#ok<AGROW>
            trueChainResid_col(end+1,1) = trueChainResid; %#ok<AGROW>
            afChainResid_col(end+1,1) = afChainResid; %#ok<AGROW>

            trueResSeq_col(end+1,1) = resseqTrue(r); %#ok<AGROW>
            afResSeq_col(end+1,1) = resseqAF(r); %#ok<AGROW>

            trueResName_col(end+1,1) = trueResName; %#ok<AGROW>
            afResName_col(end+1,1) = afResName; %#ok<AGROW>
            trueAA_col(end+1,1) = trueAA; %#ok<AGROW>
            afAA_col(end+1,1) = afAA; %#ok<AGROW>
            sameAA_col(end+1,1) = sameAA; %#ok<AGROW>

            atomName_col(end+1,1) = atomName; %#ok<AGROW>

            isBackbone = ismember(atomName, backboneAtoms);
            isCA = atomName == "CA";

            if isCA
                atomClass = "CA";
            elseif isBackbone
                atomClass = "backbone";
            else
                atomClass = "sidechain";
            end

            atomClass_col(end+1,1) = atomClass; %#ok<AGROW>

            trueAtomIndex_col(end+1,1) = trow.atomIndex; %#ok<AGROW>
            afAtomIndex_col(end+1,1) = arow.atomIndex; %#ok<AGROW>
            trueAtomKey_col(end+1,1) = string(trow.atomKey); %#ok<AGROW>
            afAtomKey_col(end+1,1) = string(arow.atomKey); %#ok<AGROW>

            xT_col(end+1,1) = trow.x; %#ok<AGROW>
            yT_col(end+1,1) = trow.y; %#ok<AGROW>
            zT_col(end+1,1) = trow.z; %#ok<AGROW>

            xA_col(end+1,1) = arow.x; %#ok<AGROW>
            yA_col(end+1,1) = arow.y; %#ok<AGROW>
            zA_col(end+1,1) = arow.z; %#ok<AGROW>

            isCA_col(end+1,1) = isCA; %#ok<AGROW>
            isBackbone_col(end+1,1) = isBackbone; %#ok<AGROW>
        end
    end

    if atomCounter == 0
        error('No matched all-atom pairs were produced.');
    end

    atom_meta = table( ...
        atomRow_col, residuePairIndex_col, alignedPairPos_col, ...
        trueResid_col, afResid_col, trueChainResid_col, afChainResid_col, ...
        trueResSeq_col, afResSeq_col, ...
        trueResName_col, afResName_col, trueAA_col, afAA_col, sameAA_col, ...
        atomName_col, atomClass_col, ...
        trueAtomIndex_col, afAtomIndex_col, trueAtomKey_col, afAtomKey_col, ...
        xT_col, yT_col, zT_col, xA_col, yA_col, zA_col, ...
        isCA_col, isBackbone_col, ...
        'VariableNames', { ...
            'atomRow', 'residuePairIndex', 'alignedPairPos', ...
            'trueResid', 'afResid', 'trueChainResid', 'afChainResid', ...
            'trueResSeq', 'afResSeq', ...
            'trueResName', 'afResName', 'trueAA', 'afAA', 'sameAA', ...
            'atomName', 'atomClass', ...
            'trueAtomIndex', 'afAtomIndex', 'trueAtomKey', 'afAtomKey', ...
            'x_true', 'y_true', 'z_true', 'x_af', 'y_af', 'z_af', ...
            'isCA', 'isBackbone'});

    coords_true_all = [atom_meta.x_true, atom_meta.y_true, atom_meta.z_true];
    coords_af_all = [atom_meta.x_af, atom_meta.y_af, atom_meta.z_af];

    ca_atom_rows = find(atom_meta.isCA);
    backbone_rows = find(atom_meta.isBackbone);

    residue_summary = table( ...
        sum_residuePairIndex, sum_alignedPairPos, ...
        sum_trueResid, sum_afResid, sum_trueResName, sum_afResName, ...
        sum_trueAA, sum_afAA, sum_sameAA, ...
        sum_nTrueAtoms, sum_nAFAtoms, sum_nCommonAtoms, sum_nSelectedAtoms, sum_status, ...
        'VariableNames', { ...
            'residuePairIndex', 'alignedPairPos', ...
            'trueResid', 'afResid', 'trueResName', 'afResName', ...
            'trueAA', 'afAA', 'sameAA', ...
            'nTrueAtoms', 'nAFAtoms', 'nCommonAtoms', 'nSelectedAtoms', 'status'});

    out = struct();
    out.coords_true_all = coords_true_all;
    out.coords_af_all = coords_af_all;
    out.atom_meta = atom_meta;
    out.ca_atom_rows = ca_atom_rows;
    out.backbone_rows = backbone_rows;
    out.residue_summary = residue_summary;
    out.n_residue_pairs = nResidues;
    out.n_atom_pairs = size(coords_true_all, 1);
    out.n_ca_pairs = numel(ca_atom_rows);
    out.n_backbone_pairs = numel(backbone_rows);
    out.truePdbPath = string(truePdbPath);
    out.afPdbPath = string(afPdbPath);
    out.trueChain = string(trueChain);
    out.afChain = string(afChain);
    out.statsTrue = statsTrue;
    out.statsAF = statsAF;
end

% =========================================================================
function c = normalize_chain(chainID)
    c = string(chainID);
    c = strtrim(c);
    if c == ""
        c = " ";
    end
end

% =========================================================================
function aa = aa3_to_1(resName)
    r = upper(strtrim(char(string(resName))));
    switch r
        case 'ALA', aa = 'A';
        case 'ARG', aa = 'R';
        case 'ASN', aa = 'N';
        case 'ASP', aa = 'D';
        case 'CYS', aa = 'C';
        case 'GLN', aa = 'Q';
        case 'GLU', aa = 'E';
        case 'GLY', aa = 'G';
        case 'HIS', aa = 'H';
        case 'ILE', aa = 'I';
        case 'LEU', aa = 'L';
        case 'LYS', aa = 'K';
        case 'MET', aa = 'M';
        case 'PHE', aa = 'F';
        case 'PRO', aa = 'P';
        case 'SER', aa = 'S';
        case 'THR', aa = 'T';
        case 'TRP', aa = 'W';
        case 'TYR', aa = 'Y';
        case 'VAL', aa = 'V';
        case 'SEC', aa = 'U';
        case 'PYL', aa = 'O';
        otherwise, aa = 'X';
    end
end
