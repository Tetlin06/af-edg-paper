function T = parse_sifts_uniprot_mapping(siftsXmlPath, pdbID, chainID, uniprotID, varargin)
% PARSE_SIFTS_UNIPROT_MAPPING
% -------------------------------------------------------------------------
% Parse a PDBe/SIFTS XML file into a residue-level PDB -> UniProt map.
%
% Example:
%   T = parse_sifts_uniprot_mapping('data/SIFTS/7k36.xml','7K36','I','Q5VSL9');
%
% Output table columns:
%   PDB_ID, Chain_ID, PDB_resNum, PDB_iCode, PDB_resNum_raw, PDB_resName,
%   UniProt_ID, UniProt_resNum, UniProt_resName,
%   PDBe_resNum, PDBe_resName, Observed
%
% Notes:
% - Observed=false when SIFTS says PDB dbResNum='null' or has a
%   residueDetail annotation 'Not_Observed'.
% - uniprotID may include an isoform suffix (e.g. Q5VSL9-2). Matching is
%   done against the base accession (Q5VSL9), because SIFTS often stores
%   canonical UniProt accessions.

    ip = inputParser;
    ip.addParameter('KeepObservedOnly', false, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('Verbose', true, @(x)islogical(x)||isnumeric(x));
    ip.parse(varargin{:});
    opt = ip.Results;

    siftsXmlPath = string(siftsXmlPath);
    pdbID = upper(string(pdbID));
    chainID = char(string(chainID));
    uniprotID = string(uniprotID);
    wantedUniProtBase = base_uniprot(uniprotID);

    if ~isfile(siftsXmlPath)
        error('SIFTS XML file not found: %s', siftsXmlPath);
    end

    % xmlread cannot read .gz directly. If needed, unzip to a temporary folder.
    cleanupObj = [];
    xmlPathForRead = char(siftsXmlPath);
    if endsWith(lower(siftsXmlPath), '.gz')
        tmpDir = tempname;
        mkdir(tmpDir);
        gunzip(char(siftsXmlPath), tmpDir);
        files = dir(fullfile(tmpDir, '*.xml'));
        if isempty(files)
            error('Could not unzip XML from: %s', siftsXmlPath);
        end
        xmlPathForRead = fullfile(tmpDir, files(1).name);
        cleanupObj = onCleanup(@() cleanup_dir(tmpDir)); 
    end

    doc = xmlread(xmlPathForRead);

    residueNodes = doc.getElementsByTagName('residue');
    if residueNodes.getLength == 0
        residueNodes = doc.getElementsByTagNameNS('*', 'residue');
    end
    nNodes = residueNodes.getLength;

    PDB_ID = strings(0,1);
    Chain_ID = strings(0,1);
    PDB_resNum = zeros(0,1);
    PDB_iCode = strings(0,1);
    PDB_resNum_raw = strings(0,1);
    PDB_resName = strings(0,1);
    UniProt_ID = strings(0,1);
    UniProt_resNum = zeros(0,1);
    UniProt_resName = strings(0,1);
    PDBe_resNum = zeros(0,1);
    PDBe_resName = strings(0,1);
    Observed = false(0,1);

    for k = 0:(nNodes-1)
        rNode = residueNodes.item(k);

        pdbeResNumRaw = get_attr(rNode, 'dbResNum');
        pdbeResNameRaw = get_attr(rNode, 'dbResName');

        pdbeResNum = str2double(pdbeResNumRaw);
        if isnan(pdbeResNum)
            pdbeResNum = NaN;
        end

        crossRefs = rNode.getElementsByTagName('crossRefDb');
        if crossRefs.getLength == 0
            crossRefs = rNode.getElementsByTagNameNS('*', 'crossRefDb');
        end

        pdbRef = [];
        uniRef = [];

        for j = 0:(crossRefs.getLength-1)
            cNode = crossRefs.item(j);
            src = string(get_attr(cNode, 'dbSource'));

            if strcmpi(src, 'PDB')
                thisPdb = upper(string(get_attr(cNode, 'dbAccessionId')));
                thisChain = char(string(get_attr(cNode, 'dbChainId')));

                if (strlength(pdbID) == 0 || strcmpi(thisPdb, pdbID)) && strcmpi(thisChain, chainID)
                    pdbRef = cNode;
                end

            elseif strcmpi(src, 'UniProt')
                thisUni = string(get_attr(cNode, 'dbAccessionId'));
                if strlength(wantedUniProtBase) == 0 || strcmpi(base_uniprot(thisUni), wantedUniProtBase)
                    uniRef = cNode;
                end
            end
        end

        if isempty(pdbRef) || isempty(uniRef)
            continue;
        end

        pdbRaw = string(get_attr(pdbRef, 'dbResNum'));
        pdbName = string(get_attr(pdbRef, 'dbResName'));
        uniAcc = string(get_attr(uniRef, 'dbAccessionId'));
        uniNum = str2double(get_attr(uniRef, 'dbResNum'));
        uniName = string(get_attr(uniRef, 'dbResName'));

        if isnan(uniNum)
            continue;
        end

        [pdbNum, pdbICode] = parse_pdb_resnum(pdbRaw);

        notObserved = strcmpi(pdbRaw, 'null') || isnan(pdbNum);

        detailNodes = rNode.getElementsByTagName('residueDetail');
        if detailNodes.getLength == 0
            detailNodes = rNode.getElementsByTagNameNS('*', 'residueDetail');
        end
        for d = 0:(detailNodes.getLength-1)
            txt = string(strtrim(char(detailNodes.item(d).getTextContent)));
            if strcmpi(txt, 'Not_Observed')
                notObserved = true;
            end
        end

        isObserved = ~notObserved;
        if logical(opt.KeepObservedOnly) && ~isObserved
            continue;
        end

        PDB_ID(end+1,1) = pdbID; 
        Chain_ID(end+1,1) = string(chainID); 
        PDB_resNum(end+1,1) = pdbNum; 
        PDB_iCode(end+1,1) = string(pdbICode); 
        PDB_resNum_raw(end+1,1) = pdbRaw; 
        PDB_resName(end+1,1) = pdbName; 
        UniProt_ID(end+1,1) = uniAcc; 
        UniProt_resNum(end+1,1) = uniNum; 
        UniProt_resName(end+1,1) = uniName; 
        PDBe_resNum(end+1,1) = pdbeResNum; 
        PDBe_resName(end+1,1) = string(pdbeResNameRaw); 
        Observed(end+1,1) = isObserved; 
    end

    T = table(PDB_ID, Chain_ID, PDB_resNum, PDB_iCode, PDB_resNum_raw, ...
        PDB_resName, UniProt_ID, UniProt_resNum, UniProt_resName, ...
        PDBe_resNum, PDBe_resName, Observed);

    if ~isempty(T)
        T = sortrows(T, {'Chain_ID','UniProt_resNum','PDB_resNum'});
    end

    if logical(opt.Verbose)
        fprintf('\n[SIFTS PARSE]\n');
        fprintf('  File: %s\n', siftsXmlPath);
        fprintf('  PDB=%s | Chain=%s | requested UniProt=%s\n', pdbID, string(chainID), uniprotID);
        fprintf('  Rows parsed: %d | Observed: %d | Not observed: %d\n', ...
            height(T), nnz(T.Observed), height(T)-nnz(T.Observed));
        if height(T) > 0
            obs = T(T.Observed,:);
            if ~isempty(obs)
                fprintf('  Observed UniProt range: %d..%d\n', min(obs.UniProt_resNum), max(obs.UniProt_resNum));
                fprintf('  Observed PDB range:     %s..%s\n', string(obs.PDB_resNum_raw(1)), string(obs.PDB_resNum_raw(end)));
            end
        end
    end
end

% =====================================================================
function s = get_attr(node, name)
    if node.hasAttribute(name)
        s = char(node.getAttribute(name));
    else
        s = '';
    end
end

% =====================================================================
function [num, icode] = parse_pdb_resnum(raw)
    raw = char(string(raw));
    raw = strtrim(raw);
    icode = '';
    num = NaN;

    if isempty(raw) || strcmpi(raw, 'null')
        return;
    end

    tok = regexp(raw, '^(-?\d+)([A-Za-z]?)$', 'tokens', 'once');
    if isempty(tok)
        num = str2double(raw);
        if isnan(num)
            return;
        end
    else
        num = str2double(tok{1});
        if numel(tok) >= 2
            icode = tok{2};
        end
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
function cleanup_dir(d)
    if exist(d, 'dir')
        rmdir(d, 's');
    end
end
