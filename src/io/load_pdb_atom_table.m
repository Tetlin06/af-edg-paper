function T = load_pdb_atom_table(pdbPathOrID, varargin)
% LOAD_PDB_ATOM_TABLE
% -------------------------------------------------------------------------
% Shared PDB atom loader for the EDG / AlphaFold pipeline.
%
% This is the single source of truth for reading atoms from PDB files.
% It uses MATLAB pdbread, then converts the atom records into a clean table.
%
% This function is intentionally FULL-ATOM capable.
% CA-only logic should live in a wrapper, not here.
%
% Outputs:
%   T : table with one row per atom.
%
% Columns:
%   modelIndex   : model number used from PDBStruct.Model
%   recordType   : "ATOM" or "HETATM"
%   sourceIndex  : index inside the original pdbread Atom/HeterogenAtom array
%   atomSerial   : PDB atom serial number, if available
%   atomName     : atom name, e.g. "CA", "N", "C", "O", "CB"
%   altLoc       : alternate location indicator, "" if blank
%   resName      : residue name, e.g. "ALA"
%   chainID      : chain ID, "" if blank
%   resSeq       : numeric author residue number
%   iCode        : insertion code, "" if blank
%   resid        : residue ID string, e.g. "70", "70A"
%   chainResid   : chain + residue ID, e.g. "A:70", "A:70A"
%   x,y,z        : atom coordinates
%   occupancy    : PDB occupancy, if available
%   tempFactor   : PDB temperature factor
%                  For AlphaFold PDBs, this column stores pLDDT.
%   element      : element symbol, if available
%   charge       : atom charge, if available
%
% Requirements:
%   MATLAB Bioinformatics Toolbox, because this uses pdbread.
%
% Examples:
%   T = load_pdb_atom_table('data/AFDB/P52799.pdb');
%   T = load_pdb_atom_table('data/AFDB/P52799.pdb', 'ChainID', 'A');
%   T = load_pdb_atom_table('data/PDB/7K36.pdb', 'ChainID', 'I');
%
% Optional name-value arguments:
%   'ChainID'        : [] default, no chain filtering.
%                      Use 'A', 'I', etc. to keep one chain.
%   'ModelIndex'     : 1 default.
%   'IncludeHetero'  : false default. If true, also include HETATM records.
%   'AltLocPolicy'   : "all" default.
%                      "all"        keeps all alternate locations.
%                      "blank_or_A" keeps blank altLoc or altLoc A.
%                      "blank"      keeps only blank altLoc.
%                      "A"          keeps only altLoc A.
%   'Verbose'        : false default.

    % ---------------- options ----------------
    ip = inputParser;
    ip.FunctionName = mfilename;

    ip.addParameter('ChainID', [], @(x) isempty(x) || ischar(x) || isstring(x));
    ip.addParameter('ModelIndex', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1 && floor(x) == x);
    ip.addParameter('IncludeHetero', false, @(x)islogical(x) || isnumeric(x));
    ip.addParameter('AltLocPolicy', "all", @(x)ischar(x) || isstring(x));
    ip.addParameter('Verbose', false, @(x)islogical(x) || isnumeric(x));

    ip.parse(varargin{:});
    opt = ip.Results;

    modelIndex = double(opt.ModelIndex);
    includeHetero = logical(opt.IncludeHetero);
    altLocPolicy = lower(string(opt.AltLocPolicy));
    verbose = logical(opt.Verbose);

    validAltPolicies = ["all", "blank_or_a", "blank", "a"];
    if ~any(altLocPolicy == validAltPolicies)
        error('AltLocPolicy must be "all", "blank_or_A", "blank", or "A".');
    end

    % ---------------- toolbox check ----------------
    if exist('pdbread', 'file') == 0
        error(['pdbread not found. It requires MATLAB Bioinformatics Toolbox. ', ...
               'Install Bioinformatics Toolbox before using load_pdb_atom_table.']);
    end

    % ---------------- resolve path ----------------
    [pdbPath, cleanupObj] = resolve_pdb_path(pdbPathOrID); %#ok<NASGU>

    if verbose
        fprintf('[load_pdb_atom_table] Reading: %s\n', pdbPath);
    end

    % ---------------- read PDB ----------------
    PDB = pdbread(char(pdbPath));

    % ---------------- get selected model ----------------
    if isfield(PDB, 'Model') && ~isempty(PDB.Model)
        if modelIndex > numel(PDB.Model)
            error('Requested ModelIndex=%d, but PDB contains only %d model(s).', ...
                modelIndex, numel(PDB.Model));
        end

        modelStruct = PDB.Model(modelIndex);
    else
        if modelIndex ~= 1
            error('PDB has no Model array, so only ModelIndex=1 is valid.');
        end

        modelStruct = PDB;
    end

    atomRecords = get_struct_array_field(modelStruct, ...
        {'Atom', 'Atoms'}, struct([]));

    heteroRecords = struct([]);
    if includeHetero
        heteroRecords = get_struct_array_field(modelStruct, ...
            {'HeterogenAtom', 'HetAtom', 'HeteroAtom', 'HeterogenAtoms'}, struct([]));
    end

    % ---------------- storage arrays ----------------
    modelIndex_col = zeros(0,1);
    recordType_col = strings(0,1);
    sourceIndex_col = zeros(0,1);
    atomSerial_col = zeros(0,1);
    atomName_col = strings(0,1);
    altLoc_col = strings(0,1);
    resName_col = strings(0,1);
    chainID_col = strings(0,1);
    resSeq_col = zeros(0,1);
    iCode_col = strings(0,1);
    resid_col = strings(0,1);
    chainResid_col = strings(0,1);
    x_col = zeros(0,1);
    y_col = zeros(0,1);
    z_col = zeros(0,1);
    occupancy_col = zeros(0,1);
    tempFactor_col = zeros(0,1);
    element_col = strings(0,1);
    charge_col = strings(0,1);

    % ---------------- append ATOM and optionally HETATM ----------------
    append_records(atomRecords, "ATOM");

    if includeHetero
        append_records(heteroRecords, "HETATM");
    end

    % ---------------- assemble table ----------------
    T = table( ...
        modelIndex_col, ...
        recordType_col, ...
        sourceIndex_col, ...
        atomSerial_col, ...
        atomName_col, ...
        altLoc_col, ...
        resName_col, ...
        chainID_col, ...
        resSeq_col, ...
        iCode_col, ...
        resid_col, ...
        chainResid_col, ...
        x_col, ...
        y_col, ...
        z_col, ...
        occupancy_col, ...
        tempFactor_col, ...
        element_col, ...
        charge_col, ...
        'VariableNames', { ...
            'modelIndex', ...
            'recordType', ...
            'sourceIndex', ...
            'atomSerial', ...
            'atomName', ...
            'altLoc', ...
            'resName', ...
            'chainID', ...
            'resSeq', ...
            'iCode', ...
            'resid', ...
            'chainResid', ...
            'x', ...
            'y', ...
            'z', ...
            'occupancy', ...
            'tempFactor', ...
            'element', ...
            'charge' ...
        } ...
    );

    % ---------------- optional chain filter ----------------
    if ~isempty(opt.ChainID)
        wantedChain = clean_string(opt.ChainID);
        T = T(strcmpi(T.chainID, wantedChain), :);
    end

    % ---------------- optional altLoc filter ----------------
    switch altLocPolicy
        case "all"
            % keep everything

        case "blank_or_a"
            T = T((T.altLoc == "") | (upper(T.altLoc) == "A"), :);

        case "blank"
            T = T(T.altLoc == "", :);

        case "a"
            T = T(upper(T.altLoc) == "A", :);
    end

    if isempty(T)
        error('No atom records found after filtering. File=%s', pdbPath);
    end

    if verbose
        fprintf('[load_pdb_atom_table] Rows returned: %d\n', height(T));
        fprintf('[load_pdb_atom_table] ATOM rows: %d | HETATM rows: %d\n', ...
            nnz(T.recordType == "ATOM"), nnz(T.recordType == "HETATM"));

        caCount = nnz(upper(T.atomName) == "CA");
        fprintf('[load_pdb_atom_table] CA atoms: %d\n', caCount);
    end

    % =====================================================================
    % Nested helper: append atom-like records
    % =====================================================================
    function append_records(records, recordType)
        if isempty(records)
            return;
        end

        for a = 1:numel(records)
            rec = records(a);

            atomSerial = to_double(get_first_field(rec, ...
                {'AtomSerNo', 'atomSerNo', 'SerialNo', 'serialNo', 'serial'}, NaN));

            atomName = clean_string(get_first_field(rec, ...
                {'AtomName', 'atomName'}, ""));

            altLoc = clean_string(get_first_field(rec, ...
                {'altLoc', 'AltLoc', 'altloc'}, ""));

            resName = upper(clean_string(get_first_field(rec, ...
                {'resName', 'ResName', 'resname'}, "")));

            chainID = clean_string(get_first_field(rec, ...
                {'chainID', 'ChainID', 'chainId', 'chain'}, ""));

            resSeq = to_double(get_first_field(rec, ...
                {'resSeq', 'ResSeq', 'resseq'}, NaN));

            iCode = clean_string(get_first_field(rec, ...
                {'iCode', 'ICode', 'icode'}, ""));

            x = to_double(get_first_field(rec, {'X', 'x'}, NaN));
            y = to_double(get_first_field(rec, {'Y', 'y'}, NaN));
            z = to_double(get_first_field(rec, {'Z', 'z'}, NaN));

            occupancy = to_double(get_first_field(rec, ...
                {'occupancy', 'Occupancy'}, NaN));

            tempFactor = to_double(get_first_field(rec, ...
                {'tempFactor', 'TempFactor', 'BFactor', 'bFactor'}, NaN));

            element = upper(clean_string(get_first_field(rec, ...
                {'element', 'Element'}, "")));

            charge = clean_string(get_first_field(rec, ...
                {'charge', 'Charge'}, ""));

            resid = make_resid(resSeq, iCode);

            if chainID == ""
                chainResid = ":" + resid;
            else
                chainResid = chainID + ":" + resid;
            end

            modelIndex_col(end+1,1) = modelIndex; 
            recordType_col(end+1,1) = recordType; 
            sourceIndex_col(end+1,1) = a; 
            atomSerial_col(end+1,1) = atomSerial; 
            atomName_col(end+1,1) = atomName; 
            altLoc_col(end+1,1) = altLoc; 
            resName_col(end+1,1) = resName; 
            chainID_col(end+1,1) = chainID; 
            resSeq_col(end+1,1) = resSeq; 
            iCode_col(end+1,1) = iCode; 
            resid_col(end+1,1) = resid; 
            chainResid_col(end+1,1) = chainResid; 
            x_col(end+1,1) = x; 
            y_col(end+1,1) = y; 
            z_col(end+1,1) = z; 
            occupancy_col(end+1,1) = occupancy; 
            tempFactor_col(end+1,1) = tempFactor; 
            element_col(end+1,1) = element; 
            charge_col(end+1,1) = charge; 
        end
    end
end

% =========================================================================
function [pdbPath, cleanupObj] = resolve_pdb_path(pdbPathOrID)
% Resolve a path or AFDB/PDB-like ID to an actual PDB file.
% If the file is gzipped, unzip to a temporary folder and return that path.

    cleanupObj = [];

    raw = string(pdbPathOrID);
    raw = strtrim(raw);

    if strlength(raw) == 0
        error('Empty PDB path or ID.');
    end

    candidates = strings(0,1);
    candidates(end+1,1) = raw;

    if ~endsWith(lower(raw), ".pdb") && ~endsWith(lower(raw), ".pdb.gz")
        candidates(end+1,1) = fullfile('data', 'AFDB', raw + ".pdb");
        candidates(end+1,1) = fullfile('data', 'AFDB', raw + ".pdb.gz");
        candidates(end+1,1) = fullfile('data', 'PDB',  raw + ".pdb");
        candidates(end+1,1) = fullfile('data', 'PDB',  raw + ".pdb.gz");
    else
        if endsWith(lower(raw), ".pdb") && ~isfile(raw)
            candidates(end+1,1) = raw + ".gz";
        end
    end

    pdbPath = "";

    for k = 1:numel(candidates)
        if isfile(candidates(k))
            pdbPath = candidates(k);
            break;
        end
    end

    if strlength(pdbPath) == 0
        msg = sprintf('PDB file not found. Tried:\n');
        for k = 1:numel(candidates)
            msg = sprintf('%s  %s\n', msg, candidates(k));
        end
        error('%s', msg);
    end

    if endsWith(lower(pdbPath), ".gz")
        tmpDir = tempname;
        mkdir(tmpDir);

        files = gunzip(char(pdbPath), tmpDir);

        if isempty(files)
            error('Could not gunzip PDB file: %s', pdbPath);
        end

        pdbPath = string(files{1});
        cleanupObj = onCleanup(@() cleanup_dir(tmpDir));
    end
end

% =========================================================================
function cleanup_dir(d)
    if exist(d, 'dir')
        try
            rmdir(d, 's');
        catch
            % Do not fail the loader because temp cleanup failed.
        end
    end
end

% =========================================================================
function records = get_struct_array_field(s, names, defaultValue)
    records = defaultValue;

    for k = 1:numel(names)
        name = names{k};
        if isfield(s, name)
            records = s.(name);
            return;
        end
    end
end

% =========================================================================
function val = get_first_field(s, names, defaultValue)
    val = defaultValue;

    for k = 1:numel(names)
        name = names{k};
        if isfield(s, name)
            val = s.(name);
            return;
        end
    end
end

% =========================================================================
function s = clean_string(x)
    if isempty(x)
        s = "";
        return;
    end

    try
        s = string(x);
    catch
        s = "";
        return;
    end

    if isempty(s) || ismissing(s(1))
        s = "";
        return;
    end

    s = strtrim(s(1));

    if ismissing(s)
        s = "";
    end
end

% =========================================================================
function x = to_double(v)
    if isempty(v)
        x = NaN;
        return;
    end

    if isnumeric(v)
        x = double(v);
        if isempty(x)
            x = NaN;
        else
            x = x(1);
        end
        return;
    end

    x = str2double(string(v));

    if isempty(x) || ismissing(x)
        x = NaN;
    else
        x = double(x(1));
    end
end

% =========================================================================
function resid = make_resid(resSeq, iCode)
    if isnan(resSeq)
        resid = "";
        return;
    end

    if abs(resSeq - round(resSeq)) < 1e-9
        base = string(sprintf('%d', round(resSeq)));
    else
        base = string(sprintf('%g', resSeq));
    end

    iCode = clean_string(iCode);

    if iCode == ""
        resid = base;
    else
        resid = base + extractBetween(iCode, 1, 1);
    end
end