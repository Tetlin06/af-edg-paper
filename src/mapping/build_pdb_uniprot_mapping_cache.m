function T_all = build_pdb_uniprot_mapping_cache(targetCSV, varargin)
% BUILD_PDB_UNIPROT_MAPPING_CACHE
% -------------------------------------------------------------------------
% Build a local residue-level PDB -> UniProt mapping CSV from official SIFTS
% XML files for every row in sheets/targets.csv.
%
% Example:
%   T = build_pdb_uniprot_mapping_cache('sheets/targets.csv');
%
% Expected targetCSV columns:
%   AFDB ID, PDB ID, Chain ID
%
% Output:
%   sheets/pdb_uniprot_residue_map.csv
%
% Notes:
% - This function does NOT hand-create residue mappings. It parses SIFTS.

    ip = inputParser;
    ip.addParameter('SiftsDir', fullfile('data','SIFTS'), @(s)ischar(s)||isstring(s));
    ip.addParameter('OutCSV', fullfile('sheets','pdb_uniprot_residue_map.csv'), @(s)ischar(s)||isstring(s));
    ip.addParameter('DownloadMissing', false, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('ForceDownload', false, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('SkipMissing', true, @(x)islogical(x)||isnumeric(x));
    ip.addParameter('AFDBID', "", @(s)ischar(s)||isstring(s));
    ip.addParameter('Verbose', true, @(x)islogical(x)||isnumeric(x));
    ip.parse(varargin{:});
    opt = ip.Results;

    targetCSV = string(targetCSV);
    siftsDir = string(opt.SiftsDir);
    outCSV = string(opt.OutCSV);
    onlyAFDB = string(opt.AFDBID);

    downloadMissing = logical(opt.DownloadMissing);
    forceDownload = logical(opt.ForceDownload);
    skipMissing = logical(opt.SkipMissing);
    verbose = logical(opt.Verbose);

    if forceDownload
        downloadMissing = true;
    end

    if ~isfile(targetCSV)
        error('targets.csv not found: %s', targetCSV);
    end

    if ~exist(siftsDir, 'dir')
        mkdir(siftsDir);
    end

    outDir = fileparts(outCSV);
    if strlength(outDir) > 0 && ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    targets = readtable(targetCSV, 'FileType','text', 'VariableNamingRule','preserve');

    required = {'AFDB ID','PDB ID','Chain ID'};
    for k = 1:numel(required)
        if ~any(strcmp(required{k}, targets.Properties.VariableNames))
            error('targets.csv is missing required column: %s', required{k});
        end
    end

    if strlength(onlyAFDB) > 0
        afdbCol = string(targets.("AFDB ID"));
        targets = targets(strcmpi(strtrim(afdbCol), strtrim(onlyAFDB)), :);

        if isempty(targets)
            error('AFDBID=%s was requested, but no matching row exists in %s.', onlyAFDB, targetCSV);
        end
    end

    T_all = table();
    skipped = strings(0,1);

    for r = 1:height(targets)
        afdbID  = scalar_text(targets.("AFDB ID"), r);
        pdbID   = upper(scalar_text(targets.("PDB ID"), r));
        chainID = scalar_text(targets.("Chain ID"), r);
        uniBase = base_uniprot(afdbID);

        if verbose
            fprintf('\n[MAPPING CACHE] %s | PDB=%s chain=%s | UniProt base=%s\n', ...
                afdbID, pdbID, chainID, uniBase);
        end

        siftsPath = find_sifts_xml(pdbID, siftsDir, verbose);

        if strlength(siftsPath) == 0 && downloadMissing
            siftsPath = download_sifts_xml(pdbID, siftsDir, verbose);
        end

        if strlength(siftsPath) == 0
            msg = sprintf('No local SIFTS XML found for PDB=%s. Skipping this target. Expected %s.xml or %s.xml.gz in %s.', ...
                pdbID, lower(pdbID), lower(pdbID), siftsDir);

            if skipMissing
                warning('%s', msg);
                skipped(end+1,1) = sprintf('%s | PDB=%s chain=%s', afdbID, pdbID, chainID); 
                continue;
            else
                error('%s', msg);
            end
        end

        T = parse_sifts_uniprot_mapping(siftsPath, pdbID, chainID, uniBase, ...
            'KeepObservedOnly', false, ...
            'Verbose', verbose);

        if isempty(T)
            warning('No real SIFTS mapping rows found for AFDB=%s | PDB=%s chain=%s. Skipping this target.', ...
                afdbID, pdbID, chainID);
            skipped(end+1,1) = sprintf('%s | PDB=%s chain=%s | no rows parsed', afdbID, pdbID, chainID); 
            continue;
        end

        AFDB_ID = repmat(afdbID, height(T), 1);
        T = addvars(T, AFDB_ID, 'Before', 1);

        T_all = [T_all; T]; 
    end

    if isempty(T_all)
        error(['No mapping rows were generated. This means no usable SIFTS files were found/parsed.\n' ...
               'Download the needed real SIFTS XML files into %s, or run this function with DownloadMissing=true.'], siftsDir);
    end

    writetable(T_all, outCSV);

    if verbose
        fprintf('\n[MAPPING CACHE DONE]\n');
        fprintf('  Real mapping rows written: %d\n', height(T_all));
        fprintf('  Output: %s\n', outCSV);

        if ~isempty(skipped)
            fprintf('  Skipped targets without usable SIFTS: %d\n', numel(skipped));
            for k = 1:numel(skipped)
                fprintf('    - %s\n', skipped(k));
            end
        end
    end
end

% =====================================================================
function siftsPath = find_sifts_xml(pdbID, siftsDir, verbose)
    pdbLower = lower(char(string(pdbID)));

    localXml = fullfile(siftsDir, sprintf('%s.xml', pdbLower));
    localGz  = fullfile(siftsDir, sprintf('%s.xml.gz', pdbLower));

    if isfile(localXml)
        siftsPath = string(localXml);
        if verbose
            fprintf('  Using local SIFTS XML: %s\n', siftsPath);
        end
        return;
    end

    if isfile(localGz)
        siftsPath = string(localGz);
        if verbose
            fprintf('  Using local SIFTS XML.GZ: %s\n', siftsPath);
        end
        return;
    end

    siftsPath = "";
end

% =====================================================================
function siftsPath = download_sifts_xml(pdbID, siftsDir, verbose)
    pdbLower = lower(char(string(pdbID)));

    localGz = fullfile(siftsDir, sprintf('%s.xml.gz', pdbLower));

    urls = strings(2,1);
    urls(1) = sprintf('https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/%s.xml.gz', pdbLower);

    splitCode = pdbLower(2:3);
    urls(2) = sprintf('https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/%s/%s.xml.gz', splitCode, pdbLower);

    lastErr = [];
    for u = 1:numel(urls)
        try
            if verbose
                fprintf('  Downloading real SIFTS XML.GZ: %s\n', urls(u));
            end
            websave(localGz, urls(u));
            siftsPath = string(localGz);
            return;
        catch ME
            lastErr = ME;
        end
    end

    if ~isempty(lastErr) && verbose
        fprintf(2, '  Could not download SIFTS for %s. Last error: %s\n', pdbID, lastErr.message);
    end

    siftsPath = "";
end

% =====================================================================
function x = scalar_text(col, r)
    if iscell(col)
        x = string(col{r});
    else
        x = string(col(r));
    end
end

% =====================================================================
function b = base_uniprot(acc)
    b = string(acc);
    b = regexprep(b, '^AF-', '');
    b = regexprep(b, '-F\d+$', '');
    b = regexprep(b, '-.*$', '');
end
