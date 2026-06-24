function m = evaluate_procrustes_metrics(predCoords, trueCoords)
% EVALUATE_PROCRUSTES_METRICS
% Best-reflection Procrustes CA-RMSD + GDT-TS + GDT-HA.
%
% Usage:
%   m = evaluate_procrustes_metrics(predCA, trueCA);
%
% Output:
%   m.rmsd
%   m.gdt_ts       thresholds [1 2 4 8]
%   m.gdt_ha       thresholds [0.5 1 2 4]
%   m.errors       per-residue CA errors after best alignment
%   m.aligned      aligned predicted coordinates

    if size(predCoords,2) ~= 3 || size(trueCoords,2) ~= 3
        error('Inputs must be N x 3 coordinate matrices.');
    end

    if size(predCoords,1) ~= size(trueCoords,1)
        error('predCoords and trueCoords must have the same number of rows.');
    end

    if any(~isfinite(predCoords(:))) || any(~isfinite(trueCoords(:)))
        error('Coordinates contain NaN or Inf.');
    end

    % Try normal rigid alignment.
    [~, Z_no] = procrustes(trueCoords, predCoords, ...
        'Scaling', false, ...
        'Reflection', false);

    err_no = sqrt(sum((Z_no - trueCoords).^2, 2));
    rmsd_no = sqrt(mean(err_no.^2));

    % Try reflected rigid alignment.
    [~, Z_ref] = procrustes(trueCoords, predCoords, ...
        'Scaling', false, ...
        'Reflection', true);

    err_ref = sqrt(sum((Z_ref - trueCoords).^2, 2));
    rmsd_ref = sqrt(mean(err_ref.^2));

    % Keep only the better result.
    if rmsd_no <= rmsd_ref
        Z = Z_no;
        err = err_no;
        rmsd = rmsd_no;
    else
        Z = Z_ref;
        err = err_ref;
        rmsd = rmsd_ref;
    end

    % GDT-TS: loose global-distance score.
    gdt_ts_thresholds = [1 2 4 8];
    gdt_ts = mean(arrayfun(@(t) mean(err <= t), gdt_ts_thresholds));

    % GDT-HA: stricter high-accuracy global-distance score.
    gdt_ha_thresholds = [0.5 1 2 4];
    gdt_ha = mean(arrayfun(@(t) mean(err <= t), gdt_ha_thresholds));

    m = struct();
    m.rmsd = rmsd;
    m.gdt_ts = gdt_ts;
    m.gdt_ha = gdt_ha;
    m.errors = err;
    m.aligned = Z;
    m.gdt_ts_thresholds = gdt_ts_thresholds;
    m.gdt_ha_thresholds = gdt_ha_thresholds;
end