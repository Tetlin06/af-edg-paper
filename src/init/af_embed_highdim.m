function [P0, Q] = af_embed_highdim(Pinit, Rk, seed, jitter_scale)
% AF_EMBED_HIGHDIM
% -------------------------------------------------------------------------
% Embed AlphaFold 3D coordinates into a higher-dimensional space using
% P0 = Pinit * Q, where Q is 3×Rk with orthonormal rows.
%
% If jitter_scale = 0, this is exactly the original AF_rank embedding.
% If jitter_scale > 0, add small random jitter only in the dimensions
% orthogonal to the embedded AF 3D subspace.
%
% IMPORTANT:
% - This function does NOT recenter Pinit.
% - This function does NOT recenter P0.
% - main.m should already pass centered AlphaFold coordinates.
% - Only the jitter itself is centered, so jitter_scale = 0 preserves the
%   exact old AF_rank baseline.

    if nargin < 2 || isempty(Rk)
        Rk = 10;
    end

    if Rk < 3
        error('Rk must be >= 3. Got Rk=%d', Rk);
    end

    if nargin >= 3 && ~isempty(seed)
        rng(seed);
    end

    if nargin < 4 || isempty(jitter_scale)
        jitter_scale = 0;
    end

    if size(Pinit,2) ~= 3
        error('Pinit must be N×3. Got %dx%d', size(Pinit,1), size(Pinit,2));
    end

    % ---------------------------------------------------------------------
    % Original exact AF_rank embedding
    % ---------------------------------------------------------------------
    A = randn(Rk, 3);
    [U, ~] = qr(A, 0);   % U is Rk×3, U'*U = I3
    Q = U.';             % Q is 3×Rk, Q*Q' = I3

    P0 = Pinit * Q;

    % ---------------------------------------------------------------------
    % Optional jitter in the extra-dimensional complement
    % ---------------------------------------------------------------------
    if jitter_scale > 0
        % RMS radius of the already-centered AlphaFold structure.
        % This gives a protein-size-relative jitter scale.
        af_rms = sqrt(mean(sum(Pinit.^2, 2)));

        % Project random noise away from the embedded AF 3D subspace.
        % Q.' * Q projects onto the 3D AF subspace inside Rk dimensions.
        % eye(Rk) - Q.' * Q projects onto the unused extra dimensions.
        Proj_perp = eye(Rk) - Q.' * Q;

        Z = randn(size(P0));
        Z_perp = Z * Proj_perp;

        % Center only the jitter, not the AF embedding.
        Z_perp = Z_perp - mean(Z_perp, 1);

        % Scale jitter so its RMS radius equals jitter_scale * af_rms.
        z_rms = sqrt(mean(sum(Z_perp.^2, 2)));

        if z_rms > 0
            Z_perp = Z_perp * ((jitter_scale * af_rms) / z_rms);
        end

        P0 = P0 + Z_perp;
    end
end