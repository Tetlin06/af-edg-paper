function P0 = make_paper_initialization(cfg, Pinit3, DistSq, Weight)
%MAKE_PAPER_INITIALIZATION Initializations used by the paper scripts only.
    method = string(cfg.init_method);

    switch method
        case "Random_rand"
            rng(cfg.init_seed, 'twister');
            P0 = rand(size(DistSq, 1), cfg.opts.rank);

        case "Random_randn"
            rng(cfg.init_seed, 'twister');
            P0 = randn(size(DistSq, 1), cfg.opts.rank);

        case "Floyd"
            [~, P0, keep] = floyd_warshall_init(DistSq, Weight, 3);
            if numel(keep) ~= size(DistSq, 1)
                error('Floyd initialization trimmed nodes. Paper runner expects a connected graph.');
            end

        case "AF_3D"
            P0 = Pinit3;

        case "AF_rank"
            P0 = af_embed_highdim(Pinit3, cfg.opts.rank, cfg.init_seed, cfg.af_rank_jitter);

        otherwise
            error('Unknown initialization: %s', method);
    end
end
