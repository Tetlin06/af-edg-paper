function [result, captured] = paper_run_cfg(cfg)
%PAPER_RUN_CFG Run one config while capturing console chatter from helpers.
    captured = evalc('result = edg_run_one(cfg);');
end
