function [opts, lsopts] = paper_default_solver_opts()
%PAPER_DEFAULT_SOLVER_OPTS Solver settings used by the paper-only scripts.
    opts = struct();
    opts.r           = 10000;
    opts.printenergy = 0;
    opts.printerror  = 0;
    opts.rank        = 10;
    opts.maxit       = 10000;
    opts.tol         = 1e-5;
    opts.lamda       = opts.r;

    lsopts = struct();
    lsopts.maxit = 30;
    lsopts.xtol  = 1e-8;
    lsopts.gtol  = 1e-8;
    lsopts.ftol  = 1e-10;
    lsopts.alpha = 1e-3;
    lsopts.rho   = 1e-4;
    lsopts.sigma = 0.1;
    lsopts.eta   = 0.8;
end
