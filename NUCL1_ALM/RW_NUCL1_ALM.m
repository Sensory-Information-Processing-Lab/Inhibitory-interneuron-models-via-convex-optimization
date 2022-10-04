function [L, S, niter, gamma]  = RW_NUCL1_ALM(X, gamma, tol, maxiter, ...
                                       alpha, beta, rw_iter, varargin)

if length(varargin) == 0
    fileid = 1;
else
    fileid = varargin{1};
end

for kk = 1:rw_iter
    [L, S, niter]  = NUCL1_ALM(X, gamma, tol, maxiter, fileid);
    gamma = alpha./(sum(abs(S), 2) + beta);
    fprintf(fileid, 'Finished Iteration %d\n', kk);
end

end