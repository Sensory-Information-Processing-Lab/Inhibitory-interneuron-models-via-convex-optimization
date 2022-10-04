function [L, S, niter]  = NUCL1_ALM(M,lambda,tol,maxiter,varargin)
% 
% Finds the Principal Component Pursuit solution 
%
% minmize         ||L||_* + lambda || S ||_1 
%
% subject to      L + Phi * S = M
%
% using an augmented Lagrangian approach
%
% Usage:  [L,S,iter]  = RPCA_ALM(M,lambda,tol,maxiter,func)
%
% Inputs:
%
% M       - input matrix of size n1 x n2 
%
% lambda  - parameter defining the objective functional 
% 
% tol     - algorithm stops when ||M - L - Phi * S||_F <= delta ||M||_F 
%
% maxiter - maximum number of iterations
%
% varargin - fileid
%
% Outputs: 
% 
% L        - low-rank component
% 
% S        - sparse component
%
% niter     - number of iterations to reach convergence

% Created: January 2011

if length(varargin) == 0
    fileid = 1;
else
    fileid = varargin{1};
end

n = size(M); 

if nargin < 4 || isempty(maxiter),  maxiter = 500;           end
if nargin < 3 || isempty(tol),      tol = 1e-6;              end
if nargin < 2 || isempty(lambda),   lambda = 1/sqrt(max(n)); end

Frob_norm = norm(M,'fro');
two_norm = svds(M,1,'L');
one_norm = sum(abs(M(:)));
inf_norm = max(abs(M(:)));

mu_inv = 4*one_norm/prod(n);

% Kicking
k = min(min(floor(mu_inv/two_norm)), min(floor(lambda*mu_inv/inf_norm)));
Y = k*M; 

% Main loop
sv = 10;
S = zeros(n); 
L = S;

for k=1:maxiter,
    % Shrink entries
    S = vector_shrink(M - L + mu_inv*Y, lambda*mu_inv);
    
    % Shrink singular values 
    [L,r] = matrix_shrink(M - S + mu_inv*Y, mu_inv,sv);
    if r < sv
        sv = min(r + 1, min(n));
    else
        sv = min(r + round(0.05*min(n)), min(n));
    end
    
    stopCriterion = norm(M-L-S,'fro')/Frob_norm;
    
    if mod(k, 1) == 0
        fprintf(fileid, ['iter: %3d, rank(L) %3d, |S|_0: %6d, ' ...
                         'stopCriterion %d \n'],k,r,sum(S(:)~= ...
                                                        0),stopCriterion);
    end    
    
    % Check convergence
    if (stopCriterion < tol)
        break
    end
    
    % Update dual variable
    Y = Y + (M-L-S)/mu_inv;         
end

niter = k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxilliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y,r] = matrix_shrink(X,tau,sv)
 
m = min(size(X));

if choosvd(m,sv) == 1,
    [U, S, V] = lansvd(X,sv,'L');
else
    [U,S,V] = svd(X,0); 
end
sigma = diag(S); r = sum(sigma > tau);
if r == 0,
    Y = zeros(size(X));
else
    Y = U(:,1:r)*sparse(diag(sigma(1:r) - tau))*V(:,1:r)';
end


function Y = vector_shrink(X,tau)

if numel(tau) == 1
    Y = sign(X).* max(abs(X) - tau, 0);
elseif numel(tau) == size(X, 1)
    Y = sign(X).* max(abs(X) - tau(:)*ones(1, size(X, 2)), 0);
else
    error('something is wrong!!')
end









