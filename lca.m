function [a varargout] = lca(s, Phi, pd, ph, varargin)
% LCA
% Mengchen Zhu. 

% Usage:
%   a = lca(s, Phi, pd, ph)
%   [a u] = lca(s, Phi, pd, ph, 'init', [u_init a_init])
%   [a u rl2e sparsity] = lca(s, Phi, pd, ph, 'init', [u_init
%   a_init], 'energy', 1)

% inputs:
%   s: Input image. *Should be 1D*: e.g. s = reshape(squeeze(stimulus(m,n,:,:)),
%     [], 1);
%   Phi: Dictionary. Column major. *Should be unit norm*.
%   pd: 
%     basis_size: whether double or normal dictionary size
%     lca pars: tau; delta; typically delta/tau = 0.1; determines the
%       step size.
%     threshold; incorporates lambda; the larger the threshold, the
%       sparser.
%     itr: number of iterations
%     thr_func: type of thresholding function
%   ph:
%     display: whether to display results
%     h_lca: handle to figure
%   varargin:
%     'init':
%       u_init: initial conditions for u;
%       a_init: initial conditions for a;
%     'energy':
%       whether to include RL2E and sparsity calculation.
%     'approx':
%       The approximation to the recurrent matrix Phi'*Phi - I.
% outputs:
%   a: the response
%   varargout:
%     u: the hidden variable.
%     rl2e, sparsity, energy: relative l2 error, sparsity, and the energy function w.r.t. iteration.


%% Parse inputs
u_init = zeros(size(Phi,2),1);
a_init = u_init; %Thresholded u.
options = struct('init', [u_init a_init] ,'energy', 0, 'approx', NaN);

options = parse_inputs(options, varargin);


%% LCA
u = options.init(:,1);
a = options.init(:,2);
   
if isnan(options.approx)
    % Default recurrent matrix
    approx = Phi'*Phi - eye(size(Phi, 2));
else
    approx = options.approx;
end

for n = 1:pd.itr
    if options.energy
        rl2e(n) = norm(s-Phi * a)/norm(s);
        sparsity(n) = nnz(a);
        energy(n) = 0.5*norm(s-Phi * a)^2 + pd.threshold*sum(abs(a));
    end
    %Euler's method to solve difference eqn
    u = u*(1-pd.delta/pd.tau) + pd.delta/pd.tau*(Phi' * s - approx * a);
    if isnan(u)
        error('NaN encountered in LCA.')
    end
            
    switch pd.thr_func
      case 'block'
        % The block l2 norm is the amplitude of (complex) u.
        u_real = u(1:(length(u)/2));
        u_img = u((length(u)/2 + 1):end);
        u_cplx = complex(u_real, u_img);
        u_block = [abs(u_cplx); abs(u_cplx)];
        a = T_lambda(u, pd, 'u_block', u_block);     
      otherwise
        a = T_lambda(u, pd);     
    end
end

%% Outputs
if options.energy
    varargout = {u, rl2e, sparsity, energy};
else
    varargout = {u};
end

if ph.display_lca
    figure(ph.h_lca);
    subplot(2,2,[1 2])
    stem(a)
    title(pd.itr)
    
    %Validation: regenerate the input image

    s_r = Phi * a;
    s_r = reshape(s_r, sqrt(size(s,1)),sqrt(size(s,1)));
    s = reshape(s, sqrt(size(s,1)),sqrt(size(s,1)));

    subplot(2,2,3)
    %imagesc(s_r, pars.graylevel)
    imagesc(s_r)
    colormap gray
    title('lca')
    axis equal
    subplot(224)
    %imagesc(s, pars.graylevel)
    imagesc(s)
    colormap gray
    title('real')
    axis equal
end
