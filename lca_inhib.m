function [a varargout] = lca_inhib(s, Phi, U_p, U_n, Sigma, V_p, ...
                                   V_n, S_p, G_E, pars, varargin)
% LCA with inhibitory interneuron dynamics


%% Parse inputs
u_init = zeros(size(Phi,2),1);
a_init = u_init; 
a_I1 = ones(size(Sigma,1), 1);
a_I2 = a_I1;
% $$$ a_I1 = u_init;
% $$$ a_I2 = u_init;
a_S = u_init;
% $$$ options = struct('init', [u_init a_init a_I1 a_I2 a_S] ,'energy', ...
% $$$                  0);
options = struct('init', [u_init a_init] ,'energy', 0, 'inhibTauScale', ...
                 1);

options = parse_inputs(options, varargin);


%% LCA
u = options.init(:,1);
a = options.init(:,2);
% $$$ a_I1 = options.init(:,3);
% $$$ a_I2 = options.init(:,4);
% $$$ a_S = options.init(:,5);

for n = 1:pars.itr
    if options.energy
        rl2e(n) = norm(s-Phi * a)/norm(s);
        sparsity(n) = nnz(a);
        energy(n) = 0.5*norm(s-Phi * a)^2 + pars.threshold*sum(abs(a));
    end
    %Euler's method to solve difference eqn
% $$$     u = u*(1-pars.delta/pars.tau) + pars.delta/pars.tau*(Phi' * s - (U_p * Sigma * V_p' + (-U_n) ...
% $$$                                                       * Sigma * (-V_n') + S_p ) * a  ...
% $$$                                                       + (eye(size(G_E)) ...
% $$$                                                       - G_E) * a);
    c_inhib = options.inhibTauScale*(pars.delta/pars.tau); % delta/tau
% $$$     a_I1 = a_I1 * (1-c_inhib) + c_inhib * (Sigma * V_p' * a);
% $$$     a_I2 = a_I2 * (1-c_inhib) + c_inhib * (Sigma * (-V_n') * ...
% $$$                                                a);
% $$$     a_S  = a_S * (1-c_inhib) + c_inhib * a;

    % exponential decay to the solution
    a_I1 = a_I1 * exp(-n*c_inhib) + Sigma * V_p' * a;
    a_I2 = a_I2 * exp(-n*c_inhib) + Sigma * (-V_n') * a;
    a_S  = a_S * exp(-n*c_inhib) + a;

    u = u*(1-pars.delta/pars.tau) + pars.delta/pars.tau*(Phi' * s - (U_p * a_I1 + (-U_n) ...
                                                 * a_I2 + S_p * a_S) ...
                                                 + (eye(size(G_E)) ...
                                                    - G_E) * a);
    if isnan(u)
        error('NaN encountered in LCA.')
    end
            
    a = T_lambda(u, pars);     
end

%% Outputs
if options.energy
    varargout = {u, rl2e, sparsity, energy};
else
    varargout = {u};
end
