function interneuron_decompose(basic_cell, decomp, approx, ...
                               time_stamp, dir_decompose, data_decompose, ...
                               varargin)
% Decompose the recurrent matrices
% Mengchen Zhu

%% Parse inputs
% Default parameter values
options = struct('rpca_file', NaN, 'S_threshold', NaN);

options = parse_inputs(options, varargin);

%% Local parameters
pd_decompose = options;
pd_decompose.decomp = decomp;
pd_decompose.approx = approx;
pd_decompose.gamma_rpca = NaN;
pd_decompose.S_sparsity = NaN;
pd_decompose.S_nonzero_rows = NaN;
pd_decompose.L_nuclearnorm = NaN;

switch decomp
  case 'svd'
    % SVD; note here to make the V here the equivalent of V in the
    % report111009, we switched the sequence of output.
    [V, Sigma, U] = svd(basic_cell.basis' * basic_cell.basis, 0);
    % Decompose matrices into positive and negative components
    Phi_p = basic_cell.basis;
    Phi_p(Phi_p < 0) = 0 ;
    Phi_n = basic_cell.basis;
    Phi_n(Phi_n > 0) = 0 ;
    V_p = V;
    V_p(V_p < 0) = 0;
    V_n = V;
    V_n(V_n > 0) = 0;    
    
  case {'arpca'}
    load(options.rpca_file);
    pd_decompose.gamma_rpca = gamma_rpca;
    
    [U_L, S_L, V_L] = svd(L, 0);
    
        
    diag_S_L = diag(S_L);    
    pd_decompose.L_nuclearnorm = sum(diag_S_L);
    

    switch approx
      case 'approx_RPCA'
        S(abs(S)<options.S_threshold) = 0; 
% $$$         S_L = diag([diag_S_L(1:options.SVD_terms); zeros(size(S_L,1) ...
% $$$                                                          - ...
% $$$                                                          options ...
% $$$                                                          .SVD_terms,1)]);
% $$$         approx_matrix = U_L * S_L * V_L' + S - eye(size(S_L)); 
    end


    
    % The projection to the inhibitory subpopulation of the sparse
    % interneurons.
    S_p = S;
    S_p(S_p < 0) = 0;
    
    pd_decompose.S_sparsity = numel(find(S_p ~= 0));
    % The indices of the non-zero rows
    nonzero_idx = find(all(S_p == 0,2) == 0);
    pd_decompose.S_nonzero_rows = numel(nonzero_idx);
    % The indices of the non-zero columns
    nonzero_idx = find(all(S_p == 0,1) == 0);
    pd_decompose.S_nonzero_cols = numel(nonzero_idx);
    
  case 'rpca_Phi'
    % TODO: need to pass in the gamma parameter.
    [L_temp S_temp iter] = inexact_alm_rpca(basic_cell.basis', gamma_rpca);
    L = L_temp';
    S = S_temp';
    
    [U_L, S_L, V_L] = svd(L, 0);
    % The projection to the inhibitory subpopulation of the low rank interneurons.
    V_L_p = V_L;
    V_L_p(V_L_p < 0) = 0;
    
    % The projection to the inhibitory subpopulation of the sparse
    % interneurons.
    S_p = S;
    S_p(S_p < 0) = 0;

end

switch decomp
  case 'svd'
    save(data_decompose, 'pd_decompose', 'Sigma', 'V', 'U');
  case {'rpca', 'rpca_Phi', 'arpca'} 
    rpca_file = options.rpca_file;
    % Do not store S since it's in rpca file
    save(data_decompose, 'pd_decompose', 'U_L', 'S_L', 'V_L');
  otherwise
    save(data_decompose, 'pd_decompose');
end

%% Log
% check if the log file already exists.
log_file = fullfile(dir_decompose, 'decompose_log.csv');
write_log_file(log_file, time_stamp, {pd_decompose}, 'sort_pars', true);

