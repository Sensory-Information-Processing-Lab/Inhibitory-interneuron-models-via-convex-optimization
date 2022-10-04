clc
clear all
close all

addpath NUCL1_ALM
%% Load basis
pd_basis = struct;
pd_basis.basis_size =  'double';
load('data/basic16x4oc50Ki_basis_size_double.mat');

%% Preprocessing
recurr_mat = basic_cell.basis' * basic_cell.basis;

decomp = {'svd'; 'arpca'};
decomp_choice = 'arpca';

dir_interneuron = fullfile('.', 'data');
if (~exist(dir_interneuron,'dir'))
    mkdir(dir_interneuron);
end
%% Analyse
% Decompose recurrent matrices
% time stamp up to second resolution in iso8610 format.
time_stamp = datestr(now, 30);
rpca_file = fullfile(dir_interneuron, 'rpca_rw_alm_gamma_rpca_250_0_250_tol_0p01_alpha_2p5_beta_0p01');


data_decompose = fullfile(dir_interneuron, ['decompose_', time_stamp, ...
                    '.mat']);

if (~exist(data_decompose))
    switch decomp_choice
      case {'arpca'}
        interneuron_decompose(basic_cell, decomp_choice, ...
                              'approx_RPCA', time_stamp, dir_interneuron, ...
                              data_decompose, 'rpca_file', rpca_file);
      case 'svd'
        interneuron_decompose(basic_cell, decomp_choice, 'no', ...
                              time_stamp, dir_interneuron, data_decompose);
    end    
else
    disp('Using old decomposition data!!!')
end


%% rL2E and sparsity

pd_stim = struct;
pd_stim.contrast = 'more';

% arpca
data_decompose = fullfile(dir_interneuron, ['decompose_', time_stamp, ...
                    '.mat']);
load(data_decompose);
decomp_choice = 'arpca';
SVD_terms = 110;
rl2e_sparsity(time_stamp, dir_interneuron, data_decompose, decomp_choice, ...
              basic_cell, pd_basis, {{'os1_s', 'os2_s'}, pd_stim, ...
                    rfpos}, 'lca_pars', 'dynamic_25', 'with_S', true, ...
              'whitened', 'whitened', 'SVD_terms', SVD_terms);

%% arpca with inhib dynamics

time_stamp_dyn = datestr(now, 30);
rl2e_sparsity(time_stamp_dyn, dir_interneuron, data_decompose, decomp_choice, ...
              basic_cell, pd_basis, {{'os1_s', 'os2_s'}, pd_stim, ...
                    rfpos}, 'lca_pars', 'dynamic_100', 'with_S', true, ...
              'whitened', 'whitened', 'SVD_terms', SVD_terms, ...
              'inhib_dyn', true, 'inhibTauScale', 1);
%% Map the RFs
map_RF(time_stamp, basic_cell, pd_basis, decomp_choice, 'lca_pars', 'dynamic_25', ...
       'whitened', 'whitened', 'stim_contrast', 1, 'SVD_terms', SVD_terms);
