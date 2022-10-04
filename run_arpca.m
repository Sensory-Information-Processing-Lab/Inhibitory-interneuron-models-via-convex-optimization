%% Run RPCA
clc
clear all
close all

addpath NUCL1_ALM
%% Load basis
pd_basis = struct;
pd_basis.basis_size =  'double';
load('basic16x4oc50Ki_basis_size_double.mat');

%% Preprocessing
recurr_mat = basic_cell.basis' * basic_cell.basis;

%% Adaptive RPCA: X = L + S, where S is column sparse.

dir_interneuron = fullfile('.', 'data');
write_file = false;
if write_file
    % open file for writing
    fileid = fopen(fullfile(dir_interneuron, 'alm_rpca_log.log'), 'a+');
else
    fileid = 1; % std
end
gamma_rpca = 0.038*ones(size(basic_cell.basis,2), 1);
alpha = 2.5;
beta = 0.01;
rw_iter = 15;
tol = 0.01;
maxiter = 300;
L = []; S = [];

fprintf(fileid, ['___________________________ \n gamma %f, tol ' ...
                 '%f, alpha %f, beta %f.\n'], gamma_rpca(1), tol, alpha, beta);

[L, S, ~, gamma_rpca]  = RW_NUCL1_ALM(recurr_mat, gamma_rpca, tol, maxiter, ...
                                      alpha, beta, rw_iter, fileid);
rpca_file = fullfile(dir_interneuron, ['rpca_rw_alm_', ...
                    pars2str(gamma_rpca), '_', pars2str(tol), '_', ... 
                    pars2str(alpha), '_', pars2str(beta), '.mat']);
save(rpca_file, 'L', 'S', 'gamma_rpca');

