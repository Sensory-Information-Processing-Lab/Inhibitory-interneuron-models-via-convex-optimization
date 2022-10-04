function rl2e_sparsity(time_stamp, dir_decompose, data_decompose, ...
                       decomp_choice, basic_cell, pd_basis, ...
                       opt_stim_pars, varargin)
% Analyze the rL2E (relative l2 errors), sparsity, and energy function in response to
% optimal gratings.
% opt_stim_pars: {{1st_pass_lca_par 2nd_pass_lca_par}, pd_stim, rfpos} 
% varargin
%   with_S: Whether to include the sparse matrix in the
%   approximation.  

% is_perturb: whether to add random noise to the weight matrix

% This is to test that with small numbers of iterations, LCA has low
% RL2E and relatively high sparsity (i.e. the initial convergence is
% fast). Also used to study the approximation accuracy when
% interneurons are introduced (through e.g. RPCA).

% Mengchen Zhu

%% Parse inputs
% Default parameter values
options = struct('lca_pars', 'os1_s', 'with_S', true, 'whitened', ...
                 'whitened', 'SVD_terms', 100, 'plot_convergence', ...
                 false, 'is_perturb', false, 'rng', 100, 'perturb', ...
                 0.005, 'multiplicative', false, 'inhib_dyn', false, ...
                 'inhibTauScale', 1);

options = parse_inputs(options, varargin);

% $$$ %% Load RPCA decomposition
% $$$ dir_decompose = fullfile('..','data', 'interneuron', 'decompose');
% $$$ data_decompose = fullfile(dir_decompose, ['decompose_', time_stamp, ...
% $$$                     '.mat']);

load(data_decompose);

if (strcmp(decomp_choice, 'rpca') || (strcmp(decomp_choice, 'arpca')))
    load(pd_decompose.rpca_file);
end

%% Load optimal stimulus parameters
% First pass pars
% Local control paramters
pd_opt_sti_1p = struct;
ph_opt_sti_1p = struct;

% Load 1st pass paramaters for lca
[pd_opt_sti_1p ph_opt_sti_1p] = load_lca_par(pd_opt_sti_1p, ...
                                             ph_opt_sti_1p, ...
                                             opt_stim_pars{1}{1});
% Stimulus contrast 
pd_opt_sti_1p.contrast_idx = 4;
% double or normal basis size
pd_opt_sti_1p.basis_size = pd_basis.basis_size;

% Second pass pars
pd_opt_sti_2p = struct;
ph_opt_sti_2p = struct;

% Load 2nd pass paramaters for lca
[pd_opt_sti_2p ph_opt_sti_2p] = load_lca_par(pd_opt_sti_2p, ...
                                             ph_opt_sti_2p, ...
                                             opt_stim_pars{1}{2});

% Parameters from the second pass to display
pdisp_2p = struct;
pdisp_2p.itr_2p = pd_opt_sti_2p.itr;

load(fullfile('data', ['opt_par_', ...
                    pars2str(opt_stim_pars{2}), '_', ...
                    pars2str(pd_opt_sti_1p), '_', pars2str(pdisp_2p), ...
                    '.mat']));

%% Local Parameters
% Local control paramters
pd_cvg = options;
% Local hidden control paramters
ph_cvg = struct;

% Load lca parameters

[pd_cvg ph_cvg] = load_lca_par(pd_cvg, ph_cvg, options.lca_pars);

% Copy over double or normal basis size
pd_cvg.basis_size = pd_basis.basis_size;

% Identifies the threshold used in opt pars. 
pd_cvg.thr_2p = pd_opt_sti_2p.threshold;

% $$$ % Use static or dynamic drifting
% $$$ drifting = {'static'; 'dynamic'};
% $$$ pd_cvg.drifting = drifting{2};

pd_cvg.approx = pd_decompose.approx;
switch decomp_choice
  case 'svd'
    pd_cvg.with_S = NaN;
    pd_cvg.gamma_rpca = NaN;
    pd_cvg.S_threshold = NaN;
    diag_S = diag(Sigma);    
    
    if options.is_perturb
        rng(options.rng);
        if options.multiplicative
            U = U.*normrnd(1,options.perturb, size(U));
            V = V.*normrnd(1,options.perturb, size(V));            
        else
            U = U+normrnd(0,options.perturb, size(U));
            V = V+normrnd(0,options.perturb, size(V));
        end
    end
        
    approx_matrix = U * diag([diag_S(1:pd_cvg.SVD_terms); 
                        zeros(size(Sigma,1)-pd_cvg.SVD_terms,1)]) * V'- eye(size(Sigma));
  case {'rpca', 'arpca'}
    switch pd_cvg.approx
      case 'approx_RPCA'
        % Make sure the approx strings match
        % Threshold for pruning the sparse; NaN indicates no threshold;
        pd_cvg.S_threshold = pd_decompose.S_threshold;
        pd_cvg.gamma_rpca = pd_decompose.gamma_rpca;
        diag_S_L = diag(S_L);    
        
        if pd_cvg.with_S
            if ~isnan(pd_cvg.S_threshold)
                S(abs(S)<pd_cvg.S_threshold) = 0;
            end
            
            if options.is_perturb
                rng(options.rng);
                if options.multiplicative
                    U_L = U_L.*normrnd(1,options.perturb, size(U_L));
                    V_L = V_L.*normrnd(1,options.perturb, size(V_L));
                else
                    U_L = U_L+normrnd(0,options.perturb, size(U_L));
                    V_L = V_L+normrnd(0,options.perturb, ...
                                      size(V_L));
                end
                non_zero_rows = find(any(S,2));
                rng(options.rng);
                if options.multiplicative
                    randmat = normrnd(1,options.perturb, ...
                                      [size(non_zero_rows,1), ...
                                       size(S,2)]);
                    for row = 1:size(non_zero_rows)
                        S(non_zero_rows(row),:) = S(non_zero_rows(row),:) .* ...
                            randmat(row,:);
                    end
                else
                    randmat = normrnd(0,options.perturb, [size(non_zero_rows,1), size(S,2)]);
                    for row = 1:size(non_zero_rows)
                        S(non_zero_rows(row),:) = S(non_zero_rows(row),:) + ...
                            randmat(row,:);
                    end
                end
            end
            
            approx_matrix = U_L * diag([diag_S_L(1:pd_cvg.SVD_terms); ...
                                zeros(size(S_L,1)  - pd_cvg.SVD_terms,1)]) ...
                * V_L' + S - eye(size(S_L));        
        else
            approx_matrix = U_L * diag([diag_S_L(1:pd_cvg.SVD_terms); ...
                                zeros(size(S_L,1)  - pd_cvg.SVD_terms,1)]) ...
                * V_L' - eye(size(S_L));
        end
        if options.inhib_dyn
            U_p = U_L;
            U_p(U_p <= 0) = 0;
            U_n = U_L;
            U_n(U_n > 0) = 0;
            V_p = V_L;
            V_p(V_p <= 0) = 0;
            V_n = V_L;
            V_n(V_n > 0) = 0;
            Sigma = diag([diag_S_L(1:pd_cvg.SVD_terms); ...
                                zeros(size(S_L,1)  - ...
                                      pd_cvg.SVD_terms,1)]);
            S_p = S;
            S_p(S_p <= 0) = 0;
            S_n = S;
            S_n(S_n > 0) = 0;
            
            G_E = U_n * Sigma * V_p' + U_p * Sigma * V_n' + S_n;
            
            % Remove the zero eigenvalues
            U_n = U_n(:, 1:pd_cvg.SVD_terms);
            U_p = U_p(:, 1:pd_cvg.SVD_terms);
            Sigma = Sigma(1:pd_cvg.SVD_terms, 1:pd_cvg.SVD_terms);
            V_n = V_n(:, 1:pd_cvg.SVD_terms);
            V_p = V_p(:, 1:pd_cvg.SVD_terms);
            
% $$$             % Remove the zero columns
% $$$             idxSp = find(any(S'));
% $$$             S_p = S_p(idxSp, :);


            
% $$$             approx_matrix = U_p * Sigma * V_p' + (-U_n) * Sigma * ...
% $$$                 (-V_n') + S_p + G_E;
% $$$             
% $$$             approx_mat2 = U_L * diag([diag_S_L(1:pd_cvg.SVD_terms); ...
% $$$                                 zeros(size(S_L,1)  - pd_cvg.SVD_terms,1)]) ...
% $$$                 * V_L' + S;
            
        end
        
      otherwise
        pd_cvg.with_S = NaN;
        pd_cvg.SVD_terms = NaN;
        pd_cvg.gamma_rpca = NaN;
        pd_cvg.S_threshold = NaN;
    end
end
% Copy over hidden parameters. This includes for example the
% oversample switch, whether to display lca, and the oversampled
% radius.
ph_cvg = ph_opt_sti_2p;

%% Load stimulus
load(fullfile('data', strcat(pd_cvg.whitened,'_stimulus_', ...
                                   pars2str(opt_stim_pars{2}), '.mat')));

% the data file with the time stamp
dir_rl2e_data = fullfile('data',  'rl2e_sparsity');
if (~exist(dir_rl2e_data,'dir'))
    mkdir(dir_rl2e_data);
end

data_rl2e = fullfile(dir_rl2e_data, strcat(time_stamp, '.mat'));


dir_rl2e = fullfile('plots',   'rl2e_sparsity');
if (~exist(dir_rl2e,'dir'))
    mkdir(dir_rl2e);
end


% $$$ % time stamp up to second resolution in iso8610 format.
% $$$ time_stamp = datestr(now, 30);

% Make the subfolder with the time stamp
subdir_rl2e = fullfile(dir_rl2e, time_stamp);
if (~exist(subdir_rl2e,'dir'))
    mkdir(subdir_rl2e);
end


%% Run the experiment
for p = 1:size(opt_par, 2) %For every basis with optimum stimulus
                           %calculated.
% $$$ for p = 1:10
    % Parfor variables
    a_temp = zeros(size(ph_stim.stim_contrast, 2), ...
                   size(ph_cvg.radius,2), size(ph_stim.stim_phase, 2));
    
    idx_basis = opt_par(1,p);
    disp (strcat('Evaluating response for basis ', ...
                 num2str(idx_basis), ', count ',num2str(p)));
    center = round([opt_stim_pars{3}(2, idx_basis), opt_stim_pars{3}(1, idx_basis)]);
    % Generate the array of stimuli. Use the optimal
    % parameters. Use the same contrast as used in opt par.
    masked_stimulus = gen_masked_stim(stimulus(pd_opt_sti_1p.contrast_idx, ...
                                               opt_par(2,p), ...
                                               opt_par(3,p), ...
                                               opt_par(4,p),:,:), ...
                                      opt_par(5,p), center, ph_cvg);
    img = reshape(squeeze(masked_stimulus), [], 1);
    if (((strcmp(decomp_choice, 'rpca') || (strcmp(decomp_choice, 'arpca'))) && (strcmp(pd_cvg.approx, ...
                                                 'approx_RPCA') || ...
                                          strcmp(pd_cvg.approx, ...
                                                 'prune_PhiPhi'))) ...
        || strcmp(decomp_choice, 'svd'))
        [coef(1,:) u rl2e(1,:) sparsity(1,:) energy(1,:)] = lca(img, ...
                                                          basic_cell.basis, ...
                                                          pd_cvg, ...
                                                          ph_cvg, ...
                                                          'energy', 1);
        if options.inhib_dyn
            [coef(2,:) u rl2e(2,:) sparsity(2,:) energy(2,:)] = ...
                lca_inhib(img, basic_cell.basis, U_p, U_n, Sigma, ...
                          V_p, V_n, S_p, G_E, pd_cvg, 'energy', 1, ...
                          'inhibTauScale', options.inhibTauScale);
        else
            [coef(2,:) u rl2e(2,:) sparsity(2,:) energy(2,:)] = lca(img, ...
                                                          basic_cell.basis, ...
                                                          pd_cvg, ...
                                                          ph_cvg, ...
                                                          'energy', ...
                                                          1, 'approx', ...
                                                          approx_matrix);
        end
        % Record the steady state values
        sparsity_population(:, p) = sparsity(:, end);
        rl2e_population(:, p) = rl2e(:, end);
        energy_population(:, p) = energy(:, end);
        l2_diff_population(p) = norm(coef(1)-coef(2));
        disp(sparsity_population(:, p));
        disp(rl2e_population(:, p));

    elseif ((strcmp(decomp_choice, 'rpca') || strcmp(decomp_choice, 'arpca')) && strcmp(pd_cvg.approx, ...
                                                    'no'))
        [coef u rl2e sparsity energy] = lca(img, basic_cell.basis, pd_cvg, ...
                                            ph_cvg, 'energy', 1);
    else
        error('Not supported.')
    end
    
    % Plot the RL2E, sparsity, and the energy
    if options.plot_convergence
        h = figure;
        set(h, 'Position', [0,0, 300, 300]);
        subplot(3,1,1)
        plot([1:pd_cvg.itr], rl2e);
        xlabel('iterations')
        ylabel('rL2E')
        legend('original','approx', 'Location', 'NorthEast')
        
        subplot(3,1,2)
        plot([1:pd_cvg.itr], sparsity);
        hold all
        xlabel('iterations')
        ylabel('sparsity')    
        legend('original','approx', 'Location', 'NorthEast')
        
        subplot(3,1,3)
        plot([1:pd_cvg.itr], energy);
        hold all
        xlabel('iterations')
        ylabel('energy function')    
        legend('original','approx', 'Location', 'NorthEast')


        set(h, 'Color','w')
        export_fig(h, fullfile(subdir_rl2e, strcat(time_stamp, '_idx_', ...
                                                   num2str(p), '.pdf')));
        close(h)
    end
end

save(data_rl2e, 'sparsity_population', 'rl2e_population', ...
     'l2_diff_population', 'energy_population', 'pd_cvg');

% $$$ if strcmp(pd_cvg.approx, 'approx_RPCA') | strcmp(pd_cvg.approx, ...
% $$$                                              'prune_PhiPhi')
    
%% Analyze the population metric
% $$$ load(data_rl2e);
% $$$     time_stamp = '20111116T124005';
% $$$     load(fullfile(dir_rl2e_data, strcat(time_stamp,'.mat')));
% $$$     subdir_rl2e = fullfile(dir_rl2e, time_stamp);


h1 = figure;
set(h1, 'Position', [0, 0, 300, 300]);
scatter(sparsity_population(1,:), sparsity_population(2,:));
hold on
max_scatter = max(max(sparsity_population));
xlabel('# of active elements, original', 'FontSize', 14);
ylabel('# of active elements, approx', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 max_scatter]);
ylim([0 max_scatter]);
plot([0:max_scatter], [0:max_scatter],'k')
set(gcf, 'Color', 'w');
export_fig(h1, fullfile(subdir_rl2e, strcat('sparsity_', time_stamp, ...
                                           '.pdf')));  
close(h1)

h1 = figure;
set(h1, 'Position', [0, 0, 300, 300]);
scatter(rl2e_population(1,:), rl2e_population(2,:));
hold on
min_scatter = min(min(rl2e_population));
max_scatter = max(max(rl2e_population));
% $$$ plot([0:1], [0:1],'k')
% $$$ ylim(xlim)
xlim([min_scatter max_scatter]);
ylim([min_scatter max_scatter]);
plot(xlim, ylim, 'k');
% $$$ xlim([0 1]);
% $$$ ylim([0 1]);
xlabel('rL2E, original', 'FontSize', 14);
ylabel('rL2E, approx', 'FontSize', 14);
set(gca, 'FontSize', 14);
set(gcf, 'Color', 'w');
export_fig(h1, fullfile(subdir_rl2e, strcat('rl2e_', time_stamp, ...
                                           '.pdf')));  
close(h1)

h1 = figure;
set(h1, 'Position', [0, 0, 300, 300]);
scatter(energy_population(1,:), energy_population(2,:));
hold on
plot([0:max(max(energy_population))], [0:max(max(energy_population))],'k')
xlim([0 max(max(energy_population))+10]);
ylim([0 max(max(energy_population))+10]);
% $$$ axis equal
xlabel('Energy function, original', 'FontSize', 14);
ylabel('Energy function, approx', 'FontSize', 14);
% $$$ title(decomp_choice)
set(gca, 'FontSize', 14);
set(gcf, 'Color', 'w');
export_fig(h1, fullfile(subdir_rl2e, strcat('energy_', time_stamp, ...
                                           '.pdf')));  
close(h1)

rel_energy_diff = (energy_population(2,:)-energy_population(1,:))./ ...
               energy_population(1,:);
disp('Relative energy difference:')
disp([num2str(mean(rel_energy_diff)), '+-', num2str(var(rel_energy_diff))])

h1 = figure;
set(h1, 'Position', [0, 0, 300, 300]);
hist(l2_diff_population);
xlabel('|a - a_{approx}|_2');
ylabel('# of units');
set(gcf, 'Color', 'w');
export_fig(h1, fullfile(subdir_rl2e, strcat('l2_diff_', time_stamp, ...
                                           '.pdf')));  
close(h1)

% $$$ end

%% Log
log_file = fullfile(dir_rl2e,'rl2e_log.csv');
write_log_file(log_file, time_stamp, {pd_stim pd_opt_sti_1p pdisp_2p pd_cvg});
