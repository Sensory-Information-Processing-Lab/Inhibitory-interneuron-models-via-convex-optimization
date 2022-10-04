function map_RF(time_stamp, basic_cell, pd_basis, decomp, varargin)
% Use spots to map out RF.
% Inputs:
%   basic_cell, pd_basis: basis parameters
%   decomp: decomposition methods: {'svd';'false'; 'rpca';
%   'grammian'; 'global'; 'rpca_Phi'};


%   * Choose appropriate stimulus contrast.
%   * For rpca, choose the gamma; for global, choose the c.

% Mengchen Zhu

%% Parse inputs
% Default parameter values
options = struct('lca_pars', 'os1_s', 'whitened', 'whitened', ...
                 'stim_contrast', 1, 'SVD_terms', NaN, 'c_global', NaN);

options = parse_inputs(options, varargin);

%% Load RPCA decomposition
dir_decompose = fullfile('data');
data_decompose = fullfile(dir_decompose, ['decompose_', time_stamp, ...
                    '.mat']);
load(data_decompose);

%% Local parameters
% Local control paramters
pd_mapRF = options;
% Local hidden control paramters
ph_mapRF = struct;

% Load lca parameters
[pd_mapRF ph_mapRF] = load_lca_par(pd_mapRF, ph_mapRF, options.lca_pars);

% Copy over double or normal basis size
pd_mapRF.basis_size = pd_basis.basis_size;

% Load stimulus; 
if strcmp(pd_mapRF.whitened, 'whitened')
    load(fullfile('data','whitened_spot_stim.mat'));
end

pd_mapRF.decomp = decomp;
pd_mapRF.approx = pd_decompose.approx;

switch pd_mapRF.decomp
  case  {'rpca', 'arpca'}
    % Load the RPCA decomposition.
    load(pd_decompose.rpca_file);
    pd_mapRF.gamma_rpca = pd_decompose.gamma_rpca;
  case 'rpca_Phi'
    pd_mapRF.gamma_rpca = 0.2;
  otherwise
    pd_mapRF.gamma_rpca = NaN;
end

switch pd_mapRF.approx
  case 'approx_RPCA'
    pd_mapRF.S_threshold = pd_decompose.S_threshold;
end

% Display LCA reconstruction?
ph_mapRF.display_lca = 0;
if ph_mapRF.display_lca
    ph_mapRF.h_lca = figure;
end

% First check if the directory already exists.
dir_mapRF = fullfile('data', 'mapRF');
if (~exist(dir_mapRF,'dir'))
    mkdir(dir_mapRF);
end
% the data file with the time stamp
data_mapRF = fullfile(dir_mapRF, ['mapRF_', time_stamp, '.mat']);


%% Calculate decomposition
switch pd_mapRF.decomp
  case 'svd'
    % SVD; note here to make the V here the equivalent of V in the
    % report111009, we switched the sequence of output.
    [V, Sigma, U] = svd(basic_cell.basis', 0);
    % Decompose matrices into positive and negative components
    Phi_p = basic_cell.basis;
    Phi_p(Phi_p < 0) = 0 ;
    Phi_n = basic_cell.basis;
    Phi_n(Phi_n > 0) = 0 ;
    V_p = V;
    V_p(V_p < 0) = 0;
    V_n = V;
    V_n(V_n > 0) = 0;    
  case {'rpca', 'arpca'}
    if strcmp(pd_mapRF.approx, 'approx_RPCA')
        diag_S_L = diag(S_L);    
        S(abs(S)<pd_mapRF.S_threshold) = 0; 
        S_L = diag([diag_S_L(1:pd_mapRF.SVD_terms); zeros(size(S_L,1) ...
                                                          - pd_mapRF .SVD_terms,1)]);
        approx_matrix = U_L * S_L * V_L' + S - eye(size(S_L)); 
    end
    
    % The projection to the inhibitory subpopulation of the low rank interneurons.
    V_L_p = V_L;
    V_L_p(V_L_p < 0) = 0;

    % The projection to the inhibitory subpopulation of the sparse
    % interneurons.
    S_p = S;
    S_p(S_p < 0) = 0;

    switch pd_mapRF.decomp
      case 'arpca'
        % In the case of arpca, the sparse connections to the
        % interneurons are just identity matrix, with nonzero rows
        % matching the S matrix.
        temp_zeros = zeros(size(S,1),1);
        nonzero_idx = find(all(S_p == 0,2) == 0);
        temp_zeros(nonzero_idx) = 1;
        S_p = diag(temp_zeros);        
    end     
            
            
  case 'rpca_Phi'
    [L_temp S_temp iter] = inexact_alm_rpca(basic_cell.basis', ...
                                            pd_mapRF.gamma_rpca);
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

% $$$ mapRF = zeros(size(basis_cell.basis, 2), size(basic_cell.basis, 1));
%% Map the RF
% Calculate response for each of the spot stimulus
for stim_idx = 1:size(basic_cell.basis, 1)
% $$$ for stim_idx = 1:size(basic_cell.basis, 1)
    disp(strcat('Mapping RF for stimulus', num2str(stim_idx)));
    % Map out the positive part of the RF using negative stimuli.
    img_p = pd_mapRF.stim_contrast * ...
            reshape(squeeze(stimulus(stim_idx,:,:)), [], 1);
    if strcmp(pd_mapRF.approx, 'approx_RPCA')
        mapRF_p(:, stim_idx) = lca(img_p, basic_cell.basis, pd_mapRF, ...
                                   ph_mapRF, 'approx', approx_matrix);            
    else
        mapRF_p(:, stim_idx) = lca(img_p, basic_cell.basis, pd_mapRF, ...
                                   ph_mapRF);                
    end
    switch pd_mapRF.decomp
      case 'svd'
        % Calculate the "response" of the two
        % populations of inhibitory cells.                    
        inhib_response_1_p(:, stim_idx) = Sigma^2*(-V_n') * ...
            mapRF_p(:, stim_idx);
        inhib_response_2_p(:, stim_idx) = Sigma^2 * V_p' * ...
            mapRF_p(:, stim_idx);                    
      case {'rpca', 'arpca'}
        % Here we only look at one of the two inhibitory
        % populations? (111107)
        rpca_inhib_response_p(:, stim_idx) = S_L * V_L_p' * mapRF_p(:, ...
                                                          stim_idx);
        rpca_sparse_p(:, stim_idx) = S_p * mapRF_p(:, stim_idx);
      case 'grammian'
        grammian_inhib_response_p(:, stim_idx) = basic_cell.basis * ...
            mapRF_p(:, stim_idx);
      case 'global'
        global_inhib_response_p(stim_idx) = sum(mapRF_p(:, stim_idx), ...
                                                1);
      case 'rpca_Phi'
        rpca_inhib_response_p(:, stim_idx) = S_L * V_L_p' * mapRF_p(:, ...
                                                          stim_idx);
        rpca_sparse_p(:, stim_idx) = S_p * mapRF_p(:, stim_idx);
    end        
    
    % Map out the negative part of the RF using negative stimuli.
    img_n = - pd_mapRF.stim_contrast * ...
          reshape(squeeze(stimulus(stim_idx,:,:)), [], 1);
    mapRF_n(:, stim_idx) = lca(img_n, basic_cell.basis, pd_mapRF, ...
                                  ph_mapRF); 
    switch pd_mapRF.decomp
      case 'svd'
        % Calculate the "response" of the two
        % populations of inhibitory cells.                    
        inhib_response_1_n(:, stim_idx) = Sigma^2*(-V_n') * ...
            mapRF_n(:, stim_idx);
        inhib_response_2_n(:, stim_idx) = Sigma^2 * V_p' * ...
            mapRF_n(:, stim_idx); 
      case {'rpca', 'arpca'}
        rpca_inhib_response_n(:, stim_idx) = S_L * V_L_p' * mapRF_n(:, ...
                                                          stim_idx);
        rpca_sparse_n(:, stim_idx) = S_p * mapRF_n(:, stim_idx);
      case 'grammian'
        grammian_inhib_response_n(:, stim_idx) = basic_cell.basis * ...
            mapRF_n(:, stim_idx);
      case 'global'
        global_inhib_response_n(stim_idx) = sum(mapRF_n(:, stim_idx), 1);
      case 'rpca_Phi'
        rpca_inhib_response_n(:, stim_idx) = S_L * V_L_p' * mapRF_n(:, ...
                                                          stim_idx);
        rpca_sparse_n(:, stim_idx) = S_p * mapRF_n(:, stim_idx);
    end        
end

% Save the RF. Extract the response of bases in the range
% basis_range.
mapRF = mapRF_p - mapRF_n;

% $$$ pd_mapRF.completed = true;
% $$$ % New parameter values
% $$$ par_values = [{time_stamp}, struct2cell(pd_mapRF)'];
% $$$ log_cell(end, :) = par_values;
% $$$ cell2csv(fullfile(dir_mapRF, log_file), log_cell);

switch pd_mapRF.decomp
 case 'svd'
    inhib_response_1 = inhib_response_1_p - inhib_response_1_n;
    inhib_response_2 = inhib_response_2_p - inhib_response_2_n;
    save(data_mapRF, 'mapRF', 'inhib_response_1', 'inhib_response_2', ...
         'pd_mapRF', 'ph_mapRF');
  case {'rpca', 'rpca_Phi', 'arpca'} 
    rpca_inhib_response = rpca_inhib_response_p - ...
        rpca_inhib_response_n;
    rpca_sparse = rpca_sparse_n - rpca_sparse_p;
    save(data_mapRF, 'mapRF', 'rpca_inhib_response', 'pd_mapRF', 'ph_mapRF', ...
         'rpca_sparse');
  case 'grammian'
    grammian_inhib_response = grammian_inhib_response_p - ...
        grammian_inhib_response_n;
    save(data_mapRF, 'mapRF', 'grammian_inhib_response', 'pd_mapRF', 'ph_mapRF');
  case 'global'
    global_inhib_response = global_inhib_response_p - global_inhib_response_n;
    save(data_mapRF, 'mapRF', 'global_inhib_response', 'pd_mapRF', 'ph_mapRF');
  otherwise
    save(data_mapRF, 'mapRF', 'pd_mapRF', 'ph_mapRF');
end

%% Plot the RF
% First check if the directory already exists.
dir_mapRF_plot = fullfile('plots', 'mapRF');
if (~exist(dir_mapRF_plot,'dir'))
    mkdir(dir_mapRF_plot);
end

% Make subdirectory with time stamp
sub_dir = fullfile(dir_mapRF_plot, time_stamp);
if (~exist(sub_dir,'dir'))
    mkdir(sub_dir);
end



img_sz = sqrt(size(basic_cell.basis, 1)); %The width and height of
                                          %each basis
if strcmp(pd_mapRF.approx, 'approx_RPCA')
    basis_range = 1:pd_mapRF.SVD_terms;
else
    basis_range = 1:100;
end

disp_sz = [10 length(basis_range)/10];

% subset of FF RFs 
h = figure;
img = basis2img2(mapRF(basis_range, :)', [img_sz img_sz], disp_sz);
imagesc(img)
colormap gray
axis image
axis off
set(gcf, 'Color', 'w');
export_fig(h, fullfile(sub_dir, strcat('FF_RF_', time_stamp, '.pdf')));    
close(h)



h2 = figure;
img_dict = basis2img2(basic_cell.basis(:, basis_range), [img_sz, ...
                    img_sz], disp_sz);
imagesc(img_dict)
colormap gray
axis image
axis off
set(gcf, 'Color', 'w');
export_fig(h2, fullfile(sub_dir, strcat('dict_', time_stamp, '.pdf')));        
close(h2)


switch(pd_mapRF.decomp)
  case 'svd'
    h = figure;
    img = basis2img2(inhib_response_1(basis_range, :)', [img_sz, ...
                        img_sz], disp_sz);
    imagesc(img)
    colormap gray
    axis image
    axis off
    set(gcf, 'Color', 'w');
    export_fig(h, fullfile(sub_dir, strcat('svd_inhib_RF_1_', time_stamp, ...
                                           '.pdf')));    
    close(h)
   
    h = figure;
    img = basis2img2(inhib_response_2(basis_range, :)', [img_sz, ...
                        img_sz], disp_sz);
    imagesc(img)
    colormap gray
    axis image
    axis off
    set(gcf, 'Color', 'w');
    export_fig(h, fullfile(sub_dir, strcat('inhib_RF_2_', time_stamp, ...
                                           '.pdf')));    
    close(h)

  case {'rpca', 'rpca_Phi', 'arpca'}
    % Plot the low rank RFs.
    h = figure;
    img = basis2img2(rpca_inhib_response(basis_range, :)', [img_sz, ...
                        img_sz], disp_sz);
    imagesc(img)
    colormap gray
    axis image
    axis off
    set(gcf, 'Color', 'w');
    export_fig(h, fullfile(sub_dir, strcat('rpca_lowrank_RF_', time_stamp, ...
                                           '.pdf')));    
    close(h)

    
    % Plot sparse RFs.
    basis_range = 1:100;
    h = figure;
    % The indices of the non-zero rows
    nonzero_idx = find(all(S_p == 0,2) == 0);
    % Display the first 100 RFs in a panel.
    sparse_RF_subset = rpca_sparse(nonzero_idx, :);
    sparse_RF_panel = basis2img2(sparse_RF_subset(basis_range, :)', ...
                                 [img_sz img_sz], ...
                                 [sqrt(length(basis_range)), ...
                        sqrt(length(basis_range))], 1);
    imagesc(sparse_RF_panel)
    axis equal; axis off
    set(gcf, 'Color', 'w');
    export_fig(h, fullfile(sub_dir, strcat('sparse_RF_', time_stamp, '.pdf')));
    close(h)

% $$$     % testing: does the FF RFs match those of the sparse RF?
% $$$     h = figure;
% $$$     img = basis2img2(mapRF(nonzero_idx, :)', [img_sz img_sz], disp_sz);
% $$$     imagesc(img)
% $$$     colormap gray
% $$$     axis image
% $$$     axis off
% $$$     set(gcf, 'Color', 'w');
% $$$     export_fig(h, fullfile(sub_dir, strcat('FF_S_RF_', time_stamp, '.pdf')));    
% $$$     close(h)


  case 'grammian'
    h = figure;
    img = basis2img2(grammian_inhib_response(basis_range, :)', [img_sz, ...
                        img_sz], disp_sz);
    imagesc(img)
    colormap gray
    axis image
    axis off
    set(gcf, 'Color', 'w');
    export_fig(h, fullfile(sub_dir, strcat('grammian_inhib_RF_', ...
                                           time_stamp, '.pdf')));    
    close(h)
    
  case 'global'
    h = figure;
    img = reshape(global_inhib_response, img_sz, img_sz);
    imagesc(img)
    colormap gray
    axis image
    axis off
    set(gcf, 'Color', 'w');
    export_fig(h, fullfile(sub_dir, strcat('global_inhib_RF_', ...
                                           time_stamp, '.pdf')));
    close(h);
end


%% Log
% check if the log file already exists.
log_file = fullfile(dir_mapRF, 'mapRF_log.csv');
write_log_file(log_file, time_stamp, {pd_mapRF}, 'sort_pars', true);

% $$$ log_file = 'mapRF_log.csv';
% $$$ if(~exist(fullfile(dir_mapRF, log_file), 'file'))
% $$$     % If not, generate the file and the first row (parameter
% $$$     % names): strip out the field names as column headings
% $$$     log_cell = par_name;
% $$$     cell2csv(fullfile(dir_mapRF, log_file), log_cell);
% $$$ else
% $$$     % If exists, open it and read to a cell.
% $$$     log_cell = csv2cell(fullfile(dir_mapRF, log_file), 'fromfile'); ...
% $$$     % If the par_name does not agree with the old one, add a new
% $$$     % line of parameter names and new dimensions.
% $$$     if size(log_cell, 2)<size(par_name,2)
% $$$         log_cell = [log_cell, cell(size(log_cell,1), (size(par_name,2) ...
% $$$                                                       - size(log_cell,2))); ...
% $$$                     par_name];
% $$$     end
% $$$ end
% $$$ % Generate the folder name from the time stamp; log the time; 
% $$$ % Write the time stamp and the parameters to the last row. 
% $$$ log_cell = [log_cell; par_values];
% $$$ % Write the new log_file
% $$$ cell2csv(fullfile(dir_mapRF, log_file), log_cell);

% $$$ matlabpool close
