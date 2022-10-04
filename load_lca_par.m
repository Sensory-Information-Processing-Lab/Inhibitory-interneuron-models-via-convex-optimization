function [pd ph] = load_lca_par(pd,ph,opt)
% Load parameters for computing and displaying LCA

switch opt
  case 'os1'
    % optimal stimulation 1st pass: sparse and fast
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.1; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 100;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'os2'
    % optimal stimulation 2nd pass: sparse
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.1; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 1000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'os1_ns'
    % optimal stimulation 1st pass; non-sparse
    pd.tau = 1; 
    pd.delta = 0.05; 
    pd.threshold = 0.005; 
    pd.itr = 100;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'os2_ns'
    % optimal stimulation 2nd pass: non-sparse
    pd.tau = 1; 
    pd.delta = 0.05; 
    pd.threshold = 0.005; 
    pd.itr = 1000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'os1_ns_small_step'
    % optimal stimulation 1st pass; non-sparse
    pd.tau = 1; 
    pd.delta = 0.025; 
    pd.threshold = 0.005; 
    pd.itr = 1000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'os2_ns_small_step'
    % optimal stimulation 2nd pass: non-sparse
    pd.tau = 1; 
    pd.delta = 0.025; 
    pd.threshold = 0.005; 
    pd.itr = 4000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'os1_s'
    % sparse; 
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.5; 
    pd.itr = 100;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'os2_s'
    % sparse; 
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.5; 
    pd.itr = 1000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'os1_s_200'
    % sparse; 
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.5; 
    pd.itr = 200;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_ns'
    % surround suppression; non-sparse
    pd.tau = 1; 
    pd.delta = 0.05; 
    pd.threshold = 0.005; 
    pd.itr = 2000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_ns_li'
    % surround suppression; non-sparse; less iterations
    pd.tau = 1; 
    pd.delta = 0.05; 
    pd.threshold = 0.005; 
    pd.itr = 200;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_ns_400'
    % surround suppression; non-sparse; less iterations
    pd.tau = 1; 
    pd.delta = 0.05; 
    pd.threshold = 0.005; 
    pd.itr = 400;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_ns_600'
    % surround suppression; non-sparse; less iterations
    pd.tau = 1; 
    pd.delta = 0.05; 
    pd.threshold = 0.005; 
    pd.itr = 600;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_ns_800'
    % surround suppression; non-sparse; less iterations
    pd.tau = 1; 
    pd.delta = 0.05; 
    pd.threshold = 0.005; 
    pd.itr = 800;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_s'
    % surround suppression; sparse
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.5; 
    pd.itr = 3000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_s_li'
    % surround suppression; sparse; less iterations
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.5; 
    pd.itr = 1000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_s_si'
    % surround suppression; sparse; small iterations
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.5; 
    pd.itr = 500;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_s_ssi'
    % surround suppression; sparse; very small iterations
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.5; 
    pd.itr = 100;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_ss'
    % surround suppression; very sparse
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 1; 
    pd.itr = 3000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_sss'
    % surround suppression; very very sparse
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 2; 
    pd.itr = 3000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'ss_m'
    % surround suppression; medium sparse
    pd.tau = 1; %2?
    pd.delta = 0.05; %?
    pd.threshold = 0.02; 
    pd.itr = 3000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'os2_more_itr'
    % optimal stimulation 2nd pass: sparse
    pd.tau = 1; %2?
    pd.delta = 0.1; %?
    pd.threshold = 0.1; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 2000;
    % soft-thresholding
    pd.thr_func = 'soft';
  case 'dynamic'
    % use parameters in rozell2008
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.1; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 50;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_hard'
    % use parameters in rozell2008
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.1; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 50;
    % soft-thresholding
    pd.thr_func = 'hard'; 
  case 'dynamic_0p05_hard'
    % use parameters in rozell2008
    pd.tau = 0.1;
    pd.delta = 0.01; %?
    pd.threshold = 0.05; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 100;
    % soft-thresholding
    pd.thr_func = 'hard'; 
  case 'dynamic_0p05'
    pd.tau = 0.1;
    pd.delta = 0.01; %?
    pd.threshold = 0.05; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 50;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_25'
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.1;
    pd.itr = 25;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_125_small_delta'
    pd.tau = 10;
    pd.delta = 0.2; %?
    pd.threshold = 0.1;
    pd.itr = 125;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_2500_small_delta'
    pd.tau = 10;
    pd.delta = 0.01; %?
    pd.threshold = 0.1;
    pd.itr = 2500;
    pd.thr_func = 'soft'; 
  case 'dynamic_50'
    % use parameters in rozell2008
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.1; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 50;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_100'
    % use parameters in rozell2008
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.1; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 100;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_200'
    % use parameters in rozell2008
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.1; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 200;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_thr_0p0005'
    % use parameters in rozell2008
    pd.tau = 10;
    pd.delta = 0.25; %?
    pd.threshold = 0.0005; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 800;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_s'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.5; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 50;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_s_25'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.5; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 25;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_s_200'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.5; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 200;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_s_300'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.5; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 300;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_s_1000'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.5; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 1000;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_s_small_delta'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 0.1; %?
    pd.threshold = 0.5; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 50;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_ns'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 0.5; %?
    pd.threshold = 0.005; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 100;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_ns_400'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 0.5; %?
    pd.threshold = 0.005; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 400;
    % soft-thresholding
    pd.thr_func = 'soft'; 
  case 'dynamic_ns_600'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 0.5; %?
    pd.threshold = 0.005; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 600;
    % soft-thresholding
    pd.thr_func = 'soft'; 
case 'dynamic_ns_700'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 0.5; %?
    pd.threshold = 0.005; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 700;
    % soft-thresholding
    pd.thr_func = 'soft'; 

case 'dynamic_ns_800'
    % dynamic sparse
    pd.tau = 10;
    pd.delta = 0.5; %?
    pd.threshold = 0.005; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 800;
    % soft-thresholding
    pd.thr_func = 'soft'; 

  case 'dynamic_25_block'
    % use parameters in rozell2008
    pd.tau = 10;
    pd.delta = 1; %?
    pd.threshold = 0.1; %0.0005 default; changed to 0.05
                          %later. Current value for whitened images.
    pd.itr = 25;
    % soft-thresholding
    pd.thr_func = 'block'; 

end

% $$$ % Normalize the gray level for display
% $$$ ph.graylevel = [-2,2];
