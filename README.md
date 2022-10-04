# Inhibitory-interneuron-models-via-convex-optimization
Code to produce the results in the paper "Modeling biologically realistic inhibitory interneurons in sensory coding models" in PLoS Computational Biology (2015) by Zhu and Rozell. Code by M. Zhu.

Follow these steps to reproduce results from the paper "Modeling inhibitory interneurons in efficient sensory coding models", Mengchen Zhu and Christopher Rozell, PLoS CB 2015.

1. Decompose the recurrent matrix
   run_arpca.m
   The adaptive robust PCA algorithm is stochastic and can produce different solutions each time you run it. It does invariably find a "good" solution (low rank and sparse) in a couple of adaptive iterations. To reproduce the results in the paper exactly, load the saved decomposition rpca_rw_alm_gamma_rpca_250_0_250_tol_0p01_alpha_2p5_beta_0p01.mat.

2. interneuron_sim.m
   This is the main file. The recurrent matrices are first decomposed into excitatory and inhibitory components (interneuron_decompose.m), the performance of the interneuron network is measured in rl2e_sparsity.m, and the RFs of the interneurons are mapped in map_RF.m.
