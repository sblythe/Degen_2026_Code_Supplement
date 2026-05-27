# Overview
This repository includes the code used in Degen et al., 2026 for running H3K27 methylation simulations and parameter sweeps.

Instructions on how to reproduce the modeling-related results presented in the main text and Modeling Supplement are included in this document.

Please contact degeneleanor@gmail.com or shelby.blythe@northwestern.edu with any questions or comments.

# Set up
MATLAB R2022a was used in the writing and running of all scripts. All scripts can be run after the installation of MATLAB. No auxiliary MATLAB toolboxes are required.

# H3K27me simulation functions
Scripts that run simulations call the functions in `~/simulationFunctions/`. You will need to add this folder to your MATLAB path. In this folder, **`simulateK27.m`** takes E(z) and Pho live imaging measurements and a structure of parameter values as inputs and simulates E(z) activity along a nucleosome array. This function calls **`methylationFunction.m`** where the decision about whether a methylation reaction occurs is made.

# Generating figures

## Figure 2C & C'
These plots present the maternal and paternal H3K27 methylation stimulation results for the primary model discussed in the main text. To generate them, you will need to run the following scripts:

- **`simulate_K27me_final.m`** loads the E(z) and Pho concentration estimates required for running the simulations, runs maternal and paternal H3K27 methylation simulations by calling **`simulateK27.m`**, and then plots and saves the simulation results.
- **`plot_100xSimResults.m`** loads the simulation results produced by **`simulate_K27me3_final.m`** in the first code chunk (right now, the script loads the simulation results in `~/22-Mar-2026 18:40:31_mat/` and `~/22-Mar-2026 18:41:58_pat/`, the ones used for the paper but can be recreated), and in the second code chunk, creates the plots of Fig. 2C & C'. Note: the third code chunk is related to Fig. 2D, as discussed below. 
