# Overview
This repository includes the code used in Degen et al., 2026 for running H3K27 methylation simulations and parameter sweeps.

Instructions on how to reproduce the modeling-related results presented in the main text and Modeling Supplement are included in this document.

Please contact degeneleanor@gmail.com or shelby.blythe@gmail.com with any questions or comments.

# Set up
MATLAB R2022a was used in the writing and running of all scripts. All scripts can be run after the installation of MATLAB. No auxiliary MATLAB toolboxes are required.

# H3K27me simulation functions
Scripts that run simulations call the functions in `~/simulationFunctions/`. You will need to add this folder to your MATLAB path. In this folder, **`simulateK27.m`** takes E(z) and Pho live imaging measurements and a structure of parameter values as inputs and simulates E(z) activity along a nucleosome array. This function calls **`methylationFunction.m`** where the decision about whether a methylation reaction occurs is made.

# Generating figures

## Fig. 2C & C'
These plots present the maternal and paternal H3K27 methylation stimulation results for the primary model discussed in the main text. To generate them, you will need to run the following scripts:

- **`simulate_K27me_final.m`** loads the E(z) and Pho concentration estimates required for running the simulations, runs maternal and paternal H3K27 methylation simulations by calling **`simulateK27.m`**, and then plots and saves the results.
- **`plot_100xSimResults.m`** loads the results produced by **`simulate_K27me3_final.m`** in the first code block (right now, the script loads the simulation results stored in `~/22-Mar-2026 18:40:31_mat/` and `~/22-Mar-2026 18:41:58_pat/`, the ones used for the paper but can be recreated), and in the second code block, creates the plots of Fig. 2C & C'. Note: the third code chunk is related to Fig. 2D, as discussed below.

## Fig. 2D, Fig. SM7 A, Fig. SM7 B
These plots show the maternal vs. paternal contributions to H3K27me1 (Fig. SM7 A of Modeling Supplement), H3K27me2 (Fig. SM7 B of Modeling Supplement), and H3K27me3 (Fig. 2D of main text). 

To generate them, run the first code chunk of **`plot_100xSimResults.m`** to load your simulation results (produced by **`simulate_K27me_final.m`** or included in in `~/22-Mar-2026 18:40:31_mat/` and `~/22-Mar-2026 18:41:58_pat/` as discussed above). Run the third code block to calculate and plot the maternal vs. paternal contributions to H3K27 methylation.

## Fig. S1
This figure shows the results of H3K27 methylation simulations that incorporate non-allosteric, nucleation-independent E(z) activity. To produce the plotted simulations, run the first code block of **`nonAllo_nucIndep_model.m`** to load the requisite E(z) and Pho concentration estimates, and then the second code block to run the simulation and plots the results. Right now, the script is configured to run the simulation with the initial conditions of the paternal chromosomes. If you'd like to run a maternal simulation, comment out line 30, and uncomment line 29.

## Fig. SM1 A, Fig. SM1 B
These plots show the results of evaluating the model's output across different initial conditions (varying H3K27me2/3), with 0% of histones starting with me1 (Fig. SM1 A) or 25% of histones starting with me1 (Fig. SM1 B). To run the initial condition sweep and plot the results, run **`IC_sweep.m`**. Right now the code is configured to perform the sweep with 0% of histones starting with me1. If you'd like to change the me1 condition, change line 33 to a value between 0 and 1. If you want to bypass running the sweep and just plot the results, run the first code block of **`IC_sweep.m`** to define the parameters, load either `allOutputs_IC_0me1.mat` or `allOutputs_IC_0.25me1.mat`, and then run the third code block of the script.

## Fig. SM3 A, Fig. SM3 B
These plots report the process of fitting the model's baseline methylation probability, V, to H3K27me1 (Fig. SM3 A of Modeling Supplement) and H3K27me3 (Fig. SM3 B of Modeling Supplement) measurements. To create the plots, first a parameter sweep for V was performed. The results of this sweep are saved in `~/V_sweep/V_sweep_final_19-Mar-2026 12:06:47/`, and can be recreated by running **`V_sweep.m`**. If running the parameter sweep, specify the path to the directory where you'd like to save the results. To plot the comparison between V sweep simulation results and experimental measurements (saved in `~/ChIP_data_for_fitting/ChIPseq_NC14_maxes_3.mat` and `~/ChIP_data_for_fitting/ChIPseq_NC14_mins_1.mat'`), run **`plotSweep_V.m`**.

## Fig. SM4
This figure shows the results of a parameter sweep for the model's allosteric stimulation parameters. You can perform the parameter sweep by running the first and second code blocks of **`allo_sweep.m`**, and plot the results by running the third code block. Note: when performing the sweep, make sure V is set at the desired value in line 25. The results of the allostery sweeps when V = 7.35 x 10<sup>-4</sup> (`allOutputs_allo.mat`) and V = 14.7 x 10<sup>-4</sup> (`allOutputs_allo_doubleV.mat`) are saved in `~/Kd_sweep/`. If you'd like to bypass performing the sweep, run the first code block of **`allo_sweep.m`** to define the parameters, load `allOutputs_allo.mat` or `allOutputs_allo_doubleV.mat`, and then run the third code block of the script. 

## Fig. SM5
This figure shows the results of a parameter sweep for the model's dissociation constants. You can perform the parameter sweep by running the first and second code blocks of **`Kd_sweep.m`**, and plot the results by running the third code block. Note: when performing the sweep, make sure V is set at the desired value in line 25. The results of the Kd sweeps when V = 7.35 x 10<sup>-4</sup> (`allOutputs_Kd.mat`) and V = 14.7 x 10<sup>-4</sup> (`allOutputs_Kd_doubleV.mat`) are saved in `~/Kd_sweep/`. If you'd like to bypass performing the sweep, run the first code block of **`Kd_sweep.m`** to define the parameters, load `allOutputs_Kd.mat` or `allOutputs_Kd_doubleV.mat`, and then run the third code block of the script.



