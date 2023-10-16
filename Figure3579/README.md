The file "make_plots.R" reads in saved results from the folder "resAD" and generates figures 3,5,7, and 9 in the paper.

The code to generate the results in the folder "resAD" are found in the files "run_gamma_eps_big.R", "run_gamma_eps_small.R", and "run_binomial_eps.R". Helper functions for these files are stored in "GCS_evalfuns.R" and "GCS_splitfuns.R". 

All other files (the .sh) files were used to run the simulations on a cluster, and may need to be updated to run the code locally or on different cluster systems. 