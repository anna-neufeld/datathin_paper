#!/bin/bash
#$ -cwd

Rscript run_gamma_eps_big.R --simname $1 --nreps $2
