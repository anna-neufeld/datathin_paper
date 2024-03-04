#!/bin/bash
#$ -cwd

Rscript run_binomial_eps.R --simname $1 --nreps $2
