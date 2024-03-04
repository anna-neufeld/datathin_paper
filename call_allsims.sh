#!/bin/bash
#$ -cwd

Rscript run_allsims.R --simname $1 --nreps $2
