#!/bin/bash

njobs=40

qsub -q w-bigmem.q -e iotrash/ -o iotrash/ -l h='biostat-b12|biostat-b13|biostat-b14|biostat-b18' -l h_vmem=10G -t 1-$njobs -tc 50 ./call_sim_binomial.sh $1 $2
