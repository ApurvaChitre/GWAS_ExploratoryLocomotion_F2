#!/bin/bash

#PBS -S /bin/bash
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=3
#PBS -j oe
#PBS -q condo

export R_LIBS=/home/aschitre/R_libs:$R_LIBS
module load R



Rscript conditional_analysis.R $PBS_O_WORKDIR $expname $trait $chr


