#!/bin/bash
#PBS -l select=1:ncpus=1:mem=40gb
#PBS -l walltime=05:00:00
#PBS -o output.file
#---------------------------------------------
date

module add lang/r/4.0.3-bioconductor-gcc

cd $PBS_O_WORKDIR

Rscript VZ_summary_mvMR_BF_function.R

Rscript VZ_summary_mvMR_SSS_function.R

Rscript Kett_SA.R

date
