#!/bin/bash -l

# SCC Project
#$ -P npbssmic

# Number of cores
#$ -pe omp 8                 # Stan runs 4 chains -> 4 cores

# Memory per core
#$ -l mem_per_core=18G

# Time limit for job
#$ -l h_rt=120:00:00

# Job name 
#$ -N bayes_histo
#$ -j y
#$ -m ea

echo "Start: $(date)  Host: $(hostname)  Slots: $NSLOTS"

module load R/4.5.2
export RENV_ACTIVATE_PROJECT=FALSE
Rscript bayesian_stats_fit_parenchyma_shells_brms_histology.R

echo "End: $(date)"
