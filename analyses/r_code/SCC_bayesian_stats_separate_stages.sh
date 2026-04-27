#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Request a whole node with 28 cores and at least 384 GB of RAM.
# Specify number of cores
#$ -pe omp 8
# Specify memory per core
#$ -l mem_per_core=18G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=120:00:00

# Name of job
#$ -N rbayes

# Combine output/error files into single file
#$ -j y

module load R
Rscript bayesian_stats_separate_stages.R

