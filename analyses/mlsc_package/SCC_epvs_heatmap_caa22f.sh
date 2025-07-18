#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Specify number of cores
#$ -pe omp 28
#$ -l mem_per_core=13G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=120:00:00

# Name of job
#$ -N caa22f

# Combine output/error files into single file
#$ -j y

module load matlab/2024b
matlab -nodisplay -nosplash -r "epvs_density_caa22f; exit;"

