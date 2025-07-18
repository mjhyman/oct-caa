#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Specify number of cores
#$ -pe omp 18

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=120:00:00

# Name of job
#$ -N caa25f

# Combine output/error files into single file
#$ -j y

module load matlab/2024b
matlab -nodisplay -nosplash -r "run('epvs_density_caa25f.m'); exit;"

