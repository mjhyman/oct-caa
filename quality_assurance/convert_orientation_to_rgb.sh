#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=48:00:00

# Name of job
#$ -N ori2rgb

# Run this for CAA22 front (dense EPVS)
#$ -pe omp 16
#$ -l mem_per_core=8G

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

# Combine output/error files into single file
#$ -j y

module load matlab/2024b
matlab -nodisplay -r 'convert_orientation_to_rgb; exit'