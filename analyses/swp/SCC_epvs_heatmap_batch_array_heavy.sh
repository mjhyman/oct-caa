#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=120:00:00

# Name of job
#$ -N swp_p2

# set the task ID
#$ -t 8-9

# Run this for CAA22 front (dense EPVS)
#$ -pe omp 28
#$ -l mem_per_core=18G

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="
echo "Starting task number $SGE_TASK_ID"

# Combine output/error files into single file
#$ -j y

module load matlab/2024b
matlab -nodisplay -batch epvs_density_batch $SGE_TASK_ID

