#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=120:00:00

# Name of job
#$ -N epvs_cluster

# set the task ID
#$ -t 1-9

# Compute Settings
#$ -pe omp 16
#$ -l mem_per_core=16G

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
matlab -nodisplay -batch epvs_clustering_main $SGE_TASK_ID

