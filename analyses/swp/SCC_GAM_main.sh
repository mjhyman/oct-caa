#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=120:00:00

# Name of job
#$ -N gam_oct

# Compute Settings (entire node w/ 512 GB)
#$ -pe omp 36

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="
echo "Starting task number $SGE_TASK_ID"

# Combine output/error files into single file
#$ -j y

module load matlab/2025a
matlab -nodisplay -r 'swp_mus_ret_gam_main; exit' 

