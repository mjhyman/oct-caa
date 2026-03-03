#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=120:00:00

# Name of job
#$ -N swp_vox

### set the task IDs for CAA25o
#$ -t 7

# Entire Node w/ 1TB
# -pe omp 36

# Compute Settings (entire node w/ 512 GB)
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

module load matlab/2025a
matlab -nodisplay -batch epvs_density_batch $SGE_TASK_ID

