#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=72:00:00

# Name of job
#$ -N hswp_batch

# Launch a job for each stain to run in parallel (3 stains in total)
#$ -t 1-3

# Entire node w/ 512 GB RAM
#$ -pe omp 28
#$ -l mem_per_core=9G

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

# Combine output/error files into single file
#$ -j y

module load matlab/2024b
matlab -nodisplay -nojvm -nodesktop -batch histology_size_weighted_proximity $SGE_TASK_ID
