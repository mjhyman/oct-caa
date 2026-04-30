#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=240:00:00

# Name of job
#$ -N histo_gam

# Entire node
#$ -pe omp 16
#$ -l mem_per_core=16G

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

# Combine output/error files into single file
#$ -j y

module load matlab/2025a
matlab -nodisplay -nojvm -nodesktop -r "histology_swp_stain_gam; exit"
