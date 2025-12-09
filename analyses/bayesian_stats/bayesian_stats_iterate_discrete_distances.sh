#!/bin/bash -l
#$ -P npbssmic
#$ -l h_rt=72:00:00
#$ -N bayesian
#$ -j y
#$ -m beas
# Compute resource requirements
#$ -pe omp 16

# Load the miniconda module
module load miniconda
# Initialize conda in shell
mamba activate pymc_spline

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

# Run python script and pass distance variable
python bayesian_stats_iterate_discrete_distances.py
