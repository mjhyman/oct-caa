#!/bin/bash -l
#$ -P npbssmic
#$ -t 1-12
#$ -l h_rt=72:00:00
#$ -N bayes_array
#$ -j y
#$ -m beas
# Compute resource requirements
#$ -pe omp 8
#$ -l mem_per_core=8G

# Load the miniconda module
module load miniconda
# Initialize conda in shell
mamba activate pymc_spline

# Define array of distances (should match --array range length)
DISTANCES=(40 80 120 160 200 240 280 320 360 400 440 480)

# Subtract 1 from the task ID to make it zero-indexed
# since the cluster requires a one-indexed number
TASK_ID_ZERO=$((SGE_TASK_ID - 1))
DISTANCE=${DISTANCES[$TASK_ID_ZERO]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "Distance = $DISTANCE (microns)"
echo "=========================================================="

# Run python script and pass distance variable
python bayesian_stats_lognormal_at_each_distance.py --distance $DISTANCE