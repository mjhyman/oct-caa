#!/bin/bash -l
#$ -P npbssmic
#$ -l h_rt=48:00:00
#$ -N nnunet
#$ -j y
#$ -m beas
# Compute resource requirements
#$ -pe omp 4
#$ -l gpus=1
# The installed Pytorch version requires a minimum of version 7.0
#$ -l gpu_c=7.0

# Load the miniconda module
module load miniconda
# Initialize conda in shell
mamba activate nnunet_pytorch_129

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

# Run python script and pass distance variable
nnUNetv2_train 001 3d_fullres all
