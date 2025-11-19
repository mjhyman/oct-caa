To run the monte carlo sampling script, perform the following steps:

1) If running locally, install miniconda. If running on the SCC, load the miniconda module (module load miniconda)
2) Open the console and navigate to the "analyses" subfolder
3) Build the anaconda environment (conda env create -f pymc_spline.yml)
4) Activate the environment (conda activate pymc_spline)
5) Navigate to the Python script (*/analyses/bayesian_stats_distance_model_scc.py) and modify the respective directory paths to point at the spreadsheet.
6) If running this from the Python IDE, then simply run the script.
7) If submitting this as a batch script on the SCC, then submit a batch job (qsub [path_to_repository]/analyses/bayesian_stats_distance_model_scc.sh)
