# Monte Carlo Sampling Script

**Instructions for running the Monte Carlo sampling script.**

## Steps

1. **Install Miniconda locally**  
   If working locally, [download Miniconda](https://docs.conda.io/en/latest/miniconda.html) and follow installation instructions.

   **On the SCC:** Load the miniconda module:
   ```bash
   module load miniconda
   ```
2. Open the console and navigate to the "analyses" subfolder
  ```bash
  cd [path_to_repository]/analyses
  ```
3. Build the anaconda environment
   ```bash
   conda env create -f pymc_spline.yml
   ```
4. Activate the environment
```bash
conda activate pymc_spline
```

5. Navigate to the Python script (*/analyses/bayesian_stats_distance_model_scc.py) and modify the respective directory paths to point at the spreadsheet.
6. If running this from the Python IDE, then simply run the script.
7. If submitting this as a batch script on the SCC, then submit a batch job
```bash
qsub [path_to_repository]/analyses/bayesian_stats_distance_model_scc.sh
```
