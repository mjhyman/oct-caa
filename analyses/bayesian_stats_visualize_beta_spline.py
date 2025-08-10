import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import os
from patsy import dmatrix

# -----------------------
# Configuration
# -----------------------
sheet_names = ["scattering", "retardance", "orientation"]

posterior_dir = ("/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/statistics/"
                 "spline_model_outputs/scattering_frontal")

# Spline Basis
fname = "spline_basis_scattering_frontal.pkl"
with open(os.path.join(posterior_dir,fname), "rb") as f:
    spline_basis = pickle.load(f)

# Summary
fname = "summary_scattering_frontal.pkl"
with open(os.path.join(posterior_dir,fname), "rb") as f:
    summary = pickle.load(f)

# Trace
fname = "trace_scattering_frontal.pkl"
with open(os.path.join(posterior_dir,fname), "rb") as f:
    trace_frontal = pickle.load(f)

# XSpline
fname = "xspline_scattering_frontal.pkl"
with open(os.path.join(posterior_dir,fname), "rb") as f:
    xspline = pickle.load(f)
