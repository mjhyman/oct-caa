import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Spline
from patsy import dmatrix
import statsmodels.api as sm

# Gaussian Process
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel

# Set seed and generate data
np.random.seed(0)
x = np.sort(np.random.rand(50))  # 50 points in [0, 1]
y = np.sin(2 * np.pi * x) + np.random.normal(0, 0.1, size=x.shape)

x_pred = np.linspace(0, 1, 200)

# ------------------------------------------
# 1. Cubic Splines using patsy + statsmodels
# ------------------------------------------
spline_basis = dmatrix("bs(x, df=6, degree=3, include_intercept=True)", {"x": x}, return_type='dataframe')
spline_model = sm.OLS(y, spline_basis).fit()

spline_pred_basis = dmatrix("bs(x_pred, df=6, degree=3, include_intercept=True)", {"x_pred": x_pred}, return_type='dataframe')
spline_pred = spline_model.predict(spline_pred_basis)

# ------------------------------------------
# 2. Gaussian Process Regression
# ------------------------------------------
kernel = ConstantKernel(1.0) * RBF(length_scale=0.2) + WhiteKernel(noise_level=0.01)
gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=0.0)
gp.fit(x.reshape(-1, 1), y)

y_gp_pred, y_gp_std = gp.predict(x_pred.reshape(-1, 1), return_std=True)

# ------------------------------------------
# Plotting
# ------------------------------------------
plt.figure(figsize=(12, 6))
plt.scatter(x, y, color='black', label='Data', zorder=3)
plt.plot(x_pred, np.sin(2 * np.pi * x_pred), '--', color='gray', label='True Function')

# Spline plot
plt.plot(x_pred, spline_pred, color='blue', label='Cubic Spline Fit', lw=2)

# GP plot with uncertainty
plt.plot(x_pred, y_gp_pred, color='darkorange', label='GP Mean Prediction', lw=2)
plt.fill_between(x_pred, y_gp_pred - 2*y_gp_std, y_gp_pred + 2*y_gp_std,
                 color='orange', alpha=0.3, label='GP 95% CI')

plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Spline vs Gaussian Process Regression")
plt.grid(True)
plt.tight_layout()
plt.show()
