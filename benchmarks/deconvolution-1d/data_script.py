# %% Based on CUQIpy version v0.3.0
import numpy as np
import cuqi
import matplotlib.pyplot as plt

# Check cuqi version
assert cuqi.__version__ == "0.3.0"
print ("CUQIpy version: ", cuqi.__version__)

# %% Benchmark data generation

np.random.seed = 0

TP = cuqi.testproblem.Deconvolution1D(
    dim=128,
    PSF="gauss",
    PSF_param=5,
    phantom="square",
    noise_std=0.01,
    noise_type="gaussian",
)

# Save data code (uncomment to overwrite data file)
#np.savez("data_square.npz", data=TP.data, exact=TP.exactSolution)

# %% Plot data

TP.data.plot()
plt.title("Data")

# %% Plot exact solution

TP.exactSolution.plot()
plt.title("Exact solution")

# %% Plot posterior samples for Gaussian

TP.prior = cuqi.distribution.Gaussian(np.zeros(128), 0.01, name="x")
samples = TP.UQ()
plt.title("Gaussian with delta=0.01")

# %% Plot posterior samples for GMRF

TP.prior = cuqi.distribution.GMRF(np.zeros(128), 1/0.01, name="x")
samples = TP.UQ()
plt.title("GMRF with delta=0.01")

# %% plot posterior samples for LMRF

TP.prior = cuqi.distribution.Laplace_diff(np.zeros(128), 0.01, name="x")
samples = TP.UQ()
plt.title("LMRF with delta=0.01")

# %% Plot posterior samples for CMRF

TP.prior = cuqi.distribution.Cauchy_diff(np.zeros(128), 0.01, name="x")
samples = TP.UQ()
plt.title("CMRF with delta=0.01")


# %%
