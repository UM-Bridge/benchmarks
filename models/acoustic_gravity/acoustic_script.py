import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.stats import multivariate_normal
from scipy.spatial.distance import cdist
from scipy.signal import hilbert

import pandas as pd
import tinyDA as tda
import umbridge

# Set domain size and simulation time: this info is shared between the forward and parameter models 
Lx = 20
Lz = 1 
TEnd = 10 
H = 1.5

# Each element of the list "captors" is a tuple of (x,z) coordinates. 
# Coordinates are given in the (scaled) domain: x \in [0, grid_size_x], z \in [0, grid_size_y]
sensors = [(4,0.8), (13.5, 0.6)]

# Create sorted sensor list to pass to the forward model: do not change this part 
sensors_sort = np.array(sensors, dtype=[('x',float),('y',float)])
sensors_sort = np.sort(sensors_sort, order=['x','y'])
sensors = sensors_sort.tolist()
n_sensors = len(sensors)
print(n_sensors)

# connect to the UM-Bridge model.
umbridge_model = umbridge.HTTPModel('http://localhost:4242', "forward")

# wrap the UM-Bridge model in the tinyDA UM-Bridge interface.
config={"captors":sensors}

# my_model = tda.UmBridgeModel(umbridge_model) # DG: for now without config 
my_model = tda.UmBridgeModel(umbridge_model, umbridge_config=config)

# For now the DG file does not use the arg config 
nx = umbridge_model.get_input_sizes(config)[0] #the input is the value of f(x) at each grid point
ny = umbridge_model.get_output_sizes(config)[0] #the outout is a time series of the pressure

print(f"input:{nx}, output:{ny}")

# Generate synthetic data
exact = np.zeros(nx)
# True source in the middle of the domain 
idx_middle = int(nx/2)
source_width = 4 # in length unit 
idx_width = int(0.5 * source_width / Lx * nx )
exact[idx_middle - idx_width : idx_middle + idx_width]=1
#print(exact)
d_true = my_model(exact)

# add some noise to the model output
sigma_noise = 0.01
d = d_true + np.random.normal(loc=0, scale=sigma_noise, size=ny)

# GRF parameters
length_scale = 1
decay_rate = 10

# Level-set parameter
eps = 1e-6

def rbf_covariance(grid, length_scale, variance=1.0):
    dists = cdist(grid, grid, 'euclidean')
    return variance * np.exp(-0.5 * (dists / length_scale) ** 2)

# Use this covariance structure moving forward to match theory
def matern32_covariance(grid, length_scale, variance=1.0):
    dists = cdist(grid, grid, 'euclidean')
    sqrt3 = np.sqrt(3.0)
    scaled_dists = sqrt3 * dists / length_scale
    return variance * (1.0 + scaled_dists) * np.exp(-scaled_dists)

def restriction(sampin):
    field_cut = sampin
    return (field_cut > eps).astype(int)

dx = Lx/nx
x = np.linspace(0, Lx, nx).reshape(-1, 1)
cov = matern32_covariance(x, length_scale)
cov += eps * np.eye(nx)
mean = np.zeros(nx)
my_prior = multivariate_normal(mean=mean, cov=cov)

def sigma_linear(v, b):
    return v + b

def sigma_signsensitive(v, b):
    pos = v >= 0
    raw = np.empty_like(v)
    raw[pos] = v[pos] + 1.0 / b
    raw[~pos] = (1.0 / b) * np.exp(b * v[~pos])
    return raw

def sigma_softplus(v, b):
    return np.logaddexp(0.0, b * v)

def sigma_envelope_exp(v, b):
    env = np.abs(hilbert(v))
    return np.exp(b * env)

class WassersteinLoglike:
    """Wasserstein log-likelihood.
    `sigma` is a scaling function"""

    def __init__(self, data, sigma, lam=1.0, b=1.0, n_quantiles=200):
        self.data = np.asarray(data, dtype=float).ravel()
        self.sigma = sigma
        self.lam = lam
        self.b = b
        self.n_quantiles = n_quantiles

    def _density(self, x):
        raw = self.sigma(x, self.b)
        return raw / raw.sum()

    def _quantile_function(self, density):
        cdf = np.cumsum(density); cdf = cdf / cdf[-1]
        t = np.arange(len(density), dtype=float)
        u = (np.arange(self.n_quantiles) + 0.5) / self.n_quantiles
        return np.interp(u, cdf, t)

    def loglike(self, x):
        x = np.asarray(x, dtype=float).ravel()
        p, q = self._density(x), self._density(self.data)
        Fx, Fd = self._quantile_function(p), self._quantile_function(q)
        return -self.lam * np.sum((Fx - Fd) ** 2) / self.n_quantiles

class MultiSensorLoglike:
    """likelihoods: a list of length n_sensors likelihoods, each already
    constructed with its own sensor's reference data."""
    def __init__(self, likelihoods, n_sensors):
        self.likelihoods = likelihoods
        self.n_sensors = n_sensors

    def loglike(self, x):
        x = np.asarray(x, dtype=float).ravel()   # <- use the argument, not a stored copy
        return sum(like.loglike(x[i::self.n_sensors]) for i, like in enumerate(self.likelihoods))
    
# Test for Gaussian log like
sigma = 2.0
cov_likelihood = sigma**2*np.eye(d_true.shape[0])
my_loglike_gaussian = tda.GaussianLogLike(d_true, cov_likelihood)

# Test for Wasserstain based log like
# TODO move implementation into tinyDA
my_loglike = MultiSensorLoglike(likelihoods=[
    WassersteinLoglike(d_true[0::n_sensors], sigma_signsensitive, lam=0.001, b=0.01),
    WassersteinLoglike(d_true[1::n_sensors], sigma_signsensitive, lam=0.001, b=0.01),
], n_sensors=2)

my_loglike_linear = MultiSensorLoglike(likelihoods=[
    WassersteinLoglike(d_true[0::n_sensors], sigma_linear, lam=1.0, b=0.1),
    WassersteinLoglike(d_true[1::n_sensors], sigma_linear, lam=1.0, b=1.1),
], n_sensors=2)

my_loglike_envelope = MultiSensorLoglike(likelihoods=[
    WassersteinLoglike(d_true[0::n_sensors], sigma_envelope_exp, lam=35.22, b=5.0),
    WassersteinLoglike(d_true[1::n_sensors], sigma_envelope_exp, lam=72.03, b=5.0),
], n_sensors=2)

def levelset_model(parameters):
    levelset_params = restriction(parameters)
    return my_model(levelset_params)

my_posterior = tda.Posterior(my_prior, my_loglike, levelset_model)

class CrankNicolson(tda.CrankNicolson):
    def __init__(self, scaling, adaptive):
        super().__init__(scaling=scaling, adaptive=adaptive)

    def make_proposal(self, link):
        # make a pCN proposal.
        self.scaling = min(self.scaling, 1.0 - 1e-3)
        return np.sqrt(
            1 - self.scaling**2
        ) * link.parameters + self.scaling * np.random.multivariate_normal(
            self._mean, self.C
        )

# preconditioned Crank-Nicolson
pcn_scaling = 0.15
pcn_adaptive = True
my_proposal = CrankNicolson(scaling=pcn_scaling, adaptive=pcn_adaptive)

# For testing purposes, iteration number is small for the given problem; Choose a larger number for real applications.
my_chains = tda.sample(my_posterior, my_proposal, iterations=2, n_chains=2, force_sequential=True)

import arviz as az

burnin = 0 
idata = tda.to_inference_data(my_chains, burnin=burnin)
az.to_netcdf(idata, "results_dg_10000.nc") # to store







