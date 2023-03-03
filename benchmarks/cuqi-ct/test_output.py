#!/usr/bin/env python3
""" Script to regression test the output of the deconvolution-1D benchmark. """

import argparse
import umbridge
import pytest
import cuqi
import numpy as np
from scipy.io import loadmat

# %% CUQI imports
from cuqi.model import LinearModel
from cuqi.geometry import Image2D
from cuqi.array import CUQIarray
from cuqi.distribution import Gaussian, GMRF, Laplace_diff, Cauchy_diff
from cuqi.problem import BayesianProblem


# Set up CUQIpy testproblem to compare with
print('Setting up problem')
data = np.load("data_ct.npz")
#TP = cuqi.testproblem.Deconvolution1D(
#    dim=128,
#    PSF="gauss",
#    PSF_param=5,
#    phantom=data["exact"], # Phantom was "pc".
#    noise_std=0.01,
#    noise_type="gaussian"
#)

N = 256
nv = 30

# Load the matrix for the forward model and scale by pixel width
Amat = loadmat('A' + str(N) + '_' + str(nv) + '.mat')["A"]

lower = -1.0
upper = 1.0
width = upper - lower
dx = width/N

Amat = dx*Amat

# Create the CUQIpy linear model
A = LinearModel(Amat)

# Create the visual_only Image2D geometries of image and sinogram spaces
dg = Image2D(( N,N), visual_only=True)
rg = Image2D((nv,N), visual_only=True)

# Equip linear operator with geometries
A.domain_geometry = dg
A.range_geometry = rg

# Create the CUQI data structure from vectorized image and geometry
imC = CUQIarray(data["exact"], geometry=dg)

# Specify placeholder x distribution
x = Gaussian(mean=np.zeros(A.domain_dim), 
                           cov=0.01,
                           geometry=A.domain_geometry)
       
# Choose noise std and create data distribution
s = 0.01
y = Gaussian(A@x, s**2)

# Create CUQIarray with loaded noisy sinogram data.
y_data = CUQIarray(data["y_data"], geometry=rg)

# Set up the Bayesian problem with random variable and observed data.
BP = BayesianProblem(y, x).set_data(y=y_data)

parameters = imC # Extract test parameters to evaluate on

# Parse command line arguments
parser = argparse.ArgumentParser(description='Model output test.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

# Define models for each type of prior
print('Setting up models')
model_Gaussian = umbridge.HTTPModel(args.url, "CT_Gaussian")
model_GMRF = umbridge.HTTPModel(args.url, "CT_GMRF")
model_CMRF = umbridge.HTTPModel(args.url, "CT_CMRF")
model_LMRF = umbridge.HTTPModel(args.url, "CT_LMRF")

# Model for exact solution
model_exactSolution = umbridge.HTTPModel(args.url, "CT_ExactSolution")

# Check basic model properties
print('Checking model properties')
assert model_Gaussian.get_input_sizes() == [256**2]
assert model_GMRF.get_input_sizes() == [256**2]
assert model_CMRF.get_input_sizes() == [256**2]
assert model_LMRF.get_input_sizes() == [256**2]
assert model_exactSolution.get_input_sizes() == [0]

assert model_Gaussian.get_output_sizes() == [1]
assert model_GMRF.get_output_sizes() == [1]
assert model_CMRF.get_output_sizes() == [1]
assert model_LMRF.get_output_sizes() == [1]
assert model_exactSolution.get_output_sizes() == [256**2]

# Check that the models support evaluation and gradient
print('Checking model capabilities')
assert model_Gaussian.supports_evaluate()
assert model_GMRF.supports_evaluate()
assert model_CMRF.supports_evaluate()
assert model_LMRF.supports_evaluate()
assert model_exactSolution.supports_evaluate()

assert model_Gaussian.supports_gradient()
assert model_GMRF.supports_gradient()
assert model_CMRF.supports_gradient()
assert not model_LMRF.supports_gradient()
assert not model_exactSolution.supports_gradient()

# Evaluate models and compare with testproblem
print('Evaluating models')
output_Gaussian = model_Gaussian([parameters.tolist()])[0][0]
output_GMRF = model_GMRF([parameters.tolist()])[0][0]
output_CMRF = model_CMRF([parameters.tolist()])[0][0]
output_LMRF = model_LMRF([parameters.tolist()])[0][0]
output_exactSolution = model_exactSolution([[]])[0]

# Check that the output is correct
print('Checking model output')
BP.prior = cuqi.distribution.Gaussian(np.zeros(256**2), 0.01,
                              geometry=BP.likelihood.geometry)
assert output_Gaussian == pytest.approx(BP.posterior.logpdf(parameters))

BP.prior = cuqi.distribution.GMRF(np.zeros(256**2), 1/(0.01),
                                  physical_dim=2,
                                  geometry=BP.likelihood.geometry)
assert output_GMRF == pytest.approx(BP.posterior.logpdf(parameters))

BP.prior = cuqi.distribution.Cauchy_diff(np.zeros(256**2), 0.01,
                                         physical_dim=2,
                                         geometry=BP.likelihood.geometry)
assert output_CMRF == pytest.approx(BP.posterior.logpdf(parameters))

BP.prior = cuqi.distribution.Laplace_diff(np.zeros(256**2), 0.1,
                                          physical_dim=2,
                                          geometry=BP.likelihood.geometry)
assert output_LMRF == pytest.approx(BP.posterior.logpdf(parameters))

assert np.allclose(output_exactSolution, parameters)

print('Regression testing against known norm values of output')
assert np.linalg.norm(output_Gaussian) == pytest.approx(-346295.2714778)
assert np.linalg.norm(output_GMRF) == pytest.approx(139120.2419227)
assert np.linalg.norm(output_CMRF) == pytest.approx(828.4505233852567)
assert np.linalg.norm(output_LMRF) == pytest.approx(466477.36087666)
assert np.linalg.norm(output_exactSolution) == pytest.approx(96.06430138193903)

print('All tests passed')

