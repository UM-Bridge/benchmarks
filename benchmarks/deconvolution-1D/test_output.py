#!/usr/bin/env python3
""" Script to regression test the output of the deconvolution-1D benchmark. """

import argparse
import umbridge
import pytest
import cuqi
import numpy as np

# Set up CUQIpy testproblem to compare with
print('Setting up testproblem')
data = np.load("data_square.npz")
TP = cuqi.testproblem.Deconvolution1D(
    dim=128,
    kernel="gauss",
    kernel_param=10,
    phantom=data["exact"], # Phantom was "pc".
    data=data["data"], # data was None (i.e. generated from exact phantom).
    noise_std=0.05,
    noise_type="gaussian"
)
parameters = TP.exactSolution # Extract test parameters to evaluate on

# Parse command line arguments
parser = argparse.ArgumentParser(description='Model output test.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

# Define models for each type of prior
print('Setting up models')
model_Gaussian = umbridge.HTTPModel(args.url, "Deconvolution1D_Gaussian")
model_GMRF = umbridge.HTTPModel(args.url, "Deconvolution1D_GMRF")
model_CMRF = umbridge.HTTPModel(args.url, "Deconvolution1D_CMRF")
model_LMRF = umbridge.HTTPModel(args.url, "Deconvolution1D_LMRF")

# Model for exact solution
model_exactSolution = umbridge.HTTPModel(args.url, "Deconvolution1D_ExactSolution")

# Check basic model properties
print('Checking model properties')
assert model_Gaussian.get_input_sizes() == [128]
assert model_GMRF.get_input_sizes() == [128]
assert model_CMRF.get_input_sizes() == [128]
assert model_LMRF.get_input_sizes() == [128]
assert model_exactSolution.get_input_sizes() == [0]

assert model_Gaussian.get_output_sizes() == [1]
assert model_GMRF.get_output_sizes() == [1]
assert model_CMRF.get_output_sizes() == [1]
assert model_LMRF.get_output_sizes() == [1]
assert model_exactSolution.get_output_sizes() == [128]

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
TP.prior = cuqi.distribution.Gaussian(np.zeros(128), 0.01)
assert output_Gaussian == pytest.approx(TP.posterior.logpdf(parameters))

TP.prior = cuqi.distribution.GMRF(np.zeros(128), 1/(0.01))
assert output_GMRF == pytest.approx(TP.posterior.logpdf(parameters))

TP.prior = cuqi.distribution.Cauchy_diff(np.zeros(128), 0.01)
assert output_CMRF == pytest.approx(TP.posterior.logpdf(parameters))

TP.prior = cuqi.distribution.Laplace_diff(np.zeros(128), 0.01)
assert output_LMRF == pytest.approx(TP.posterior.logpdf(parameters))

assert np.allclose(output_exactSolution, TP.exactSolution)

print('Regression testing against known norm values of output')
assert np.linalg.norm(output_Gaussian) == pytest.approx(527.3167297418795)
assert np.linalg.norm(output_GMRF) == pytest.approx(275.11317646030136)
assert np.linalg.norm(output_CMRF) == pytest.approx(623.5524285890178)
assert np.linalg.norm(output_LMRF) == pytest.approx(500.2274783053115)
assert np.linalg.norm(output_exactSolution) == pytest.approx(4.242640687119285)

print('All tests passed')

