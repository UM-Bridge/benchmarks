#!/usr/bin/env python3
""" Script to regression test the output of the Heat1D benchmark. """

import argparse
import umbridge
import pytest
import cuqi
import numpy as np
from heat1D_problem import create_forward_PDE_model, obs_everywhere_indices,\
    obs_left_half_indices, dim, N, sigma, large_noise_level, small_noise_level, prior_mean, prior_cov

# A. Set up CUQIpy testproblem to compare with
print('Setting up heat 1D problem')

# Small noise case:
data_dic_small_noise = np.load("data_small_noise.npz")
model_small_noise = create_forward_PDE_model(obs_everywhere_indices)
domain_geometry_small_noise = model_small_noise.domain_geometry
range_geometry_small_noise = model_small_noise.range_geometry
x_small_noise = cuqi.distribution.Gaussian(
    mean=prior_mean, cov=prior_cov, geometry=domain_geometry_small_noise)
y_small_noise = cuqi.distribution.Gaussian(model_small_noise(x_small_noise), sigma(
    small_noise_level, data_dic_small_noise['y_exact'])**2*np.eye(model_small_noise.range_dim),
    geometry=range_geometry_small_noise)
posterior_small_noise = cuqi.distribution.JointDistribution(
    x_small_noise, y_small_noise)(y_small_noise=data_dic_small_noise['data'])

# Large noise case:
data_dic_large_noise = np.load("data_large_noise.npz")
model_large_noise = create_forward_PDE_model(obs_left_half_indices)
domain_geometry_large_noise = model_large_noise.domain_geometry
range_geometry_large_noise = model_large_noise.range_geometry
x_large_noise = cuqi.distribution.Gaussian(
    mean=prior_mean, cov=prior_cov, geometry=domain_geometry_large_noise)
y_large_noise = cuqi.distribution.Gaussian(model_large_noise(x_large_noise), sigma(
    large_noise_level, data_dic_large_noise['y_exact'])**2*np.eye(model_large_noise.range_dim),
    geometry=range_geometry_large_noise)
posterior_large_noise = cuqi.distribution.JointDistribution(
    x_large_noise, y_large_noise)(y_large_noise=data_dic_large_noise['data'])


# B. Parse command line arguments
parser = argparse.ArgumentParser(description='Model output test.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

# C. Define umbridge models for benchmark cases
print('Setting up models')
um_model_large_noise = umbridge.HTTPModel(args.url, "Heat1DLargeNoise")
um_model_small_noise = umbridge.HTTPModel(args.url, "Heat1DSmallNoise")


# D. Define umbridge for exact solution
um_model_exactSolution = umbridge.HTTPModel(args.url, "Heat1DExactSolution")

# E. Define umbridge model for KL expansion
um_kl_expansion_coefficient2function = umbridge.HTTPModel(
    args.url, "KLExpansionCoefficient2Function")
um_kl_expansion_function2coefficient = umbridge.HTTPModel(
    args.url, "KLExpansionFunction2Coefficient")

## F. Check basic model properties
print('Checking model properties')
assert um_model_large_noise.get_input_sizes() == [20]
assert um_model_small_noise.get_input_sizes() == [20]
assert um_model_exactSolution.get_input_sizes() == [0]
assert um_kl_expansion_coefficient2function.get_input_sizes() == [20]
assert um_kl_expansion_function2coefficient.get_input_sizes() == [100]

assert um_model_large_noise.get_output_sizes() == [1]
assert um_model_small_noise.get_output_sizes() == [1]
assert um_model_exactSolution.get_output_sizes() == [100]
assert um_kl_expansion_coefficient2function.get_output_sizes() == [100] 
assert um_kl_expansion_function2coefficient.get_output_sizes() == [20]



## G. Check that the models support evaluation and gradient
print('Checking model capabilities')
assert um_model_large_noise.supports_evaluate()
assert um_model_small_noise.supports_evaluate()
assert um_model_exactSolution.supports_evaluate()
assert um_kl_expansion_coefficient2function.supports_evaluate()
assert um_kl_expansion_function2coefficient.supports_evaluate()

assert not um_model_large_noise.supports_gradient()
assert not um_model_small_noise.supports_gradient()
assert not um_model_exactSolution.supports_gradient()
assert not um_kl_expansion_coefficient2function.supports_gradient()
assert not um_kl_expansion_function2coefficient.supports_gradient()


## H. Evaluate models and compare with probelm set up in A
print('Evaluating models')
np.random.seed(0)
parameters = np.random.rand(20)
output_large_noise = um_model_large_noise([parameters.tolist()])[0][0]
output_small_noise = um_model_small_noise([parameters.tolist()])[0][0]
output_exactSolution = um_model_exactSolution([[]])[0]
output_kl_expansion_coefficient2function = um_kl_expansion_coefficient2function([parameters.tolist()])[0]
output_kl_expansion_function2coefficient = um_kl_expansion_function2coefficient([output_kl_expansion_coefficient2function])[0]

## Check that the output is correct
print('Checking model output')
assert np.allclose(output_large_noise, posterior_large_noise.logpdf(parameters))
assert np.allclose(output_small_noise, posterior_small_noise.logpdf(parameters))
assert np.allclose(output_exactSolution, data_dic_large_noise['x_exact'])
assert np.allclose(output_kl_expansion_coefficient2function, domain_geometry_large_noise.par2fun(parameters))
assert np.allclose(output_kl_expansion_function2coefficient, domain_geometry_large_noise.fun2par(output_kl_expansion_coefficient2function))

print('All tests passed')
