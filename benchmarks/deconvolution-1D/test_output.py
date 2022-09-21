#!/usr/bin/env python3
import argparse
import umbridge
import pytest
import cuqi
import numpy as np

parser = argparse.ArgumentParser(description='Model output test.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

model = umbridge.HTTPModel(args.url, "posterior")

# Load test problem (same seed as in server.py)
np.random.seed(0)
TP = cuqi.testproblem.Deconvolution1D(dim=128)

# Extract test parameters to evaluate on
parameters = TP.exactSolution

# Test Info
print(f"Model input and output: {model.get_input_sizes()}, {model.get_output_sizes()}")

# Supports evaluate?
print(f"Model supports evaluate: {model.supports_evaluate()}")

# Supports gradient?
print(f"Model supports gradient: {model.supports_gradient()}")

# Test evaluate
print("Evaluating model...")
print(model([parameters.tolist()]) == pytest.approx(TP.posterior.logpdf(parameters)))
print("Finished evaluating model.")

