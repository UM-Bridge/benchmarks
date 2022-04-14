#!/usr/bin/env python3
import argparse
import umbridge
import pytest

parser = argparse.ArgumentParser(description='Minimal HTTP model demo.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

# Set up a model by connecting to URL
model = umbridge.HTTPModel(args.url)

print(model.get_input_sizes())
print(model.get_output_sizes())

# Model evaluation with configuration parameters
output = model([[0.0,0.0]],{"level" : 0, "vtk-output" : False, "verbose" : True})[0]

print(output)

assert pytest.approx(output[0], 18) == 1847.59, "Tsunami output not as expected"
assert pytest.approx(output[1], 0.00001) == 0.0016, "Tsunami output not as expected"
assert pytest.approx(output[2], 55) == 5492.84, "Tsunami output not as expected"
assert pytest.approx(output[3], 0.000003) == 0.00032, "Tsunami output not as expected"

