#!/usr/bin/env python3

# run this script as python3 umbridge-client.py http://localhost:4242

import argparse
import umbridge
import pytest

parser = argparse.ArgumentParser(description='Minimal HTTP model demo.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

# Set up a model by connecting to URL
model = umbridge.HTTPModel(args.url, "forward")

#test get methods
output = model.get_input_sizes()
print(output)
assert pytest.approx(output[0][0]) == 8, "get_input_sizes() returns wrong value"


output = model.get_output_sizes()
print(output)
assert pytest.approx(output[0][0]) == 1, "get input sizes returns wrong value"


#test output
param = [[-2.3939786473373e-01, -8.6610045659126e-01, -2.1086275315687e-01, -9.2604304103162e-01, -6.0002531612112e-01, -5.5677423053456e-01, -7.7546408441658e-01, -7.6957620518706e-01]]
output = model(param,{"NumThreads": 10, "BasisDegree": 3, "Fidelity": 3})

print(output)
assert pytest.approx(output[0][0]) == 0.06932827462480169, "Output not as expected"


#another test, this time for the benchmark version
model = umbridge.HTTPModel(args.url, "benchmark")

param = [[-0.2,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8]]
output = model(param,{"NumThreads": 10})
print(output)
assert pytest.approx(output[0][0]) == 0.06932827462480169, "Output not as expected"

