#!/usr/bin/env python3

# first run the container as
#
# docker run -it -p 4242:4242 <image name>
#
# then run this script as python3 test_output.py http://localhost:4242

import argparse
import umbridge
import pytest

parser = argparse.ArgumentParser(description='Minimal HTTP model demo.')
parser.add_argument('url', metavar='url', type=str, help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

# Set up a model by connecting to URL
model = umbridge.HTTPModel(args.url, "forwardparabolic")

#test get methods
output = model.get_input_sizes()
print("get_input_sizes() returns "+str(output[0]))
assert pytest.approx(output[0]) == 8, "get_input_sizes() returns wrong value"


output = model.get_output_sizes()
print("get_output_sizes() returns "+str(output[0]))
assert pytest.approx(output[0]) == 1, "get input sizes returns wrong value"


#test output for default config (thread=1, p=4, fid=2)
param = [[-2.3939786473373e-01, -8.6610045659126e-01, -2.1086275315687e-01, -9.2604304103162e-01, -6.0002531612112e-01, -5.5677423053456e-01, -7.7546408441658e-01, -7.6957620518706e-01]]
output = model(param)
print("model output (quantity of interest) for default config values = "+str(output[0][0]))
assert  pytest.approx(output[0][0], abs=1e-6) == 0.06933585905253055, "Output not as expected"

#test output for another config
output = model(param,{"BasisDegree": 3, "Fidelity": 2})
print("model output (quantity of interest) = "+str(output[0][0]))
assert  pytest.approx(output[0][0], abs=1e-6) == 0.06936772917504516, "Output not as expected"
