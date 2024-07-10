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
model = umbridge.HTTPModel(args.url, "forward")

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
assert pytest.approx(output[0][0]) == 0.0693282746248043, "Output not as expected"

#test output for another config
output = model(param,{"BasisDegree": 3, "Fidelity": 3})
print("model output (quantity of interest) = "+str(output[0][0]))
assert pytest.approx(output[0][0]) == 0.06934748547844366, "Output not as expected"


#another test, this time for the benchmark version (i.e. p=4, fid=2 again, same as model default)
model_B = umbridge.HTTPModel(args.url, "benchmark")

param = [[-0.2,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8]]
output = model_B(param)
print("model output (quantity of interest) in benchmark configuration = "+str(output[0][0]))
assert pytest.approx(output[0][0]) == 0.05725269745090122, "Output not as expected"

