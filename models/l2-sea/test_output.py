#!/usr/bin/env python3

# run this script as python3 umbridge-client.py http://localhost:4242

import argparse
import umbridge
import pytest

parser = argparse.ArgumentParser(description='Minimal HTTP model demo.')
parser.add_argument('url',metavar='url',type=str,
                    help='the URL on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

# Set up a model by connecting to URL
model = umbridge.HTTPModel(args.url,"forward")

#test get method
output = model.get_input_sizes()
print(output)
assert pytest.approx(output[0]) == 16, "get_input_sizes() returns wrong value"

#test model output
param = [[0.28, -6.16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
output = model(param, {"fidelity": 7, "sinkoff":'y', "trimoff":'y'})

print(output)
assert pytest.approx(output[0][0]) == 48.9337769, "Output not as expected"
assert pytest.approx(output[0][1]) == -1.00000000, "Output not as expected"
assert pytest.approx(output[0][2]) == -1.00000000, "Output not as expected"
assert pytest.approx(output[0][3]) == -0.107635260, "Output not as expected"
assert pytest.approx(output[0][4]) == -1.72240210, "Output not as expected"


#another test, this time for the benchmark UQ version
model = umbridge.HTTPModel(args.url,"benchmark_UQ")

param = [[0.28,-6.16]]
output = model(param, {"fidelity": 7, "sinkoff":'y', "trimoff":'y'})
print(output)
assert pytest.approx(output[0][0]) == 48.9337769, "Output not as expected"
assert pytest.approx(output[0][1]) == -1.00000000, "Output not as expected"
assert pytest.approx(output[0][2]) == -1.00000000, "Output not as expected"
assert pytest.approx(output[0][3]) == -0.107635260, "Output not as expected"
assert pytest.approx(output[0][4]) == -1.72240210, "Output not as expected"

#another test, this time for the benchmark OPT version
model = umbridge.HTTPModel(args.url,"benchmark_OPT")

param = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
output = model(param, {"fidelity": 7, "sinkoff":'y', "trimoff":'y'})
print(output)
assert pytest.approx(output[0][0]) == 48.9337769, "Output not as expected"
assert pytest.approx(output[0][1]) == -1.00000000, "Output not as expected"
assert pytest.approx(output[0][2]) == -1.00000000, "Output not as expected"
assert pytest.approx(output[0][3]) == -0.107635260, "Output not as expected"
assert pytest.approx(output[0][4]) == -1.72240210, "Output not as expected"

