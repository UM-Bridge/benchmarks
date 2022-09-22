#!/usr/bin/env python3
import argparse
import umbridge
import pytest

parser = argparse.ArgumentParser(description='Model output test.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

model = umbridge.HTTPModel(args.url, "forward")

# Model l = 0

output = model([[0.0,0.0]],{"level" : 0, "vtk-output" : False, "verbose" : False})[0]
print(output)

assert pytest.approx(output[0], 18) == 1847.590542295737578, "Tsunami output not as expected"
assert pytest.approx(output[1], 0.00001) == 0.001675019006365, "Tsunami output not as expected"
assert pytest.approx(output[2], 55) == 5492.84, "Tsunami output not as expected"
assert pytest.approx(output[3], 0.000003) == 0.000329495, "Tsunami output not as expected"

# Model l = 1

output = model([[0.0,0.0]],{"level" : 1, "vtk-output" : False, "verbose" : False})[0]
print(output)

assert pytest.approx(output[0], 18) == 1710.68767768625, "Tsunami output not as expected"
assert pytest.approx(output[1], 0.00001) == 0.000693793982096, "Tsunami output not as expected"
assert pytest.approx(output[2], 55) == 5075.04, "Tsunami output not as expected"
assert pytest.approx(output[3], 0.000003) == 0.000180066, "Tsunami output not as expected"
