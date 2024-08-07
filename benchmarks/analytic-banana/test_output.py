#!/usr/bin/env python3
import argparse
import umbridge
import pytest

parser = argparse.ArgumentParser(description='Model output test.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

model = umbridge.HTTPModel(args.url, "posterior")


output = model([[0.0,0.0]])[0]
print(output)
assert pytest.approx(output[0]) == -5.534036030183454, "Output not as expected"

output = model([[2.0,0.0]])[0]
print(output)
assert pytest.approx(output[0]) == -3.3207026968501214, "Output not as expected"
