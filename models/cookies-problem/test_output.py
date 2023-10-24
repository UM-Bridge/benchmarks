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
#model = umbridge.HTTPModel(args.url)
model = umbridge.HTTPModel(args.url, "forward")
#model = umbridge.HTTPModel(args.url, "posterior")


print(model.get_input_sizes())
print(model.get_output_sizes())

#param = [[-0.2,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8]]
param = [[-2.3939786473373e-01 , -8.6610045659126e-01, -2.1086275315687e-01 , -9.2604304103162e-01 , -6.0002531612112e-01 , -5.5677423053456e-01 , 
 -7.7546408441658e-01 , -7.6957620518706e-01]]
print(param)


# Simple model evaluation
#print(model(param))
#print(model(param,{"NumThreads": 10, "BasisDegree": 4, "Fidelity": 2}))
output = model(param,{"NumThreads": 10, "BasisDegree": 4, "Fidelity": 2})

print(output)
assert pytest.approx(output[0][0]) == 0.06932827462480169, "Output not as expected"

