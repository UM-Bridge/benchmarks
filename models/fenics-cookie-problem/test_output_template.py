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
modelB = umbridge.HTTPModel(args.url, "benchmark")

#test get methods
output = model.get_input_sizes()
print("get_input_sizes() returns "+str(output[0]))
assert pytest.approx(output[0]) == 8, "get_input_sizes() returns wrong value"


output = model.get_output_sizes()
print("get_output_sizes() returns "+str(output[0]))
assert pytest.approx(output[0]) == 1, "get input sizes returns wrong value"

txt = ""

for ii in range(1,16,1):
    quad_degree=ii
    print(quad_degree,end=',')

    # #test output for another config
    # txt = "model output (quantity of interest) = "
    param = [[-2.3939786473373e-01, -8.6610045659126e-01, -2.1086275315687e-01, -9.2604304103162e-01, -6.0002531612112e-01, -5.5677423053456e-01, -7.7546408441658e-01, -7.6957620518706e-01]]
    output = model(param,{"BasisDegree": 1, "Fidelity": 1, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')
    output = model(param,{"BasisDegree": 1, "Fidelity": 2, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')
    output = model(param,{"BasisDegree": 1, "Fidelity": 3, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')
    output = model(param,{"BasisDegree": 1, "Fidelity": 4, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')


    param = [[-0.2,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8]]
    output = model(param,{"BasisDegree": 1, "Fidelity": 1, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')
    output = model(param,{"BasisDegree": 1, "Fidelity": 2, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')
    output = model(param,{"BasisDegree": 1, "Fidelity": 3, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')
    output = model(param,{"BasisDegree": 1, "Fidelity": 4, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')

    param = [[0,0,0,0,0,0,0,0]]
    output = model(param,{"BasisDegree": 1, "Fidelity": 1, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')
    output = model(param,{"BasisDegree": 1, "Fidelity": 2, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')
    output = model(param,{"BasisDegree": 1, "Fidelity": 3, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end=',')
    output = model(param,{"BasisDegree": 1, "Fidelity": 4, "quad_degree":quad_degree})
    print(txt+str(output[0][0]),end='\n')
