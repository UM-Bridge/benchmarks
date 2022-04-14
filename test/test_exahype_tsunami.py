#!/usr/bin/env python3
import argparse
import requests

parser = argparse.ArgumentParser(description='ExaHyPE Tsunami Test.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

print("Requesting input sizes...")
r = requests.get(f"{args.url}/GetInputSizes")
rInputSizes = r.json()
print(rInputSizes)

print("Requesting output sizes...")
r = requests.get(f"{args.url}/GetOutputSizes")
print(r.text)

print("Requesting info...")
r = requests.get(f"{args.url}/Info")
print(r.text)

print("Requesting evaluation")

# Build input parameter vectors of dimensions expected by model, fill with zeros for testing
inputParams = {"input": [], "config": {"level" : 0, "vtk-output" : False, "verbose" : True}}
for i in range(0,len(rInputSizes["inputSizes"])):
  inputParams["input"].append([0] * rInputSizes["inputSizes"][i])
print(inputParams)

r = requests.post(f"{args.url}/Evaluate", json=inputParams)
print(r.text)

output = r.json()["output"][0]

if abs(output[0] - 1847.59) > 18:
    raise SystemExit("Tsunami output not as expected")

if abs(output[1] - 0.0016) > 0.00001:
    raise SystemExit("Tsunami output not as expected")

if abs(output[2] - 5492.84) > 55:
    raise SystemExit("Tsunami output not as expected")

if abs(output[3] - 0.00032) > 0.000003:
    raise SystemExit("Tsunami output not as expected")

