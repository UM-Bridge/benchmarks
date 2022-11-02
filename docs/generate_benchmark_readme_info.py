import argparse
import umbridge

# Read URL from command line argument
parser = argparse.ArgumentParser(description='Minimal HTTP model demo.')
parser.add_argument('url', metavar='url', type=str,
                    help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to URL {args.url}")

print("# TODO: Model / Benchmark Name")

print("## Overview")

print("## Authors")
print("- [Jane Doe](mailto:jane.doe@uni.edu)")
print("- [John Doe](mailto:john.doe@uni.edu)")

print("## Run")
print("```")
print("docker run -it -p 4242:4242 TODO/IMAGE_NAME")
print("```")
print("")


print("## Properties")
print("")

print("Model | Description")
print("---|---")
for model_name in umbridge.supported_models(args.url):
  print(f"{model_name} | TODO: Description")
print("")

for model_name in umbridge.supported_models(args.url):
  model = umbridge.HTTPModel(args.url, model_name)

  print(f"### {model_name}")

  print(f"Mapping | Dimensions | Description")
  print("---|---|---")
  print(f"input | {model.get_input_sizes()} | TODO: INPUT DESCRIPTION")
  print(f"output | {model.get_output_sizes()} | TODO: OUTPUT DESCRIPTION")

  print("")

  print("Feature | Supported")
  print("---|---")
  print(f"Evaluate | {model.supports_evaluate()}")
  print(f"Gradient | {model.supports_gradient()}")
  print(f"ApplyJacobian | {model.supports_apply_jacobian()}")
  print(f"ApplyHessian | {model.supports_apply_hessian()}")

  print("")

  print("Config | Type | Default | Description")
  print("---|---|---|---")
  print("None | | |")

  print("")

print("## Mount directories")

print("Mount directory | Purpose")
print("---|---")
print("None |")

print("")

print("## Description")
