import umbridge

model = umbridge.HTTPModel("http://localhost:4242")

print("# Model / Benchmark Name")

print("## Overview")

print("## Authors")

print("## Run")

print("## Properties")
print("")

print("Mapping | Dimensions | Description")
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

print("Mount directory | Purpose")
print("---|---")
print("None |")

print("## Description")