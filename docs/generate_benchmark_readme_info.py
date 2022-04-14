import umbridge

model = umbridge.HTTPModel("http://localhost:4243")

model.supports_evaluate()
model.supports_gradient()
model.supports_apply_jacobian()
model.supports_apply_hessian()

print("Value | Dimensions")
print("---|---")
print(f"inputSizes | {model.get_input_sizes()}")
print(f"outputSizes | {model.get_output_sizes()}")

print("")

print("Feature | Supported")
print("---|---")
print(f"Evaluate | {model.supports_evaluate()}")
print(f"Gradient | {model.supports_gradient()}")
print(f"ApplyJacobian | {model.supports_apply_jacobian()}")
print(f"ApplyHessian | {model.supports_apply_hessian()}")
