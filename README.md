![UM-bridge_map](https://raw.githubusercontent.com/UM-Bridge/benchmarks/main/UM-bridge_map.png "UQ-Model-UM")

**UM-Bridge** (the **U**Q and **M**odel **Bridge**) provides a unified interface for numerical models that is accessible from virtually any programming language or framework. It is primarily intended for coupling advanced models (e.g. simulations of complex physical processes) to advanced statistics or optimization methods for UQ.

Many uncertainty quantification (UQ) and optimization methods treat a model as an abstract function and only interact with the model through operations such as *simple model evaluation*, *gradient evaluation* or *Jacobian action*.

![UQ-Model-UM](https://raw.githubusercontent.com/UM-Bridge/benchmarks/main/UQ-Model-UM.png "UQ-Model-UM")

The key idea of UM-Bridge is to provide the mathematical "interface" as an abstract interface in software as well. By using **HTTP** behind the scenes, a high degree of flexibility is achieved, allowing for:

* **Coupling** of codes written in arbitrary languages and frameworks
* **Accelerating** development of advanced software stacks
* **Containarization** of models &rarr; easier collaboration and &rarr; access to container-based compute resourcs in the cloud (e.g., GCP, AWS)
* **Portability** across operating systems
* **Fully reproducible** models and benchmarks.

To get started the following minimal Python code will give you a first impression:

```
import umbridge
url = "http://testmodel.linusseelinger.de"
model = umbridge.HTTPModel(url, "forward")
print(model([[100]]))
```

This passes an input to a simple 1D test model running on a remote server and prints the model's output.

See [quickstart guide](https://um-bridge-benchmarks.readthedocs.io/en/docs/quickstart.html) and [tutorial](https://um-bridge-benchmarks.readthedocs.io/en/docs/tutorial.html) for more information!

## Languages and frameworks

These tables show what languages and frameworks UM-Bridge currently provides integrations for.

"Server" refers to the model side, while "client" is the UQ / statistics / optimization side. We are happy to actively support the development of new integrations.

Language | Client (UQ) | Server (model)
---|---|---
C++ | ✓ | ✓
MATLAB | ✓ | ✗
Python | ✓ | ✓
R | ✓ | ✗
Julia | ✓ | ✗

Framework | Client (UQ) | Server (model)
---|---|---
CUQIpy | ✓ | ✓
emcee | ✓ | ✗
MUQ | ✓ | ✓
PyApprox | ✓ | ✗
PyMC | ✓ | ✗
QMCPy | ✓ | ✗
Sparse Grids MATLAB Kit | ✓ | ✗
tinyDA | ✓ | ✗
TT Toolbox | ✓ | ✗
UQPy | ✓ | ✗

## Opinions

<figure style="display: flex; align-items: center;">
  <img src="https://raw.githubusercontent.com/UM-Bridge/benchmarks/main/mikkel_lykkegaard_picture.png" alt="Person Image" style="width:100px;height:100px;margin-right:10px;border-radius:50%;">
  <figcaption>
  
> *I love Uncertainty Quantification. But I don't love fiddling around with complex numerical solver routines. UM-Bridge takes the pain away from doing UQ with complex models.*

> Dr Mikkel Bue Lykkegaard - Data Science Lead, digiLab

  </figcaption>
</figure>

## Citation

Please cite: Seelinger et al., (2023). UM-Bridge: Uncertainty quantification and modeling bridge. Journal of Open Source Software, 8(83), 4748, [url](https://doi.org/10.21105/joss.04748).

```
@article{UMBridge, doi = {10.21105/joss.04748}, url = {https://doi.org/10.21105/joss.04748}, year = {2023}, publisher = {The Open Journal}, volume = {8}, number = {83}, pages = {4748}, author = {Linus Seelinger and Vivian Cheng-Seelinger and Andrew Davis and Matthew Parno and Anne Reinarz}, title = {UM-Bridge: Uncertainty quantification and modeling bridge}, journal = {Journal of Open Source Software} }
```

## Resources

This repository hosts stand-alone reference problems for benchmarking of UQ algorithms.

The documentation for the project and for the benchmark problems within is at: [Documentation](https://um-bridge-benchmarks.readthedocs.io/en/docs/).

The project source code is [hosted on GitHub](https://github.com/UM-Bridge).

Join the project [slack channel](https://join.slack.com/t/um-bridge/shared_invite/zt-1da1ebkly-8s0YQdZUIYkJ1vws6edsAQ).
