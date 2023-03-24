UM-Bridge (the UQ and Model Bridge) provides a unified interface for numerical models that is accessible from virtually any programming language or framework. It is primarily intended for coupling advanced models (e.g. simulations of complex physical processes) to advanced statistics or optimization methods.

In many statistics / uncertainty quantification or optimization methods, the model only appears as a function mapping vectors onto vectors with some of the following:
* Simple evaluation,
* Gradient evaluation,
* Jacobian action,
* Hessian action.

The key idea of UM-Bridge is to now provide this mathematical "interface" as an abstract interface in software as well. By using HTTP behind the scenes, a high degree of flexibility is achieved, allowing for:

* Coupling of codes written in arbitrary languages and frameworks, accelerating development of advanced software stacks combining the state-of-the art of modelling with statistics / optimization.
* Containarization of models, making collaboration easier due to portability of models and separation of concerns between fields (specifically model and statistics experts).
* Unified, portable, fully reproducible and black-box benchmark problems defined software.

## Languages and frameworks

This table shows what languages and frameworks UM-Bridge currently provides integrations for. Note that "server" refers to the model side, while "client" is the uncertainty quantification / statistics / optimization side. We are happy to actively support the development of new integrations.

Language / framework | Client support | Server support
---|---|---
C++ | ✓ | ✓
MATLAB | planned | ✗
Python | ✓ | ✓
R | ✓ | ✗
MUQ | ✓ | ✓
PyMC (4.x) | ✓ | ✗
QMCPy | ✓ | ✗
Sparse Grids MATLAB Kit | planned | ✗
tinyDA | ✓ | ✗

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
