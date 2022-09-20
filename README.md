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

This table shows what languages and frameworks UM-Bridge currently provides integrations for. Note that "server" refers to the model side, while "client" is the uncertainty quantification / statistics / optimization side. We are happy to actively support the development of new integrations.

Language / framework | Client support | Server support
---|---|---
C++ | 🗸 | 🗸
Python | 🗸 | 🗸
R | 🗸 | ✗
MUQ | 🗸 | 🗸
PyMC (4.x) | 🗸 | ✗


This repository hosts stand-alone reference problems for benchmarking of UQ algorithms.

The documentation for the project and for the benchmark problems within is at: [Documentation](https://um-bridge-benchmarks.readthedocs.io/en/docs/).

The project source code is [hosted on GitHub](https://github.com/UM-Bridge).

Join the project [slack channel](https://join.slack.com/t/um-bridge/shared_invite/zt-1da1ebkly-8s0YQdZUIYkJ1vws6edsAQ).
