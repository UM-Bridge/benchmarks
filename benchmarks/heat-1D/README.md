# Heat1D

## Overview
This benchmark is built using the library [CUQIpy](https://cuqi-dtu.github.io/CUQIpy/). It defines a posterior distribution for a 1D heat inverse problem, with a Gaussian likelihood and Karhunen–Loève (KL) parameterization of the uncertain parameters. Two posteriors are available, one with data everywhere in the domain and with a noise level of $0.1\%$ and the other with data available only on the left half of the domain and with a noise level of $5\%$ [... add reference CUQIpy paper 2, which cases ...].


Plot of data and exact solution for the two noise-level cases

![Data](data.png "Data")

Credibility interval plots of posterior samples for the two noise-level cases

![Samples](samples.png "Credibility interval of samples")

## Authors
- [Nicolai A. B. Riis](mailto:nabr@dtu.dk)
- [Jakob S. Jørgensen](mailto:jakj@dtu.dk)
- [Amal M. Alghamdi](mailto:amaal@dtu.dk)

## Run
```
docker run -it -p 4243:4243 linusseelinger/benchmark-heat-1D
```

## Properties

Model | Description
---|---
Heat1DSmallNoise | Posterior distribution for the small noise case
Heat1DLargeNoise | Posterior distribution for the large noise and incomplete data case
Heat1DExactSolution | Exact solution to the 1D heat problem
KLExpansionCoefficient2Function | Map from the KL parameter space to the function space
KLExpansionFunction2Coefficient | Projection from the function space to the KL parameter space


### Heat1DSmallNoise
Mapping | Dimensions | Description
---|---|---
input | [20] | KL Coefficients $\mathbf{x}$
output | [1] | Log PDF $\pi(\mathbf{x}\mid\mathbf{b})$ of posterior for the small noise case

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
None | | |

### Heat1DLargeNoise
Mapping | Dimensions | Description
---|---|---
input | [20] | KL coefficients $\mathbf{x}$
output | [1] | Log PDF $\pi(\mathbf{x}\mid\mathbf{b})$ of posterior for the large noise and incomplete data case

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
None | | |


### Heat1DExactSolution
Mapping | Dimensions | Description
---|---|---
input | [0] | No input to be provided. 
output | [100] | Returns the exact solution $\mathbf{x}$ for the heat 1D problem (defined on the function space).

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
None | | |


### KLExpansionCoefficient2Function
Mapping | Dimensions | Description
---|---|---
input | [20] | KL coefficients $\mathbf{x}$
output | [100] | The function values on grid points: the result of the mapping from the KL coefficients space to the function space 

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
None | | |


### KLExpansionFunction2Coefficient
Mapping | Dimensions | Description
---|---|---
input | [100] | Function values on grid points
output | [20] | The KL coefficients: the result of projecting the function values onto the KL coefficients space 

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
None | | |

## Mount directories
Mount directory | Purpose
---|---
None |

## Source code

[Sources here.](./heat1D_problem.py)

## Description

