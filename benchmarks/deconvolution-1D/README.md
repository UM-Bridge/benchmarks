# Deconvolution1D

## Overview
This benchmark is based on a [1D Deconvolution test problem](https://cuqi-dtu.github.io/CUQIpy/api/_autosummary/cuqi.testproblem/cuqi.testproblem.Deconvolution1D.htm) from the library [CUQIpy](https://cuqi-dtu.github.io/CUQIpy/).

It defines a posterior distribution that characterizes the probability distribution of a clean signal determined from a noisy blurred observation of the signal.

The noise is assumed the be Gaussian, and so the likelihood function comes from a Gaussian distribution. The prior can be configured by choosing of the the multiple choices of HTTP models.

- `Gaussian`: Gaussian (Normal) distribution
- `GMRF` Gaussian Markov Random Field
- `CMRF` Cauchy Markov Random Field
- `LMRF` Laplace Markov Random Field

## Authors
- [Nicolai A. B. Riis](mailto:nabr@dtu.dk)
- [Jakob S. Jørgensen](mailto:jakj@dtu.dk)
- [Amal M. Alghamdi](mailto:amaal@dtu.dk)

## Run
```
docker run -it -p 4243:4243 linusseelinger/benchmark-deconvolution-1D
```

## Call

Using `umbridge` for python:
```python
import umbridge

# Set url and name of the HTTP model
url = "localhost" # Set this to url of docker image
name = "Deconvolution1D_Gaussian" # Switch prior here

# Define model
model = umbridge.HTTPModel(url, name)

# Get model info from server
print(model.info())
```

## Properties

Mapping | Dimensions | Description
---|---|---
input | [128] | Signal $x$
output | [1] | Log PDF $\pi$ evaluated at $x$

Feature | Supported
---|---
Evaluate | True
Gradient | True (in most cases)
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
None

Mount directory | Purpose
---|---
None |

## Description

TODO or maybe ref to CUQIpy?
