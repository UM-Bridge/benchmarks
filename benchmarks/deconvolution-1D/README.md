# Deconvolution1D

## Overview
This benchmark is based on a Deconvolution1D problem.

Gaussian prior, Gaussian likelihood, and Gaussian posterior.

## Authors
- [Nicolai A. B. Riis](mailto:nabr@dtu.dk)

## Run
```
docker run -it -p 4243:4243 linusseelinger/benchmark-analytic-banana
```

## Properties

Mapping | Dimensions | Description
---|---|---
input | [128] | Signal $x$
output | [1] | Log PDF $\pi$ evaluated at $x$

Feature | Supported
---|---
Evaluate | True
Gradient | True
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
a | double | 2.0 | Transformation parameter
b | double | 0.2 | Transformation parameter
scale | double | 1.0 | Scaling factor applied to the underlying normal distribution's variance

Mount directory | Purpose
---|---
None |

## Description

..
