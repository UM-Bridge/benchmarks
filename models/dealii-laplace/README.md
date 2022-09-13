# A Bayesian inverse problem benchmark based on the Laplace equation
## Overview

This model implements the benchmark described in [this
preprint](https://arxiv.org/abs/2102.07263). In it, the 64 inputs
correspond to an 8x8 grid of stiffness values of a membrane, which
deforms under a (known) form, and the 169 outputs correspond to the
values of the resulting deformation at a 13x13 grid of points.

The complete benchmark compares these 169 outputs with actual
measurements (obtained through independent computations that solve the
forward problem with a different method) and augments it with a prior
probability distribution on the set of 64 parameters. The result is a
posterior probability distribution whose properties (such as mean and
covariance) the benchmark probes.

The complete benchmark, along with its solution, is described in great
detail in the paper mentioned above.

## Authors

David Aristoff and Wolfgang Bangerth (Colorado State University)

## Run
```
docker run -it -p 4242:4242 linusseelinger/model-laplace:latest
```

## Properties

Mapping | Dimensions | Description
---|---|---
input | [64] | A set of 64 parameters corresponding to an 8x8 grid of stiffness values of a membrane
output | [169] | A set of 169 displacement values corresponding to the displacement of the membrane at a grid of 13x13 points

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
None | | |

Mount directory | Purpose
---|---
None |

## Description

See [this preprint](https://arxiv.org/abs/2102.07263).