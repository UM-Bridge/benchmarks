# Deconvolution1D

## Overview
This benchmark is based on a [1D Deconvolution test problem](https://cuqi-dtu.github.io/CUQIpy/api/_autosummary/cuqi.testproblem/cuqi.testproblem.Deconvolution1D.html) from the library [CUQIpy](https://cuqi-dtu.github.io/CUQIpy/). It defines a posterior distribution for a 1D deconvolution problem, with a Gaussian likelihood and four different choices of prior distributions with configurable parameters.

## Authors
- [Nicolai A. B. Riis](mailto:nabr@dtu.dk)
- [Jakob S. Jørgensen](mailto:jakj@dtu.dk)
- [Amal M. Alghamdi](mailto:amaal@dtu.dk)

## Run
```
docker run -it -p 4243:4243 linusseelinger/benchmark-deconvolution-1D
```
## Properties

Mapping | Dimensions | Description
---|---|---
input | [128] | Signal $\mathbf{x}$
output | [1] | Log PDF $\pi(\mathbf{b}\mid\mathbf{x})$

Feature | Supported
---|---
Evaluate | True
Gradient | True (in most cases)
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
delta | double | 0.05 | The prior parameter $\delta$ (see below).

Mount directory | Purpose
---|---
None |

## Description


The 1D periodic deconvolution problem is defined by the inverse problem

$$
\mathbf{b} = \mathbf{A}\mathbf{x} + \mathbf{e},
$$

where $\mathbf{b}$ is an $m$-dimensional random vector representing the observed data, $\mathbf{A}$ is an $m\times n$ matrix representing the convolution operator, $\mathbf{x}$ is an $n$-dimensional random vector representing the unknown signal, and $\mathbf{e}$ is an $m$-dimensional random vector representing the noise.

This benchmark defines a posterior distribution over $\mathbf{x}$ given $\mathbf{b}$ as

$$
\pi(\mathbf{x}\mid \mathbf{b}) \propto \pi(\mathbf{b}\mid \mathbf{x})\pi(\mathbf{x}),
$$

where $\pi(\mathbf{b}|\mathbf{x})$ is a likelihood function and $\pi(\mathbf{x})$ is a prior distribution. 

The noise is assumed to be Gaussian with a known noise level, and so the likelihood is defined via

$$
\mathbf{b} \mid \mathbf{x} \sim \mathcal{N}(\mathbf{A}\mathbf{x}, \sigma^2\mathbf{I}_m),
$$

where $\mathbf{I}_m$ is the $m\times m$ identity matrix and $\sigma=0.05$ is the noise standard deviation.

The prior can be configured by choosing of the following assumptions about $\mathbf{x}$:

- `Gaussian`: Gaussian (Normal) distribution: $x_i \sim \mathcal{N}(0, \delta)$.

- `GMRF` Gaussian Markov Random Field: $x_i-x_{i-1} \sim \mathcal{N}(0, \delta)$.

- `CMRF` Cauchy Markov Random Field: $x_i-x_{i-1} \sim \mathcal{C}(0, \delta)$.

- `LMRF` Laplace Markov Random Field: $x_i-x_{i-1} \sim \mathcal{L}(0, \delta)$.

where $\mathcal{C}$ is the Cauchy distribution and $\mathcal{L}$ is the Laplace distribution. The parameter $\delta$ is the prior parameter and is configurable (see above).

The choice of prior is specified by providing the name to the HTTP model. In this case `Deconvolution1D_Gaussian`, `Deconvolution1D_GMRF`, `Deconvolution1D_CMRF`, and `Deconvolution1D_LMRF`, respectively. See [um-bridge Clients](https://um-bridge-benchmarks.readthedocs.io/en/docs/umbridge/clients.html) for more details.

In addition to the HTTP models for the posterior, there is also an HTTP model for the exact solution to the problem. This model is called `Deconvolution1D_ExactSolution` and returns exact phantom used to generate the synthetic data when called.
