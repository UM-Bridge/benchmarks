# Heat1D: 1D heat Bayesian inverse problem using CUQIpy

## Overview
This benchmark is built using the library [CUQIpy](https://cuqi-dtu.github.io/CUQIpy/). It defines a posterior distribution for a 1D heat inverse problem, with a Gaussian likelihood and Karhunen–Loève (KL) parameterization of the uncertain parameters. Two posteriors are available, one with data available everywhere in the domain and with a noise level of $0.1\%$ and the other with data available only on the left half of the domain and with a noise level of $5\%$ [... add reference CUQIpy paper 2, which cases ...].


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
output | [100] | Returns the exact solution $\mathbf{g}^\text{exact}$ for the heat 1D inverse problem (the discretized initial heat profile).

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



This benchmark defines the posterior distribution of a Bayesian inverse problem governed by a 1D heat equation. The underlying inverse problem is to infer an initial temperature profile $g(\xi)$ at time $\tau=0$ on the unit interval $\xi \in [0,1]$ form measurements of temperature $u(\xi, \tau)$ at time $\tau=\tau^\text{max}$.

In this example, the governing partial differential equation that can be solved to find the temperature $u(\xi, \tau)$ given the initial temperature profile $g(\xi)$ is
$$
\begin{align}
    \frac{\partial u(\xi,\tau)}{\partial t} -  \frac{\partial^2 u(\xi,\tau)}{\partial \xi^2}   & = 0, \quad \xi\in[0,1], \quad 0\le \tau \le \tau^\mathrm{max}, \\
    u(0,\tau)= u(1,\tau)&= 0,\\
    u(\xi,0)&= g(\xi), 
\end{align}
$$

assuming zero source term and a constant diffusion coefficient of value 1. We discretize the system using first order finite difference method in space on a regular grid of $n$ nodes, and forward Euler in time. We denote by $\mathbf{g}$ and $\mathbf{u}$ the discretization of $\mathbf{g}$ and $\mathbf{u}$.

Furthermore, we parameterize $\mathbf{g}$ using a truncated Karhunen–Loève (KL) expansion to impose some regularity and spatial correlation and reduce the dimension of the discretized unknown parameter from $n$ to $n_\text{KL}$, where $n_\text{KL} \ll n$. For given expansion basis $\mathbf{a}_i$ for $i=1,...,n_\text{KL}$, parameterization $\mathbf{g}$ in terms of the KL expansion coefficients  $\mathbf{x}=[x_1, ..., x_{n_\text{KL}}]$  can be written as  

$$
\mathbf g =\sum_{i=1}^{n_\text{KL}} x_i   \mathbf{a}_i.
$$

In the Bayesian setting, we consider  $\mathbf{x}$  an $n_\text{KL}$-dimensional random vector representing the unknown KL coefficients. And thus the corresponding $\mathbf{g}$ and $\mathbf{u}$ are random vectors as well. We define the inverse problem as 

$$
\mathbf{b} = \mathbf{A}(\mathbf{x}) + \mathbf{e},
$$

where $\mathbf{b}$ is an $m$-dimensional random vector representing the observed data, in this case measured temperature profile at time $\tau^\text{max}$, and $\mathbf{e}$ is an $m$-dimensional random vector representing the noise. $\mathbf{A}$ is the forward operator that maps $\mathbf{x}$ to $\mathbf{b}$. Applying $\mathbf{A}$ involves solving the heat PDE above for given $\mathbf{x}$ and extracting the observations $\mathbf{b}$ from the space-time solution $\mathbf{u}$.  


This benchmark defines a posterior distribution over $\mathbf{x}$ given $\mathbf{b}$ as

$$
\pi(\mathbf{x}\mid \mathbf{b}) \propto \pi(\mathbf{b}\mid \mathbf{x})\pi(\mathbf{x}),
$$

where $\pi(\mathbf{b}|\mathbf{x})$ is a likelihood function and $\pi(\mathbf{x})$ is a prior distribution.

The noise is assumed to be Gaussian with a known noise level, and so the likelihood is defined via

$$
\mathbf{b} \mid \mathbf{x} \sim \mathcal{N}(\mathbf{A}\mathbf{x}, \sigma^2\mathbf{I}_m),
$$

where $\mathbf{I}_m$ is the $m\times m$ identity matrix and $\sigma$ defines the noise level.

Two setups of this Bayesian problem are available

- `Heat1DLargeNoise`: A posterior for which the data is available everywhere in the domain and at time $\tau^\text{max}$, with a noise level of $0.1\%$

- `Heat1DSmallNoise` A posterior for which the data is available only on the left half of the domain  and at time $\tau^\text{max}$, with a noise level of $5\%$

See [um-bridge Clients](https://um-bridge-benchmarks.readthedocs.io/en/docs/umbridge/clients.html) for more details.

In addition to the two HTTP models for the posterior, there is also an HTTP model for the exact solution to the problem. This model is called `Heat1DExactSolution` and returns exact initial heat profile used to generate the synthetic data when called. The map from the coefficients $\mathbf{x}$ to the discretized function $\mathbf{g}$ is provided via the HTTP model `KLExpansionCoefficient2Function` and the projection of $\mathbf{g}$ on the coefficient space $\mathbf{x}$ is provided by the HTTP model `KLExpansionFunction2Coefficient`. 

Using [CUQIpy](https://cuqi-dtu.github.io/CUQIpy/) this benchmark is defined in the files `heat1D_problem.py`, `data_script.py`, and `server.py` provided here.
