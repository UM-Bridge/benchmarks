# Cookie Benchmark
## Overview
This model implements the so-called 'cookies problem' or 'cookies in the oven problem' (see for reference [[Bäck et al.,2011]](https://doi.org/10.1007/978-3-642-15337-2_3), [[Ballani et al.,2015]](https://doi.org/10.1137/140960980), [[Kressner et al., 2011]](https://doi.org/10.1137/100799010)), i.e., a simplified thermal equation in which the conductivity coefficient is uncertain in 8 circular subdomains ('the cookies'), whereas it is known (and constant) in the remaining of the domain ('the oven').

The approximation is constructed using the legacy version of [FEniCs](https://fenicsproject.org/) via the Python interface.
An approximation is constructed using a mesh of quadrilaterals.
A continuous piecewise polynomial is constructed using a Lagrange finite element basis.
The basis degree is user configurable.
All the required code is included in this repository, and it is hoped that the FEniCs code can be easily adapted by a user if required.

An UM-BRIDGE server can be run using Docker with the command
``` docker run -p 4242:4242 -itd linusseelinger/cookies-problem:latest```
where the first port number is mapped through to ```4242``` which is exposed in the container.

## Authors
- [Benjamin Kent](mailto:kent@imati.cnr.it)
- [Massimiliano Martinelli](mailto:martinelli@imati.cnr.it)
- [Lorenzo Tamellini](mailto:tamellini@imati.cnr.it)

## Run 
```
docker run -p 4242:4242 -it linusseelinger/cookies-problem:latest
```
The compressed size of the container is 754.5 MB.

The Docker container serves four models.
The elliptic models are:
- forward: the forward elliptic model with configurable model inputs
- benchmark: the forward elliptic model with fixed model inputs

A parabolic formulation (i.e. $\frac{\partial u}{\partial t} + \nabla \cdot(a \nabla u) = f)$ is also considered with initial condition $u\equiv 0$.
For further details consult the implementation in ```cookiepde.py```.
- forwardparabolic: the forward parabolic model with configurable inputs
- benchmarkparabolic: the forward parabolic model with fixed model inputs

## Properties
All four models take the same parametric input and return the same QoI (the integral of the solution over the spatial subdomain $F=[0.4,0.6]^2$).

Mapping | Dimensions	| Description
--------|-------------|------------
input |	[8] |	These values modify the conductivity coefficient in the 8 cookies, each of them must be greater than -1 (software does not check that input values are valid)
output |	[1] |	The integral of the solution over the central subdomain $F$

All four models only support the UM-BRDIGE evaluate feature.

Feature	| Supported
--------|---------
Evaluate|	True
Gradient|	False
ApplyJacobian|	False
ApplyHessian|	False

### Config
The basic configuration parameters are:

Model |Dictionary Key | Default value | User control | Description
------|---------------|---------------|--------------|-------------
Both | N | 400 | integer | The number of cells in each dimension (i.e. a mesh of $N^2$ elements).
Both | Fidelity | n/a | N = 100 * config['Fidelity'] | The 'fidelity' config key takes precedent over 'N' if both defined.
Both | BasisDegree | 1 | Integer | The degree of the piecewise polynomial FE approximation

The models take the further additional parameters.
These are not documented in detail.

Model |Dictionary Key | Default value | User control | Description
------|---------------|---------------|--------------|-------------
Both | advection | 0 | Integer | A flag to add an an advection field to the problem. If set to $1$ an advection term with advection field $\vec{w}=(4(x[1]-0.5)*(1-4(x[0]-0.5)^2), -4(x[0]-0.5)*(1-4(x[1]-0.5)^2))$ is added to the PDE problem.
Both | quad_degree | 8 | Integer | The quadrature degree used to evaluate integrals in the matrix assembly.
Both | diffzero | 1.0 | Float | Defines the background diffusion field a_0
Elliptic | directsolver | 1 | Integer | Uses a direct LU solve for linear system. Set to 0 to use GM-RES.
Elliptic | pc  | "none" | "ILU" or "JACOBI" | Preconditioning for the GM-RES solver
Elliptic | tol | 1e-4 | Float | Relative tolerance for the GM-RES solver.
Parabolic | letol  | 1e-4 | Float | Local timestepping error tolerance for a simple implementation of TR-AB2 timestepping.
Parabolic | T | 10.0 | Float | Final time for timestepping approximation. The QoI is evaluated and returned for time T.

## Mount directories
Mount directory | Purpose
---             |---
None            | 

## Description

![cookies-problem](https://raw.githubusercontent.com/UM-Bridge/benchmarks/main/models/cookies-problem/cookies_domain.png "geometry of the cookies problem")

The model implements the version of the cookies problem in [[Bäck et al.,2011]](https://doi.org/10.1007/978-3-642-15337-2_3), see also e.g. [[Ballani et al.,2015]](https://doi.org/10.1137/140960980), [[Kressner et al., 2011]](https://doi.org/10.1137/100799010) for slightly different versions. With reference to the computational domain $D=[0,1]^2$ in the figure above, the cookies model consists in the thermal diffusion problem below, where $\mathbf{y}$ are the uncertain parameters discussed in the following and $\mathrm{x}$ are physical coordinates 

$$-\mathrm{div}\Big[ a(\mathbf{x},\mathbf{y}) \nabla u(\mathbf{x},\mathbf{y}) \Big] = f(\mathrm{x}), \quad \mathbf{x}\in D$$

with homogeneous Dirichlet boundary conditions and forcing term defined as

$$f(\mathrm{x}) = \begin{cases} 
100 &\text{if } \,  \mathrm{x} \in F \\
0 &\text{otherwise} 
\end{cases}$$

where $F$ is the square $[0.4, 0.6]^2$. The 8 subdomains with uncertain diffusion coefficient (the cookies) are circles with radius 0.13 and the following center coordinates:

cookie | 1   | 2   | 3   | 4   | 5   | 6   | 7   | 8   |
--     | --  | --  | --  | --  | --  | --  | --  | --  |
x      | 0.2 | 0.5 | 0.8 | 0.2 | 0.8 | 0.2 | 0.5 | 0.8 |
y      | 0.2 | 0.2 | 0.2 | 0.5 | 0.5 | 0.8 | 0.8 | 0.8 |

The uncertain diffusion coefficient is defined as

$$a = 1 + \sum_{n=1}^8 y_n \chi_n(\mathrm{x})$$ 


where $y_n>-1$ and 

$$\chi_n(\mathrm{x}) = \begin{cases} 1 &\text{inside the n-th cookie} \\ 0 &\text{otherwise} \end{cases}$$


The output of the simulation is the integral of the solution over $F$, i.e. $\Psi = \int_F u(\mathrm{x}) d \mathrm{x}$

## Implementation details
### cookiepde.py
We consider a Python based approximation of the solution to the elliptic and parabolic ``cookie'' PDE problem.
Finite element approximation is performed on a quadrilateral grid implemented using FEniCS.
The linear systems are solved using PETSc via the petsc4py interface.
See the comments in the file for further implementation details.

### umbridge-server.py
This defines the UM-BRIDGE interfaces for the test problem.
Four models are defined
- forward: the forward elliptic PDE problem. The discrete formulation is constructed via FEniCS and solved using the PETSC KSP linear solver.
- benchmark: the forward elliptic PDE problem with fixed configuration.
- forwardparabolic: the forward parabolic PDE problem. This uses a custom implementation of the TR-AB2 adaptive timestepping algorithm to solve a parabolic formulation of the problem. We use a initial condition $u\equiv 0$.
- forwardparabolic: the forward parabolic PDE problem with fixed configuration 

### test_output.py
Runs the model for three different parameter/configuration pairs and verifies against a precomputed approximation.
Used for regression testing.

### test_output_parabolic.py
The equivalent set of sets for the parabolic problem.
Note that the solutions are almost identical to the elliptic case --- the steady state solution (long time limit) is the same as the elliptic problem.

### test_output_template.py
Python script to generate QoI approximations for testing. This evaluates the quantity of interest for three different parameters at four different fidelities, with varied quadrature degree.
Results can be piped from the console for plotting
```
python3 test_output_template.py http://0.0.0.0:4242 >> results.txt
```
This does not verify results against precomputed values.