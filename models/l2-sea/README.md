# L2-Sea

## Overview
The model provides the calm-water total resistance of a destroyer-type vessel as a function of the advancing speed (Froude number) and up to 14 design variables for the shape modification. The parent vessel under investigation is the DTMB 5415, an hull-form widley used for towing tank experiments, computational fluid dynamics studies, and shape optimization. The model is available as an open-source fortran code for the solution of potential flow equations at CNR-INM, MAO Research Group repository github.com/MAORG-CNR-INM/NATO-AVT-331-L2-Sea-Benchmark.   

![L2-Sea-Model](https:/github.com/UM-Bridge/benchmarks/edit/main/models/l2-sea/figs/l2-sea_example.png "DTMB 5415 view of the wave elevation pattern and pressure field on the hull surface")

## Authors
- [Andrea Serani](mailto:andrea.serani@cnr.it)
- [Lorenzo Tamellini](mailto:lorenzo.tamellini@cnr.it)
- [Riccardo Pellegrini](mailto:riccardo.pellegrini@cnr.it)
- [Matteo Diez](mailto:matteo.diez@cnr.it)

## Run
```
docker run -it -p 4242:4242 -v ~/l2-sea_output:/output linusseelinger/model-l2-sea 
```

## Properties

Model | Description
---|---
forward | l2-sea
optimization | l2-sea

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
fidelity | integer | 1 | Fidelity level for the total resistance evaluation associated to the numerical grid discretization. Fidelity goes from 1 to 7, where 1 is highest-fidelity level (finest grid) and 7 is the lowest-fidelity level (coarsest grid).
sinkoff | character | 'y' | Enabling hydrodynamics coupling with the rigid-body equation of motions for the ship sinkage. 'n' enables, 'y' disables.
trimoff | character | 'y' | Enabling hydrodynamics coupling with the rigid-body equation of motions for the ship trim. 'n' enables, 'y' disables.

## Mount directories
Mount directory | Purpose
---|---
/output | \texttt{ASCII} files for visualization of pressure distribution along the hull \texttt{pre\textit{XXXX}.plt} and free-surface \texttt{intfr\textit{XXXX}.plt} formatted for Tecplot and Paraview, where \texttt{\textit{XXXX}} is the Froude number.

## Source code

[Model sources here.](https://github.com/UM-Bridge/benchmarks/tree/main/models/l2-sea)

## Description
This model describes the calm-water resistance of a destroyer-type vessel by potential flow. Specifically, the vessel under investigation is the DTMB 5415 (at model scale), which is a widely used benchmark for towing tank experiments, CFD studies, and hull-form optimization, considering both deterministic and stochastic formulations.

Potential flow solver is used to evaluate the hydrodynamic loads, based on the Laplacian equation
%
\begin{equation}
    \nabla^2\phi = 0
\end{equation}
%
where $\phi$ is the velocity scalar potential, satisfying $\mathbf{u}=\nabla\phi$ and $\mathbf{u}$ is the flow velocity vector. The velocity potential $\phi$ is evaluated numerically through the Dawson linearization of the potential flow equations, using the boundary element method. Finally, the total resistance is estimated as the sum of the wave and the frictional resistance: the wave resistance component is estimated by integrating the pressure distribution over the hull surface, obtained using the Bernoulli theorem
%
\begin{equation}
    \frac{p}{\rho} + \frac{\left(\nabla\phi\right)^2}{2}-gz = cost;
\end{equation}
%
the frictional resistance component is estimated using a flat-plate approximation based on the local Reynolds number.

The steady 2 degrees of freedom (sinkage and trim) equilibrium is achieved considering iteratively the coupling between the hydrodynamic loads and the rigid-body equation of motion. 

The model can exploit multiple grid discretization levels, whose details can be found in [Pellegrini et al. (2022)](https://dl.acm.org/doi/10.1145/3458817.3476150](https://www.mdpi.com/2227-7390/10/3/481).

### forward UQ problem

The benchmark pertains to the evaluation of the expected value and standard deviation of the DTMB 5415 model scale total resistance in calm water, conditional to operational and geometrical uncertain parameters.

The standard benchmark is a bi-dimensional problem considering the operational uncertainties only, the speed and the payload. The latter is associated with the hull draft. Speed is expressed by its non-dimensional counterpart, the Froude number (Fr), ranging between 0.25 and 0.41 with a triangular distribution, with the maximum at Fr=0.25 and equal to zero at Fr=0.41. For the draft, a beta distribution is defined as follows ...

For UQ method scalability testing, the benchmark can be extended in dimensionality (up to 16 dimensions) by adding the geometrical uncertainties associated with the shape modification (14 design variables are available). Also for the geometrical uncertainties, a beta distribution is used.

Mapping | Dimensions | Description
---|---|---
input | [16] | The first input is the Froude number (from 0.25 to 0.41); the second is the draft (from -6.776 to -5.544); the other 14 are the $x$ design variables for the shape modification with $-1\leq x_i \leq 1$ for $i=1,\dots,14$\\.
output | [5] | The first output is the model scale total resistance ($R_\mathrm{T}$) in Newton, whereas the other four are geometrical constraints (negative to be satisfied), related to the beam, draft, and sonar dome dimensions.

### shape optimization problem

The benchmark, developed within the activities of the NATO-AVT-331 Research Task Group on ``Goal-driven, multifidelity approaches for military vehicle system-level design'' [Beran et al. (2020)](https://arc.aiaa.org/doi/abs/10.2514/6.2020-3158), pertains to the total resistance reduction of the DTMB 5415 in calm water at fixed speed, corresponding to a Froude number (Fr) equal to 0.28. The optimization problem reads
%
\begin{eqnarray}\label{eq:5415prob}
    \begin{array}{rll}
        \mathrm{minimize}      & \Delta R (\mathbf{x}) = \frac{R (\mathbf{x})}{R_0}-1  \qquad \mathrm{with} \qquad \mathbf{x}\in\mathbb{R}^N\\
        \mathrm{subject \,  to}& L_{\rm pp}(\mathbf{x}) = L_{\rm pp_0}\\
        \mathrm{and \, to}     & \nabla(\mathbf{x}) = \nabla_0 \\
        & |\Delta B(\mathbf{x})| \leq 0.05B_0 \\
        & |\Delta T(\mathbf{x})| \leq 0.05T_0 \\
        & V(\mathbf{x})\geq V_0\\
        & -1\leq x_i \leq 1 \qquad \mathrm{with} \qquad \forall i=1,\dots, N\\
    \end{array}
\end{eqnarray}
%
where $\mathbf{x}$ are the design variables, $L_{\rm pp}$ is the length between perpendiculars, $B$ is the overall beam, $T$ is the drought, and $V$ is the volume reserved for the sonar in the bow dome. Subscript ``0'' indicates parent (original) hull values. Equality and inequality constraints for the geometry deformations are taken from \cite{grigoropoulos2017mission}. 

The shape modifications $\tilde{\boldsymbol{\gamma}}(\boldsymbol{\xi},\mathbf{x})$ are produced directly on the Cartesian coordinates $\boldsymbol{\xi}$ of the computational body surface grid $\mathbf{g}$, as per

\begin{equation}\label{eq:reducedspace}
    \mathbf{g}(\boldsymbol{\xi},\mathbf{x})=\mathbf{g}_0(\boldsymbol{\xi}) + \boldsymbol{\gamma}(\boldsymbol{\xi},\mathbf{x})
\end{equation}

where $\mathbf{g}_0$ is the original geometry and $\boldsymbol{\gamma}$ is a shape modification vector obtained by a physics-informed design-space dimensionality reduction \cite{serani2019stochastic}

\begin{equation}
    {\boldsymbol{\gamma}}(\boldsymbol{\xi},\mathbf{x}) = \sum_{k=1}^N x_k \boldsymbol{\psi}_k(\boldsymbol{\xi})
    \label{e:exp_gamma}
\end{equation}

with $\boldsymbol{\psi}$ a set of orthonormal functions, with $N=14$ the number of design variables ($\mathbf{x}$). It may be noted that the design variables and the associated shape modifications are organized in a hierarchical order, meaning that the first variables produce larger design modifications than the last ones \cite{serani2021hull}.   

The multifidelity levels are defined by the computational grid size. Specifically, the benchmark is defined with seven grid (fidelity) levels with a refinement ratio of 2$^{0.25}$

Mapping | Dimensions | Description
---|---|---
input | [16] | The first input is the Froude number (fixed to 0.28); the second is the draft (fixed to -6.16); the other 14 are the $x$ design variables for the shape modification with $-1\leq x_i \leq 1$  (with parent hull has $x_i = 0$) for $i=1,\dots,14$\\.
output | [5] | The first output is the model scale total resistance ($R_\mathrm{T}$) in Newton, whereas the other four are geometrical constraints (negative to be satisfied), related to the beam, draft, and sonar dome dimensions.
