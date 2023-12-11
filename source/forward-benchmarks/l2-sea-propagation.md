# L2-Sea propagation

## Overview
This propagation benchmark is based on the L2-Sea model. The goal is to propagate operational uncertainties to the total resistance.

![L2-Sea-Model](https://raw.githubusercontent.com/UM-Bridge/benchmarks/main/models/l2-sea/l2sea_example.png "DTMB 5415 view of the wave elevation pattern and pressure field on the hull surface")

## Authors
- [Andrea Serani](mailto:andrea.serani@cnr.it)
- [Lorenzo Tamellini](mailto:lorenzo.tamellini@cnr.it)
- [Riccardo Pellegrini](mailto:riccardo.pellegrini@cnr.it)
- [Matteo Diez](mailto:matteo.diez@cnr.it)

## Run
```
docker run -it -p 4242:4242 linusseelinger/model-l2-sea
```

## Properties

Model | Description
---|---
forward | l2-sea

### forward
Mapping | Dimensions | Description
---|---|---
input | [16] | The first input is the Froude number (from 0.25 to 0.41); the second is the draft (from -6.776 to -5.544); the other 14 are the $x$ design variables for the shape modification with $-1\leq x_i \leq 1$ for $i=1,\dots,14$\\.
output | [5] | The first output is the model scale total resistance ($R_\mathrm{T}$) in Newton, whereas the other four are geometrical constraints (negative to be satisfied), related to the beam, draft, and sonar dome dimensions.

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


The benchmark pertains to the evaluation of the expected value and standard deviation of the DTMB 5415 model scale total resistance in calm water, conditional to operational and geometrical uncertain parameters.

The standard benchmark is a bi-dimensional problem considering the operational uncertainties only, the speed and the payload. The latter is associated with the hull draft. Speed is expressed by its non-dimensional counterpart, the Froude number (Fr), ranging between 0.25 and 0.41 with a triangular distribution, with the maximum at Fr=0.25 and equal to zero at Fr=0.41. For the draft, a beta distribution is defined as follows ...

For UQ method scalability testing, the benchmark can be extended in dimensionality (up to 16 dimensions) by adding the geometrical uncertainties associated with the shape modification (14 design variables are available). Also for the geometrical uncertainties, a beta distribution is used.

