# L2-Sea

## Overview
The model provides the calm-water total resistance of a destroyer-type vessel as a function of the advancing speed (Froude number) and up to 14 design variables for the shape modification. The parent vessel under investigation is the DTMB 5415, an hull-form widley used for towing tank experiments, computational fluid dynamics studies, and shape optimization. The model is available as an open-source fortran code for the solution of potential flow equations at CNR-INM, MAO Research Group repository github.com/MAORG-CNR-INM/NATO-AVT-331-L2-Sea-Benchmark.   

## Authors
Andrea Serani, Lorenzo Tamellini, Riccardo Pellegrini, Matteo Diez

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
input | [15] | Froude number ($0.25\leq \mathrm{Fr}\leq 0.41$) and design variables vector for the shape modification composed by $N=14$ design varibales $-1\leq x_i \leq 1$ for {i=1,\dots,N}.
output | [5] | Total resistance in calm water ($R_\mathrm{T}$) and 4 geometrical constraints related to ... 

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
fidelity | integer | 7 | Fidelity level for the total resistance evaluation associated to the numerical grid discretization. Fidelity goes from 1 to 7, where 1 is highest-fidelity level (finest grid) and 7 is the lowest-fidelity level (coarsest grid).
sinkoff | character | 'y' | Enabling hydrodynamics coupling with the rigid-body equation of motions for the ship sinkage. 'n' enables, 'y' disables.
trimoff | character | 'y' | Enabling hydrodynamics coupling with the rigid-body equation of motions for the ship trim. 'n' enables, 'y' disables.

## Mount directories
Mount directory | Purpose
---|---
/output | ASCII output for visualization formatted for Tecplot and Paraview (*.plt)

## Source code

[Model sources here.](https://github.com/UM-Bridge/benchmarks/tree/main/models/l2-sea)

## Description
TODO
