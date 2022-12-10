# Tritium Diffusion

## Overview
In this benchmark, we use [Achlys](https://github.com/aurora-multiphysics/achlys) to model the macroscopic transport of tritium through fusion reactor materials using the Foster-McNabb equations. Achlys is built on top of the  [MOOSE Finite Element Framework](https://mooseframework.inl.gov/). The aim is to construct a functional emulator of the tritium desorption curve from uncertain input parameters.

## Authors
[Mikkel Lykkegaard](mailto:mikkel@digilab.co.uk)


## Run
```
docker run -it -p 4242:4242 linusseelinger/model-achlys:latest
```

## Properties

Model | Description
---|---
forward | Achlys Tritium Diffusion

### forward
Mapping | Dimensions | Description
---|---|---
input | [5] | E1, E2, E3: The detrapping energy of the traps. n1, n2: The density of the intrinsic traps.
output | [500] | Volumetric flux of tritium across the boundary as a function of time.

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

[Model sources here.](https://github.com/UM-Bridge/benchmarks/tree/main/models/achlys)

## Description
Achlys models macroscopic tritium transport processes through fusion reactor materials as described in the [Achlys documentation](https://aurora-multiphysics.github.io/achlys/module/introduction.html) and in [Delaporte-Mathurin et al (2019)](https://www.sciencedirect.com/science/article/pii/S2352179119300547).

The setup used in this particular benchmark models the experimental work of [Ogorodnikova et al. (2003)](https://www.sciencedirect.com/science/article/abs/pii/S0022311502013752)

## Benchmark
1. The model is evaluated with QMC samples in the following input space:
    - E1: $\mathcal U(0.7, 1.0)$
    - E2: $\mathcal U(0.9, 1.3)$
    - E3: $\mathcal U(1.1, 1.75)$
    - n1: $\mathcal U(5 \cdot 10^{-4}, 5 \cdot 10^{-3})$
    - n2: $\mathcal U(10^{-4}, 10^{-3})$

2. An emulator of the input-to-output map is constructed using a functional GP. The input to the emulator is the five-dimensional parameter space described above, while the output is the corresponding tritium desorption curve, including the predictive uncertainty.

## References
- Rémi Delaporte-Mathurin, Etienne A. Hodille, Jonathan Mougenot, Yann Charles, Christian Grisolia, *Finite element analysis of hydrogen retention in ITER plasma facing components using FESTIM*, Nuclear Materials and Energy, Volume 21, 2019

- O.V Ogorodnikova, J Roth, M Mayer, *Deuterium retention in tungsten in dependence of the surface conditions*, Journal of Nuclear Materials, Volumes 313–316, 2003
