# Tritium Diffusion

## Overview
We use [Achlys](https://github.com/aurora-multiphysics/achlys) to model the macroscopic transport of tritium through fusion reactor materials using the Foster-McNabb equations. Achlys is built on top of the  [MOOSE Finite Element Framework](https://mooseframework.inl.gov/).

## Authors
- [Mikkel Lykkegaard](mailto:mikkel@digilab.co.uk)
- [Anne Reinarz](mailto:anne.k.reinarz@durham.ac.uk)


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
output | [500] | Flux of tritium across the boundary as a function of time in atomic fraction.

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

Particularly, we solve the following equations:


```math
 \frac{\partial C_m}{\partial t} = \nabla \cdot \left( D \nabla C \right) - \sum_i \frac{\partial C_{t,i}}{\partial t} + \dot{S}_{\mathrm{ext}}
```
```math
\frac{\partial C_{t,i}}{\partial t} = \nu_m C_m (n_i - C_{t,i}) - \nu_i C_{t,i}
```
where $C_m$ is the concentration of mobile particles and $C_{t,i}$ is the concentration of particles at the $i$th trap type.

The flux of tritium is then assumed to be proportional to the tritium creation rate at additional extrinsic traps:
$$\frac{dN_{3}}{dt} = (1 - r) \phi \left[ \left(1-\frac{N_3}{n_{3a,max}}\right)\eta_a f(x) + \left(1-\frac{N_3}{n_{3b,max}}\right)\eta_b\theta(x) \right]$$

Please see [Hodille et al. (2015)](https://www.sciencedirect.com/science/article/pii/S0022311515300660) for more details.

The setup used in this particular benchmark models the experimental work of [Ogorodnikova et al. (2003)](https://www.sciencedirect.com/science/article/abs/pii/S0022311502013752)

## References
- Rémi Delaporte-Mathurin, Etienne A. Hodille, Jonathan Mougenot, Yann Charles, Christian Grisolia, *Finite element analysis of hydrogen retention in ITER plasma facing components using FESTIM*, Nuclear Materials and Energy, Volume 21, 2019

- E.A. Hodille, X. Bonnin, R. Bisson, T. Angot, C.S. Becquart, J.M. Layet, C. Grisolia, *Macroscopic rate equation modeling of trapping/detrapping of hydrogen isotopes in tungsten materials*, Journal of Nuclear Materials, Volume 467, Part 1, 2015

- O.V Ogorodnikova, J Roth, M Mayer, *Deuterium retention in tungsten in dependence of the surface conditions*, Journal of Nuclear Materials, Volumes 313–316, 2003
