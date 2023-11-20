# GS2 Fusion Plasma Simulation

## Overview
This model uses [GS2](https://gyrokinetics.gitlab.io/gs2/index.html) to simulate plasma in a spherical tokamak. It aims to study various modes in the plasma due to microinstabilities. The current setup terminates once an unstable mode is found. Otherwise, it will continue until reaching a fixed timestep. Therefore, the runtime varies depending on the input parameters.

## Authors
- [Chung Ming Loi](mailto:chung.m.loi@durham.ac.uk)


## Run
```
docker run -it -p 4242:4242 linusseelinger/model-achlys:latest
```

## Properties

Model | Description
---|---
forward | Plasma simulation

### forward
Mapping | Dimensions | Description
---|---|---
input | [2] | [\texttt{tprim}: normalised inverse temperature gradient, \texttt{vnewk}: normalised species-species collisionality frequency]. Both are set for electrons only. 
output | [3] | [Electron heat flux, electric field growth rate, electric field mode frequency]

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

[Model sources here.](https://github.com/UM-Bridge/benchmarks/tree/gs2/models/gs2)

## Description

