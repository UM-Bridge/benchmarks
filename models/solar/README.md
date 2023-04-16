# Solar

## Overview
The SOLAR blackbox optimization benchmarking framework.

## Authors
- [Sebastien Le Digibel](mailto:sebastien.le-digabel@polymtl.ca)
- [Anne Reinarz](mailto:anne.k.reinarz@durham.ac.uk)

## Run
```
docker run -it -p 4242:4242 annereinarz/solar
```

## Properties

Model | Description
---|---
forward | solar

### forward
Mapping | Dimensions | Description
---|---|---
input | [9] | Feasibility Point, at which the simulator is evaluated
output | [6] | Total solar energy on the receiver

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
[Model sources here.](https://github.com/UM-Bridge/TODO)

## Description

References
- Mathieu Lemyre Garneau, Modelling of a solar thermal power plant for benchmarking blackbox optimization solvers, 2015, [url](https://publications.polymtl.ca/1996/).
