# Composite material with random wrinkle

## Overview
This model implements the 3D anisotropic linear elasticity equations for a composite part with randomised wrinkle.


## Authors
[Anne Reinarz](mailto:anne.k.reinarz@durham.ac.uk)


## Run
```
docker run -it -p 4242:4242 linusseelinger/model-dune-composites:latest
```

## Properties

Model | Description
---|---
forward | Linear elasticity

### forward
Mapping | Dimensions | Description
---|---|---
input | [346000] | Coefficients of a Karhunen Loeve expansion of a typical wrinkle
output | [1] | Maximum deflection of the composite part

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
ranks | int | 2 | Number of MPI ranks (i.e. parallel processes) to be used
stack | string | "example2.csv" | Path to the stacking sequence to be run

## Mount directories
Mount directory | Purpose
---|---
None |

## Source code

[Model sources here.](https://github.com/UM-Bridge/benchmarks/tree/main/models/dune-composites)

## Description
In the simulation the composite strength of a corner part with a wrinkle is modeled. The analysis assumes standard anisotropic 3D linear elasticity.

## References
- Anne Reinarz, Tim Dodwell, Tim Fletcher, Linus Seelinger, Richard Butler, Robert Scheichl, *Dune-composites – A new framework for high-performance finite element modelling of laminates*, Composite Structures, 2018.

- Anhad Sandhu and Anne Reinarz and Tim J. Dodwell, *A Bayesian framework for assessing the strength distribution of composite structures with random defects*, Composite Structures, 2018.
}
