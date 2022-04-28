# Analytic-Donut Benchmark

## Overview
This benchmark consists of an analytically defined PDF $\pi : \mathbb{R}^2 \rightarrow \mathbb{R}$ resembling the shape of a donut.

## Authors
- [Linus Seelinger](mailto:linus.seelinger@iwr.uni-heidelberg.de)

## Run
```
docker run -it -p 4243:4243 linusseelinger/benchmark-analytic-donut
```

## Properties

Mapping | Dimensions | Description
---|---|---
input | [2] | 2D coordinates $x \in \mathbb{R}^2$
output | [1] | PDF $\pi$ evaluated at $x$

Feature | Supported
---|---
Evaluate | True
Gradient | True
ApplyJacobian | True
ApplyHessian | False

Config | Type | Default | Description
---|---|---|---
None | | |

Mount directory | Purpose
---|---
None |

## Description

The PDF $\pi$ is defined as

$$ \pi(x) := - \frac{(\| x \| - r)^2}{\sigma^2}, $$

where $r = 2.6$ and $\sigma^2 = 0.033$.

The implementation then returns the log PDF $\log(\pi(x))$.
