# Analytic-Banana Benchmark

## Authors
- [Linus Seelinger](mailto:linus.seelinger@iwr.uni-heidelberg.de)

## Overview
This benchmark consists of an analytically defined PDF $\pi : \mathbb{R}^2 \rightarrow \mathbb{R}$ resembling the shape of a banana. It is based on a transformed normal distribution. The variance may be adjusted.

## Run
```
docker run -it -p 4243:4243 linusseelinger/benchmark-analytic-banana
```

## Properties
Value | Dimensions
---|---
inputSizes | [2]
outputSizes | [1]

Feature | Supported
---|---
Evaluate | True
Gradient | False
ApplyJacobian | False
ApplyHessian | False

### Configuration

- `double a` (default: 2.0)
- `double b` (default: 0.2)
- `double scale` (default: 1.0)

### Description

- Input: 2D coordinates $x \in \mathbb{R}^2$
- Output: PDF $\pi$ evaluated at $x$
- JSON Configuration:
    - `scale`: Scaling factor applied to the underlying normal distribution's variance
    - `a`: Transformation parameter
    - `b`: Transformation parameter

## Model

We begin with a normally distributed random variable $Z \sim \mathcal{N}(\begin{pmatrix} 0 \\ 4 \end{pmatrix}, scale \begin{pmatrix} 1.0 & 0.5\\ 0.5 & 1.0 \end{pmatrix})$, and denote its PDF by $f_Z$.

In order to the normal distribution, define a transformation $T : \mathbb{R}^2 \rightarrow \mathbb{R}^2$

$$ T(x) := \begin{pmatrix} x_1 / a \\ a x_2 + a b (x_1^2 + a^2) \end{pmatrix}. $$

Finally, the benchmark log PDF is defined as

$$ log(\pi(x)) := log(f_Z(T(x))). $$
