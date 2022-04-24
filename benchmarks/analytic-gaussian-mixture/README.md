# Analytic-Gaussian-Mixture Benchmark

## Authors
- [Linus Seelinger](mailto:linus.seelinger@iwr.uni-heidelberg.de)

## Overview
This benchmark consists of an analytically defined PDF $\pi : \mathbb{R}^2 \rightarrow \mathbb{R}$ consisting of a Gaussian mixture.

## Run
```
docker run -it -p 4243:4243 linusseelinger/benchmark-analytic-gaussian-mixture
```

## Properties
Value | Dimensions
---|---
inputSizes | [2]
outputSizes | [1]

Feature | Supported
---|---
Evaluate | True
Gradient | True
ApplyJacobian | True
ApplyHessian | False

### Configuration

None

### Description

- Input: 2D coordinates $x \in \mathbb{R}^2$
- Output: PDF $\tau$ evaluated at $x$

## Model

Let
$X_1 \sim \mathcal{N}(\begin{pmatrix} -1.5 \\ -1.5 \end{pmatrix}, 0.8 I)$,
$X_2 \sim \mathcal{N}(\begin{pmatrix} 1.5 \\ 1.5 \end{pmatrix}, 0.8 I)$,
$X_3 \sim \mathcal{N}(\begin{pmatrix} -2 \\ 2 \end{pmatrix}, 0.5 I)$.
Denote by $f_{X_1}, f_{X_2}, f_{X_3}$ the corresponding PDFs.

The PDF $\pi$ is then defined as

$$ \pi(x) := \sum_{i=1}^3 f_{X_i}(x), $$

and the benchmark outputs $\log(\pi(x))$.
