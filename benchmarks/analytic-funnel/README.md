# Analytic-Funnel Benchmark

## Authors
- [Linus Seelinger](mailto:linus.seelinger@iwr.uni-heidelberg.de)

## Overview
This benchmark consists of an analytically defined PDF $\tau : \mathbb{R}^2 \rightarrow \mathbb{R}$ resembling the shape of a funnel.

## Run
```
docker run -it -p 4243:4243 linusseelinger/benchmark-analytic-funnel
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

First, define a helper function

$$ f(x,m,s) := - \frac12 \log(2 \pi) - \log(s) - \frac12 ((x-m)/s)^2. $$

Now, the output log PDF is defined as

$$ \log(\tau(x)) := f(x_1, 0, 3) + f(x_2, 0, \exp(\frac12 x_1)). $$
