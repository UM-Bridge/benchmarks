# ExaHyPE-Tsunami Benchmark

### Brief intro
In this benchmark we model the propagation of the 2011 Tohoku tsunami by solving the shallow water equations. For the numerical solution of the PDE, we apply an ADER-DG method implemented in the [ExaHyPE framework](https://www.sciencedirect.com/science/article/pii/S001046552030076X). The aim is to obtain the parameters describing the initial displacements from the data of two available buoys located near the Japanese coast

## Purpose 
Demonstrate the parallelized multi-level Markov Chain Monte Carlo (MLMCMC) method on a practically relevant large-scale application

## Run
```
docker run -it -p 4243:4243 linusseelinger/model-exahype-tsunami
```

## Technical properties

- Input size: 2
- Output size: 4
- Config JSON structure:
    - `int level`
    - `bool verbose`
    - `bool vtk_output`
- Supported endpoints (evaluate, jacobian etc)

## Detailed description on technical properties

- Input consists of x and y coordinates of a proposed tsunami origin
- Output consists of time and maximum water height at two buoy points
- JSON Configuration:
    - `level`: chooses the model level to run (see below for further details)
    - `verbose`: switches text output on/off
    - `vtk_output`: switches vtk output on/off

## Forward model
The underlying PDE model can be written in first-order hyperbolic form as
$$
    \frac{\partial}{\partial t}
    \begin{pmatrix}
    h\\hu\\hv\\ b
    \end{pmatrix} + \nabla \cdot
    \begin{pmatrix}
    hu   &   hv\\
    hu^2 & huv\\
    huv & hv^2 \\
    0 & 0\\
    \end{pmatrix}+
    \begin{pmatrix}
    0\\
    hg \, \partial_x (b+h)\\
    hg \, \partial_y (b+h)\\
    0\\
    \end{pmatrix}= 0,
$$

where 
- $h$ denotes the height of the water column, 
- $(u,v)$ the horizontal flow velocity, 
- $g$  gravity 
- $b$ denotes the bathymetry.

This benchmark creates a sequence of three models:
1. First model: 
    - bathymetry is approximated only by a depth average over the entire domain
    - pure DG discretisation of order 2
2. The second model:
    - DG discretisation with a finite volume subcell limiter allowing for wetting and drying
    - smoothed bathymetry data (Gaussian filter)
3. The third models:
    - DG discretisation with a finite volume subcell limiter allowing for wetting and drying
    - full bathymetry data

- The bathymetry data has been obtained from [GEBCO](https://www.gebco.net/)
- More details: [Reference Paper](https://dl.acm.org/doi/10.1145/3458817.3476150)

##  UQ problem

The likelihood of a given set of parameters given the simulation results is computed using weighted average of the maximal wave height and the time at which it is reached. The likelihood is given by a normal distribution $\mathcal{N}\left(\mu, \Sigma \right)$ with mean $\mu$ given by maximum waveheight $\max\{h\}$ and the time $t$ at which it is reached for the the two DART buoys 21418 and 21419 (This data can be obtained from [NDBC](https://www.ndbc.noaa.gov/). The covariance matrix $\Sigma$ depends on the level, but not the probe point.

| $\mu$   | $\Sigma$ l=0 |  $\Sigma$ l=1 |  $\Sigma$ l=2 |
|---------|--------------|---------------|---------------|
| 1.85232 | 0.15         | 0.1           | 0.1           |
| 0.6368  | 0.15         | 0.1           | 0.1           |
| 30.23   | 2.5          | 1.5           | 0.75          |
| 87.98   | 2.5          | 1.5           | 0.75          |

The prior cuts off all parameters which would lead to an initial displacement which is too close to the domain boundary. Some parameters may lead to unstable models, e.g. a parameter which initialise the tsunami on dry land, in this case we have treated the parameter as unphysical and assigned an almost zero likelihood.

The paralle MLMCMC was implemented in the [MUQ library](https://joss.theoj.org/papers/10.21105/joss.03076).
