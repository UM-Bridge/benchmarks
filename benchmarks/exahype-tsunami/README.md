# ExaHyPE-Tsunami Benchmark

In this benchmark we solve the basic shallow water equations with bathymetry source terms
(neglecting friction terms or more advanced models for with non-hydrostatic corrections).

The resulting equations can be written in first-order hyperbolic form as

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




