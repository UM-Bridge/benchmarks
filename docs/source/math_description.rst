.. _math-description:

=====================================
Mathematical abstraction in UM-Bridge
=====================================

UM-Bridge's interface can be described mathematically, which is what we will do here.

Model Evaluation
================
Let :math:`\mathcal{F}` denote the numerical model that maps the model input vector, :math:`\mathbf{x}` to
the output vector :math:`\mathbf{f(\mathbf{x})}`:

.. math::    
    \mathcal{F}\, : \,
    \mathbf{x}
    \;\longrightarrow\;
    \mathbf{f}(\mathbf{x}), \quad
    \mathbf{x} \in \mathbb{R}^d, \;
    \mathbf{f}(\mathbf{x}) \in \mathbb{R}^n.


Gradient Evaluation
===================

The ``gradient`` function evaluates the sensitivity of a scalar 
objective, :math:`L(\mathbf{f}(\mathbf{x}))`, that depends on the model output, with respect to the model input. Using the 
chain rule:

.. math::
    \nabla_{\mathbf{x}}L
    = \left(\frac{\partial \mathbf{f}}{\partial \mathbf{x}}\right)^{\!\top}
    \boldsymbol{\lambda},
    \qquad
    \boldsymbol{\lambda} = \frac{\partial L}{\partial \mathbf{f}},

where :math:`\lambda` is known as the sensitivity vector.


Applying Jacobian
=================

The ``apply_jacobian`` function evaluates the product of the model's Jacobian, :math:`J`, and a
vector, :math:`\mathbf{v}`, of the user's choice. The Jacobian of a vector-valued function 
is given by
 
.. math::
    J =
    \frac{\partial \mathbf{f}}{\partial \mathbf{x}} =
    \left[
    \begin{array}{ccc}
    \dfrac{\partial \mathbf{f}}{\partial x_1} & \cdots & \dfrac{\partial \mathbf{f}}{\partial x_d}
    \end{array}
    \right] = 
    \begin{pmatrix}
    \dfrac{\partial f_{1}}{\partial x_{1}} & \cdots &
    \dfrac{\partial f_{1}}{\partial x_{d}} \\[12pt]
    \vdots & \ddots & \vdots \\[4pt]
    \dfrac{\partial f_{n}}{\partial x_{1}} & \cdots &
    \dfrac{\partial f_{n}}{\partial x_{d}}
    \end{pmatrix}
    \in \mathbb{R}^{n \times d}.


The output of this function for a chosen :math:`\mathbf{v} \in \mathbb{R}^{d}` is then

.. math::
    \texttt{output}
    = J\,\mathbf{v}
    = \frac{\partial \mathbf{f}}{\partial \mathbf{x}}\,\mathbf{v}.

Additionally, we can use this (or vice versa) to expression the ``gradient`` function by setting 
:math:`\mathbf{v} = \mathbf{\lambda}`.  


Applying Hessian
================

This is a combination of the previous two sections: the output is still a matrix-vector product, but 
the matrix is the Hessian of an objective function. The Hessian, :math:`H`, is given by

.. math::
    H =
    \frac{\partial^2 L}{\partial \mathbf{x}\,\partial \mathbf{x}}
    = \frac{\partial}{\partial \mathbf{x}}
    \left(
        \frac{\partial \mathbf{f}}{\partial \mathbf{x}}
    \right)^{\!\top}
    \boldsymbol{\lambda},

where :math:`L` is the objective function and :math:`\mathbf{\lambda}` is the sensitivity vector as defined in the ``gradient`` 
section.

So the output for a chosen vector can be written as

.. math::
    H\,\mathbf{v}
    = \frac{\partial^2 \mathcal{L}}{\partial \mathbf{x}\,\partial \mathbf{x}}\,\mathbf{v} = 
    \left[\frac{\partial}{\partial \mathbf{x}}
    \left(
    \frac{\partial \mathbf{f}}{\partial \mathbf{x}}
    \right)^{\!\top}
    \boldsymbol{\lambda}\right]\,\mathbf{v}.
