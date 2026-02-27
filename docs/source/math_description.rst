.. _math-description:

=====================================
Mathematical abstraction in UM-Bridge
=====================================

In this section, we will describe UM-Bridge's interface mathematically. Note that both inputs and 
ouputs are required to be list of lists in actual implementation, but we only consider a single 
element within the outer list to simply notation.

Let :math:`F` denote the numerical model that maps the model input vector, :math:`\mathbf{\theta}` 
to the output vector :math:`\mathbf{F}(\mathbf{\theta})`:

.. math::    
    F\, : \,
    \mathbb{R}^n
    \;\longrightarrow\;
    \mathbb{R}^m.

Additionally, there may be an objective function :math:`L = L(\mathbf{F}(\mathbf{\theta}))`. UM-Bridge
allows the following four operations.

Model Evaluation
================
This is simply the so called forward map that takes an input 
:math:`\mathbf{\theta} = (\theta_1, \ldots, \theta_n) \in \mathbb{R}^n` and returns the output 
:math:`\mathbf{F}(\mathbf{\theta}) = (F(\mathbf{\theta})_1, \ldots, F(\mathbf{\theta})_m) \in \mathbb{R}^m`.


Gradient Evaluation
===================

The ``gradient`` function evaluates the sensitivity of a scalar 
objective. Using the chain rule:

.. math::
    \nabla_{\mathbf{\theta}}L
    = \left(\frac{\partial \mathbf{f}}{\partial \mathbf{\theta}}\right)^{\!\top}
    \boldsymbol{\lambda},
    \qquad
    \boldsymbol{\lambda} = \frac{\partial L}{\partial \mathbf{f}},

where :math:`\mathbf{\lambda}` is known as the sensitivity vector.


Applying Jacobian
=================

The ``apply_jacobian`` function evaluates the product of the model's Jacobian, :math:`J`, and a
vector, :math:`\mathbf{v}`, of the user's choice. The Jacobian of a vector-valued function 
is given by
 
.. math::
    J =
    \frac{\partial \mathbf{f}}{\partial \mathbf{\theta}} =
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
    = \frac{\partial \mathbf{f}}{\partial \mathbf{\theta}}\,\mathbf{v}.

Additionally, we can use this (or vice versa) to express the ``gradient`` function by setting 
:math:`\mathbf{v} = \mathbf{\lambda}`.  


Applying Hessian
================

This is a combination of the previous two sections: the output is still a matrix-vector product, but 
the matrix is the Hessian of an objective function. The Hessian, :math:`H`, is given by

.. math::
    H =
    \frac{\partial^2 L}{\partial \mathbf{\theta}\,\partial \mathbf{\theta}}
    = \frac{\partial}{\partial \mathbf{\theta}}
    \left(
    \frac{\partial \mathbf{f}}{\partial \mathbf{\theta}}
    \right)^{\!\top}
    \boldsymbol{\lambda} = 
    \begin{bmatrix}
    \dfrac{\partial^2 L}{\partial x_1^2} & \dfrac{\partial^2 L}{\partial x_1 \partial x_2} & \cdots & \dfrac{\partial^2 L}{\partial x_1 \partial x_n} \\[18pt]
    \dfrac{\partial^2 L}{\partial x_2 \partial x_1} & \dfrac{\partial^2 L}{\partial x_2^2} & \cdots & \dfrac{\partial^2 L}{\partial x_2 \partial x_n} \\[18pt]
    \vdots & \vdots & \ddots & \vdots \\[6pt]
    \dfrac{\partial^2 L}{\partial x_n \partial x_1} & \dfrac{\partial^2 L}{\partial x_n \partial x_2} & \cdots & \dfrac{\partial^2 L}{\partial x_n^2}
    \end{bmatrix},

where :math:`L` is the objective function and :math:`\mathbf{\lambda}` is the sensitivity vector as defined in the ``gradient`` 
section.

So the output for a chosen vector can be written as

.. math::
    H\,\mathbf{v}
    = \frac{\partial^2 \mathcal{L}}{\partial \mathbf{\theta}\,\partial \mathbf{\theta}}\,\mathbf{v} = 
    \left[\frac{\partial}{\partial \mathbf{\theta}}
    \left(
    \frac{\partial \mathbf{f}}{\partial \mathbf{\theta}}
    \right)^{\!\top}
    \boldsymbol{\lambda}\right]\,\mathbf{v}.
