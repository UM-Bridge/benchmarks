.. _math-description:

=====================================
Mathematical abstraction in UM-Bridge
=====================================

In this section, we will describe UM-Bridge's interface mathematically. Note that both inputs and 
ouputs are required to be a list of lists in the actual implementation, but we only consider a single 
element within the outer list to simply the notation from hereon.

Let :math:`F` denote the numerical model that maps the model input vector, :math:`\mathbf{\theta}` 
to the output vector :math:`\mathbf{F}(\mathbf{\theta})`:

.. math::    
    F\, : \,
    \mathbb{R}^n
    \;\longrightarrow\;
    \mathbb{R}^m.

Additionally, there may be an objective function :math:`L = L(\mathbf{F}(\mathbf{\theta}))`. 

UM-Bridge allows the following four operations.

Model Evaluation
================
This is simply the so called forward map that takes an input 
:math:`\mathbf{\theta} = (\theta_1, \ldots, \theta_n) \in \mathbb{R}^n` and returns the model output 
:math:`\mathbf{F}(\mathbf{\theta}) = (F(\mathbf{\theta})_1, \ldots, F(\mathbf{\theta})_m) \in \mathbb{R}^m`.


Gradient of the objective
=========================

The gradient function evaluates the sensitivity of the scalar objective. Using the chain rule:

.. math::
    :name: eq:1
    
    \nabla_{\mathbf{\theta}}L
    = \left(\frac{\partial \mathbf{F}}{\partial \mathbf{\theta}}\right)^{\!\top}
    \boldsymbol{\lambda},
    \qquad
    \boldsymbol{\lambda} = \frac{\partial L}{\partial \mathbf{F}},

where :math:`\mathbf{\lambda}` is known as the sensitivity vector.

Most UQ algorithms do not evaluate the full gradient vector but rather select a specific
component within the input (:math:`\theta_i`) and output vectors (:math:`F_j`). These indices are
chosen using ``inWrt`` and ``outWrt``, respectively, in the implementation. So :ref:`(1) <eq:1>` becomes

.. math::
    
    \frac{\partial L}{\partial \theta_i}
    = \frac{\partial F_j}{\partial \theta_i}
    \lambda_j,
    \qquad
    \lambda_j = \frac{\partial L}{\partial F_j},
    
where :math:`\lambda_j` is the ``sens`` argument in the code. 

Applying Jacobian to a vector
=============================

The apply Jacobian function evaluates the product of the model's Jacobian, :math:`J`, and a
vector, :math:`\mathbf{v}`, of the user's choice (``vec``). The Jacobian of a vector-valued function 
is given by
 
.. math::
    J =
    \frac{\partial \mathbf{F}}{\partial \mathbf{\theta}} =
    \left[
    \begin{array}{ccc}
    \dfrac{\partial \mathbf{F}}{\partial \theta_1} & \cdots & \dfrac{\partial \mathbf{F}}{\partial \theta_n}
    \end{array}
    \right] = 
    \begin{pmatrix}
    \dfrac{\partial F_{1}}{\partial \theta_{1}} & \cdots &
    \dfrac{\partial F_{1}}{\partial \theta_{n}} \\[12pt]
    \vdots & \ddots & \vdots \\[4pt]
    \dfrac{\partial F_{n}}{\partial \theta_{1}} & \cdots &
    \dfrac{\partial F_{n}}{\partial \theta_{n}}
    \end{pmatrix}
    \in \mathbb{R}^{m \times n}.


For a chosen :math:`\mathbf{v} \in \mathbb{R}^{n}`, this is simply

.. math::
    J\,\mathbf{v}
    = \frac{\partial \mathbf{F}}{\partial \mathbf{\theta}}\,\mathbf{v}.

Additionally, we can use this to express the gradient function by setting 
:math:`\mathbf{v} = \mathbf{\lambda}`.  

However, we don't actually assemble the full Jacobian. We apply specific indices of the Jacobian, 
:math:`J_{ji} = \frac{\partial F_j}{\partial \theta_i}`, to the vector instead. The output of this 
action is then

.. math::
    \texttt{output} =
    J_{ji}\,\mathbf{v}
    = \frac{\partial F_j}{\partial \theta_i}\,\mathbf{v},

where the the :math:`i^{th}` and :math:`j^{th}` indices coresspond to ``inWrt`` and ``outWrt``.

Applying Hessian to a vector
============================

The apply Hessian action is a combination of the previous two sections: the action is still a matrix-vector product, but 
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

where :math:`L` is the objective function and :math:`\mathbf{\lambda}` is the sensitivity vector as defined previously.

So the product of :math:`H` and the chosen vector can be written as

.. math::
    H\,\mathbf{v}
    = \frac{\partial^2 \mathcal{L}}{\partial \mathbf{\theta}\,\partial \mathbf{\theta}}\,\mathbf{v} = 
    \left[\frac{\partial}{\partial \mathbf{\theta}}
    \left(
    \frac{\partial \mathbf{f}}{\partial \mathbf{\theta}}
    \right)^{\!\top}
    \boldsymbol{\lambda}\right]\,\mathbf{v}.

Again, we don't evaluate the full Hessian in UM-Bridge. As in the apply Jacobian action, we select certain indices and
apply them the vector. Since :math:`H` contains the seconds derivative of :math:`L`, we require two indices for the input:
:math:`inWrt1` and :math:`inWrt2`. The output of this action is 

.. math::
    \texttt{output} =
    \left( \frac{\partial}{\partial \theta_i}
    \left[ \frac{\partial F_k}{\partial x_j} \, \lambda_k \right]
    \, \mathbf{v} \right).


