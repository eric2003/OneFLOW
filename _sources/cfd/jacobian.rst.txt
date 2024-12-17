Jacobian matrix
==================================

In vector calculus, the Jacobian matrix of a vector-valued function of several variables is the matrix of all its first-order partial derivatives. When this matrix is square, that is, when the function takes the same number of variables as input as the number of vector components of its output, its determinant is referred to as the Jacobian determinant. Both the matrix and (if applicable) the determinant are often referred to simply as the Jacobian in literature.

Reference
---------------------------------------------------------------------------
#. `Jacobian matrix and determinant <https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant/>`_
#. `A Gentle Introduction to the Jacobian <https://machinelearningmastery.com/a-gentle-introduction-to-the-jacobian/>`_
#. `Jacobian <https://mathworld.wolfram.com/Jacobian.html>`_
#. `Change of Variables in Multiple Integrals (Jacobians) <https://math.libretexts.org/Courses/Monroe_Community_College/MTH_212_Calculus_III/Chapter_14%3A_Multiple_Integration/14.7%3A_Change_of_Variables_in_Multiple_Integrals_(Jacobians)/>`_


The Jacobian Matrix
---------------------------------------------------------------------------
The Jacobian matrix collects all first-order partial derivatives of a multivariate function.
Specifically, consider first a function that maps :math:`n` real inputs, to a single real output:

.. math::
  \begin{align}
  f&:{R}^n\to {R}\\
  f&=f(x_{1},x_{2}, \cdots ,x_{n})
  \end{align}
  
Then, for an input vector, :math:`\mathbf{x}`, of length, :math:`n`, the Jacobian vector of size, :math:`1 × n`, can be defined as follows:

.. math::
  \begin{align}
  J&=\cfrac{\mathrm{d} f(\mathbf{x})}{\mathrm{d} \mathbf{x}} =\cfrac{\mathrm{d}f(x_{1},x_{2}, \cdots ,x_{n})}{\mathrm{d}(x_{1},x_{2}, \cdots ,x_{n})} \\
  &=\begin{bmatrix}
    \cfrac{\partial f(x_{1},x_{2}, \cdots ,x_{n})}{\partial x_{1}},
    &\cfrac{\partial f(x_{1},x_{2}, \cdots ,x_{n})}{\partial x_{2}},
    & \cdots,
    &\cfrac{\partial f(x_{1},x_{2}, \cdots ,x_{n})}{\partial x_{n}}
  \end{bmatrix}
  \end{align}
  
This formula can also be written as
 
.. math::
  \begin{align}
  J&=\cfrac{\mathrm{d} f(\mathbf{x})}{\mathrm{d} \mathbf{x}} =\cfrac{\mathrm{d}f(x_{1},x_{2}, \cdots ,x_{n})}{\mathrm{d}(x_{1},x_{2}, \cdots ,x_{n})} \\
  &=\begin{bmatrix}
    \cfrac{\partial f(\mathbf{x})}{\partial x_{1}},
    &\cfrac{\partial f(\mathbf{x})}{\partial x_{2}},
    & \cdots,
    &\cfrac{\partial f(\mathbf{x})}{\partial x_{n}}
  \end{bmatrix}
  \end{align}
  
Now, consider another function that maps :math:`n` real inputs, to :math:`m` real outputs:

.. math::
  \begin{align}
   \mathbf{f}&:{R}^n\to {R}^m\\
   \mathbf{J}&\in R^{m\times n}\\
   \mathbf{f}&=
      \begin{bmatrix}
       f_{1}(x_{1},x_{2}, \cdots ,x_{n})\\
       f_{2}(x_{1},x_{2}, \cdots ,x_{n})\\
       \cdots \\
       f_{m}(x_{1},x_{2}, \cdots ,x_{n})
      \end{bmatrix}
  \end{align}


Then, for the same input vector, :math:`\mathbf{x}`, of length, :math:`n`, the Jacobian is now a :math:`m × n` matrix, :math:`\mathbf{J}\in R^{m\times n}`, that is defined as follows:

.. math::
  \begin{align}
  \mathbf{J}&=\cfrac{\mathrm{d} \mathbf{f}(\mathbf{x})}{\mathrm{d} \mathbf{x}} =\cfrac{\mathrm{d}\mathbf{f}(x_{1},x_{2}, \cdots ,x_{n})}{\mathrm{d}(x_{1},x_{2}, \cdots ,x_{n})} \\
  &=\begin{bmatrix}
    \cfrac{\partial \mathbf{f}(\mathbf{x})}{\partial x_{1}},
    &\cfrac{\partial \mathbf{f}(\mathbf{x})}{\partial x_{2}},
    & \cdots,
    &\cfrac{\partial \mathbf{f}(\mathbf{x})}{\partial x_{n}}
  \end{bmatrix}\\
  &=\begin{bmatrix}
    \cfrac{\partial {f}_{1}(\mathbf{x})}{\partial x_{1}}& \cfrac{\partial {f}_{1}(\mathbf{x})}{\partial x_{2}} & \cdots & \cfrac{\partial {f}_{1}(\mathbf{x})}{\partial x_{n}} \\
    \cfrac{\partial {f}_{2}(\mathbf{x})}{\partial x_{1}}& \cfrac{\partial {f}_{2}(\mathbf{x})}{\partial x_{2}} & \cdots & \cfrac{\partial {f}_{2}(\mathbf{x})}{\partial x_{n}} \\
    \vdots & \vdots &  & \vdots\\
    \cfrac{\partial {f}_{m}(\mathbf{x})}{\partial x_{1}}& \cfrac{\partial {f}_{m}(\mathbf{x})}{\partial x_{2}} & \cdots & \cfrac{\partial {f}_{m}(\mathbf{x})}{\partial x_{n}} \\
  \end{bmatrix}
  \end{align}
  
-

.. math::
  \begin{align}
    \cfrac{\partial {f}_{1}(\mathbf{x})}{\partial x_{1}} & = \cfrac{\partial {f}_{1}(x_{1},x_{2}, \cdots ,x_{n})}{\partial x_{1}}\\
    \cfrac{\partial {f}_{i}(\mathbf{x})}{\partial x_{j}} & = \cfrac{\partial {f}_{i}(x_{1},x_{2}, \cdots ,x_{n})}{\partial x_{j}}
  \end{align}  