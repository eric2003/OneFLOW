Tridiagonal matrix algorithm
==================================

#. `Tridiagonal matrix algorithm <https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm/>`_
#. `Tridiagonal Matrices: Thomas Algorithm <http://www.industrial-maths.com/ms6021_thomas.pdf>`_
#. `三对角矩阵(Tridiagonal Matrices)的求法：Thomas Algorithm(TDMA) <https://www.cnblogs.com/xpvincent/archive/2013/01/25/2877411.html>`_
#. `Tridiagonal_matrix_algorithm：Thomas Algorithm(TDMA) <https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm>`_
#. `Tri-Diagonal Matrix Algorithm <http://www.thevisualroom.com/tri_diagonal_matrix.html>`_
#. `Understanding the TDMA/Thomas algorithm and its Implementation in Python <https://www.researchgate.net/publication/364389215_Understanding_the_TDMAThomas_algorithm_and_its_Implementation_in_Python/>`_
#. `Thomas Algorithm for Tridiagonal Matrix <https://matlabgeeks.weebly.com/uploads/8/0/4/8/8048228/001-thomas_algorithm_for_tridiagonal_matrix-v16.pdf>`_


In numerical linear algebra, the tridiagonal matrix algorithm, also known as the Thomas algorithm (named after Llewellyn Thomas), is a simplified form of Gaussian elimination that can be used to solve tridiagonal systems of equations.
A tridiagonal system for :math:`n` unknowns may be written as

.. math::
  a_{i}x_{i-1}+b_{i}x_{i}+c_{i}x_{i+1}=d_{i}

where :math:`a_{1}=0` and :math:`c_{n}=0`

.. math::
  \left[\begin{array}{ccccc}
  b_{1} & c_{1} & & & 0 \\
  a_{2} & b_{2} & c_{2} & & \\
  & a_{3} & b_{3} & \ddots & \\
  & & \ddots & \ddots & c_{n-1} \\
  0 & & & a_{n} & b_{n}
  \end{array}\right]\left[\begin{array}{c}
  x_{1} \\
  x_{2} \\
  x_{3} \\
  \vdots \\
  x_{n}
  \end{array}\right]=\left[\begin{array}{c}
  d_{1} \\
  d_{2} \\
  d_{3} \\
  \vdots \\
  d_{n}
  \end{array}\right]

The following is an example to deduce the solution process

.. math::
  \left[\begin{array}{cccccc}
  b_{1} & c_{1} & 0 & 0 & 0 & 0 \\
  a_{2} & b_{2} & c_{2} & 0 & 0 & 0 \\
  0 & a_{3} & b_{3} & c_{3} & 0 & 0 \\
  0 & 0 & a_{4} & b_{4} & c_{4} & 0 \\
  0 & 0 & 0 & a_{5} & b_{5} & c_{5} \\
  0 & 0 & 0 & 0 & a_{6} & b_{6}
  \end{array}\right]\left[\begin{array}{l}
  x_{1} \\
  x_{2} \\
  x_{3} \\
  x_{4} \\
  x_{5} \\
  x_{6}
  \end{array}\right]=\left[\begin{array}{l}
  d_{1} \\
  d_{2} \\
  d_{3} \\
  d_{4} \\
  d_{5} \\
  d_{6}
  \end{array}\right]
  
Step 1: Transform the matrix into an upper triangular matrix
First, the coefficient matrix in the above formula should be changed into an upper triangular matrix.  

.. math::
  b_{1}x_{1}+c_{1}x_{2}=d_{1}
  
Divide the above formula by b1:

.. math::
  x_{1}+\cfrac{c_{1}}{b_{1}}x_{2}=\cfrac{d_{1}}{b_{1}}
  
can be written as:

.. math::
  x_{1}+\gamma_{1} x_{2}=\rho_{1},\quad\gamma_{1}=\cfrac{c_{1}}{b_{1}},\quad\rho_{1}=\cfrac{d_{1}}{b_{1}} 
  
So the matrix equation can be written as:

.. math::
  \left[\begin{array}{cccccc}
  1 & \gamma_{1} & 0 & 0 & 0 & 0 \\
  a_{2} & b_{2} & c_{2} & 0 & 0 & 0 \\
  0 & a_{3} & b_{3} & c_{3} & 0 & 0 \\
  0 & 0 & a_{4} & b_{4} & c_{4} & 0 \\
  0 & 0 & 0 & a_{5} & b_{5} & c_{5} \\
  0 & 0 & 0 & 0 & a_{6} & b_{6}
  \end{array}\right]\left[\begin{array}{l}
  x_{1} \\
  x_{2} \\
  x_{3} \\
  x_{4} \\
  x_{5} \\
  x_{6}
  \end{array}\right]=\left[\begin{array}{l}
  \rho_{1} \\
  d_{2} \\
  d_{3} \\
  d_{4} \\
  d_{5} \\
  d_{6}
  \end{array}\right]  
  
second line:

.. math::
  a_{2}x_{1}+b_{2}x_{2}+c_{2}x_{3}=d_{2}  
  
Multiply the transformed first line by :math:`a_{2}`, and then subtract it from the second line to eliminate :math:`x_{1}`, and get:  

.. math::
  (b_{2}-a_{2}\gamma_{1})x_{2}+c_{2}x_{3}=d_{2}-a_{2}\rho_{1}
  
- 
 
.. math::  
  x_{2}+\cfrac{c_{2}}{(b_{2}-a_{2}\gamma_{1})}x_{3}=\cfrac{d_{2}-a_{2}\rho_{1}}{(b_{2}-a_{2}\gamma_{1})}  
  
- 
 
.. math:: 
  x_{2}+\gamma_{2}x_{3}=\rho_{2},\quad
  \gamma_{2}=\cfrac{c_{2}}{(b_{2}-a_{2}\gamma_{1})},\quad\rho_{2}=\cfrac{d_{2}-a_{2}\rho_{1}}{(b_{2}-a_{2}\gamma_{1})}\\
  
  
So the new matrix equation is:

.. math::
  \left[\begin{array}{cccccc}
  1 & \gamma_{1} & 0 & 0 & 0 & 0 \\
  0 & 1 & \gamma_{2} & 0 & 0 & 0 \\
  0 & a_{3} & b_{3} & c_{3} & 0 & 0 \\
  0 & 0 & a_{4} & b_{4} & c_{4} & 0 \\
  0 & 0 & 0 & a_{5} & b_{5} & c_{5} \\
  0 & 0 & 0 & 0 & a_{6} & b_{6}
  \end{array}\right]\left[\begin{array}{l}
  x_{1} \\
  x_{2} \\
  x_{3} \\
  x_{4} \\
  x_{5} \\
  x_{6}
  \end{array}\right]=\left[\begin{array}{l}
  \rho_{1} \\
  \rho_{2} \\
  d_{3} \\
  d_{4} \\
  d_{5} \\
  d_{6}
  \end{array}\right]
  
In the same way,  

.. math::
  x_{3}+\gamma_{3}x_{4}=\rho_{3},\quad
  \gamma_{3}=\cfrac{c_{3}}{(b_{3}-a_{3}\gamma_{2})},\quad\rho_{3}=\cfrac{d_{3}-a_{3}\rho_{2}}{(b_{3}-a_{3}\gamma_{2})}\\
  
-
  
.. math::  
  x_{4}+\gamma_{4}x_{5}=\rho_{4},\quad
  \gamma_{4}=\cfrac{c_{4}}{(b_{4}-a_{4}\gamma_{3})},\quad\rho_{4}=\cfrac{d_{4}-a_{4}\rho_{3}}{(b_{4}-a_{4}\gamma_{3})}\\
  
-
  
.. math::  
  x_{5}+\gamma_{5}x_{6}=\rho_{5},\quad
  \gamma_{5}=\cfrac{c_{5}}{(b_{5}-a_{5}\gamma_{4})},\quad\rho_{5}=\cfrac{d_{5}-a_{5}\rho_{4}}{(b_{5}-a_{5}\gamma_{4})}\\
  
-
  
.. math::
  x_{6}=\rho_{6},\quad\rho_{6}=\cfrac{d_{6}-a_{6}\rho_{5}}{(b_{6}-a_{6}\gamma_{5})}  
  
Finally, the new upper triangular matrix formula is obtained as:  

.. math::
  \left[\begin{array}{cccccc}
  1 & \gamma_{1} & 0 & 0 & 0 & 0 \\
  0 & 1 & \gamma_{2} & 0 & 0 & 0 \\
  0 & 0 & 1 & \gamma_{3} & 0 & 0 \\
  0 & 0 & 0 & 1 & \gamma_{4} & 0 \\
  0 & 0 & 0 & 0 & 1 & \gamma_{5} \\
  0 & 0 & 0 & 0 & 0 & 1
  \end{array}\right]\left[\begin{array}{l}
  x_{1} \\
  x_{2} \\
  x_{3} \\
  x_{4} \\
  x_{5} \\
  x_{6}
  \end{array}\right]=\left[\begin{array}{l}
  \rho_{1} \\
  \rho_{2} \\
  \rho_{3} \\
  \rho_{4} \\
  \rho_{5} \\
  \rho_{6}
  \end{array}\right]
  
Step 2: Solving
The reverse order of :math:`x` can be obtained as follows:  

.. math::
  \begin{array}{l}
  x_{6}=\rho_{6}\\
  x_{5}=\rho_{5}-\gamma_{5}x_{6}\\
  x_{4}=\rho_{4}-\gamma_{4}x_{5}\\
  x_{3}=\rho_{3}-\gamma_{3}x_{4}\\
  x_{2}=\rho_{2}-\gamma_{2}x_{3}\\
  x_{1}=\rho_{1}-\gamma_{1}x_{2}\\
  \end{array}
  
In general, the following relation holds  

.. math::
  \begin{array}{l}
  x_{n}=\rho_{n}\\
  x_{i}=\rho_{i}-\gamma_{i}x_{i+1};\quad i=n-1,n-2,\cdots,1
  \end{array}  
  
-
  
.. math::  
  \rho_{i}=\left\{\begin{array}{l}
  \cfrac{d_{1}}{b_{1}}\ \ \quad \quad \quad;& i=1\\
  \cfrac{d_{i}-a_{i}\rho_{i-1}}{b_{i}-a_{i}\gamma_{i-1}};&i=2,3,\cdots,n
  \end{array}\right. 
  
-
  
.. math::  
  \gamma_{i}=\left\{\begin{array}{l}
  \cfrac{c_{1}}{b_{1}}\ \ \quad \quad \quad;& i=1\\
  \cfrac{c_{i}}{b_{i}-a_{i}\gamma_{i-1}};&i=2,3,\cdots,n-1
  \end{array}\right.
  
We can rewrite this as

.. math::  
  \hat{c}_{i}=\left\{\begin{array}{l}
  \cfrac{c_{1}}{b_{1}}\ \ \quad \quad \quad;& i=1\\
  \cfrac{c_{i}}{b_{i}-a_{i}\hat{c}_{i-1}};&i=2,3,\cdots,n-1
  \end{array}\right.
  
-
  
.. math::  
  \hat{d}_{i}=\left\{\begin{array}{l}
  \cfrac{d_{1}}{b_{1}}\ \ \quad \quad \quad;& i=1\\
  \cfrac{d_{i}-a_{i}\hat{d}_{i-1}}{b_{i}-a_{i}\hat{c}_{i-1}};&i=2,3,\cdots,n
  \end{array}\right.  
  
-
  
.. math::  
  {x}_{i}=\left\{\begin{array}{l}
  \hat{d}_{n}\;\; \quad \quad \quad;& i=n\\
  \hat{d}_{i}-\hat{c}_{i}{x}_{i+1};&i=n-1,n-2,\cdots,1
  \end{array}\right.   
  
Problem 1. Solve the following Tridiagonal matrix using the Thomas algorithm:

.. math::  
  \left[\begin{array}{crccr}
  2 & -1 & 0 & 0 & 0 \\
  -1 & 2 & -1 & 0 & 0 \\
  0 & -1 & 2 & -1 & 0 \\
  0 & 0 & -1 & 2 & -1 \\
  0 & 0 & 0 & -1 & 2
  \end{array}\right]\left[\begin{array}{l}
  x_{1} \\
  x_{2} \\
  x_{3} \\
  x_{4} \\
  x_{5}
  \end{array}\right]=\left[\begin{array}{l}
  1 \\
  1 \\
  1 \\
  1 \\
  1
  \end{array}\right]
  
Results:

.. math::  
  \left[\begin{array}{l}
  x_{1} \\
  x_{2} \\
  x_{3} \\
  x_{4} \\
  x_{5}
  \end{array}\right]=\left[\begin{array}{l}
  2.5 \\
  4.0 \\
  4.5 \\
  4.0 \\
  2.5
  \end{array}\right]

  