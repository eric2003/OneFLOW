Iterative Solver 
==================================

#. `Scientific Computing || 02 Week 7 19 1 Introduction to spectral methods 10 46 <https://www.youtube.com/watch?v=ymsY8IFbOwY/>`_
#. `Solving the Discrete Poisson Equation using Jacobi, SOR, Conjugate Gradients, and the FFT <https://people.eecs.berkeley.edu/~demmel/cs267/lecture24/lecture24.html>`_
#. `Iteration methods <https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html>`_
#. `Numerical Methods for Elliptic Equations <http://www.fem.unicamp.br/~phoenics/SITE_PHOENICS/Apostilas/CFD-1_U%20Michigan_Hong/Lecture11.pdf>`_
#. `Iterative Methods: Part 2 <https://crunchingnumbers.live/2017/07/09/iterative-methods-part-2/>`_
#. `FFTW Tutorial <https://github.com/jonathanschilling/fftw_tutorial/>`_
#. `Finite difference method for 1D Poisson equation with mixed boundary conditions <https://mathematica.stackexchange.com/questions/220627/finite-difference-method-for-1d-poisson-equation-with-mixed-boundary-conditions/>`_


Iterative Solver 
--------------------------------
Iterative methods use successive approximations to obtain the most accurate solution to a linear system at every iteration step. These methods start with an initial guess and proceed to generate a sequence of approximations, in which the k-th approximation is derived from the previous ones. There are two main types of iterative methods. Stationary methods that are easy to understand and are simple to implement, but are not very effective. Non-stationary methods which are based on the idea of the sequence of orthogonal vectors. We would like to recommend a text by Barrett et al. which provides a good discussion on the implementation of iterative methods for solving a linear system of equations.

Stationary Methods: Gauss-Seidel
------------------------------------

The iterative methods work by splitting the matrix :math:`\mathbf{A}` into

.. math::
   \mathbf{A}=\mathbf{M}-\mathbf{P}

Equation 
   
.. math::
   \mathbf{A}\mathbf{u}=\mathbf{b}
   
becomes
   
.. math::
   \mathbf{M}\mathbf{u}=\mathbf{P}\mathbf{u}+\mathbf{b}
   

The solution at :math:`(k+1)` iteration is given by

.. math::
   \mathbf{M}\mathbf{u}^{(k+1)}=\mathbf{P}\mathbf{u}^{(k+1)}+\mathbf{b}
   
If we take the difference between the above two equations, we get the evolution of error as
:math:`\epsilon^{(k+1)}=\mathbf{M}^{-1}\mathbf{P}\epsilon^{(k)}`. For the solution to converge to the exact solution and the error to go to zero,
the largest eigenvalue of iteration matrix :math:`\mathbf{M}^{-1}\mathbf{P}` should be less than 1.  The approximate number of iterations required for the error
:math:`\epsilon` to go below some specified tolerance :math:`\delta` is given by

.. math::
  Q=\cfrac{\text{ln}(\delta)}{\text{ln}(\lambda_{1})} 
  
where :math:`Q` is the number of iterations and :math:`\lambda_{1}` is the maximum eigenvalue of iteration matrix :math:`\mathbf{M}^{-1}\mathbf{P}`.
If the convergence tolerance is specified to be :math:`1\times10^{-5}` then the number of iterations for convergence will be :math:`Q=1146` and 
:math:`Q=23,020` for :math:`\lambda_{1}=0.99` and :math:`\lambda_{1}=0.9995` respectively. Therefore, the maximum eigenvalue of the iteration matrix should be less for faster convergence. When we implement iterative methods, we do not form the matrix
:math:`\mathbf{M}` and :math:`\mathbf{P}`. Equation :math:`\mathbf{M}\mathbf{u}^{(k+1)}=\mathbf{P}\mathbf{u}^{(k+1)}+\mathbf{b}` is just the matrix representation of the iterative method. In matrix terms, the Gauss-Seidel method can be expressed as

.. math::
  \begin{array}{l}
  \mathbf{M}= \mathbf{D}-\mathbf{L}\\
  \mathbf{u}^{(k+1)}=(\mathbf{D}-\mathbf{L})^{-1}\mathbf{P}\mathbf{u}^{(k)}+\mathbf{M}^{-1}\mathbf{b}
  \end{array}
  
  
where :math:`\mathbf{D}`, :math:`\mathbf{L}` and :math:`\mathbf{U}` represent the diagonal, the strictly lower-triangular, and the strictly upper-triangular parts of
:math:`\mathbf{A}` respectively. We do not explicitly construct the matrices :math:`\mathbf{D}`, :math:`\mathbf{L}` and :math:`\mathbf{U}`, instead we use the update formula based on data vectors (i.e., we obtain a vector once a matrix operates on a vector).

The update formulas for Gauss-Seidel iterations if we iterate from left to right and from bottom to top can be written as

.. math::
  \begin{array}{l}
  d=\cfrac{-2}{\Delta x^2}+\cfrac{-2}{\Delta y^2}\\
  u_{i,j}^{(k+1)}=u_{i,j}^{(k)}+\cfrac{r_{i,j}}{d} \\
  r_{i,j}=f_{i,j}^{(k)}-\cfrac{u_{i-1,j}^{(k)}-2u_{i,j}^{(k)}+u_{i+1,j}^{(k)}}{\Delta x^{2}}
  -\cfrac{u_{i,j-1}^{(k)}-2u_{i,j}^{(k)}+u_{i,j+1}^{(k)}}{\Delta y^{2}}
  \end{array}
  
where we define an operator 

.. math::
  Au_{i,j}=\cfrac{u_{i-1,j}^{(k)}-2u_{i,j}^{(k)}+u_{i+1,j}^{(k)}}{\Delta x^{2}}
  +\cfrac{u_{i,j-1}^{(k)}-2u_{i,j}^{(k)}+u_{i,j+1}^{(k)}}{\Delta y^{2}}
  
Jacobi's Method
--------------------------
Iterative Methods for Solving :math:`\mathbf{A}\mathbf{x}=\mathbf{b}` - Jacobi's Method 

Two assumptions made on Jacobi Method:

1. The system given by 

.. math::
  \begin{array}{l}
  a_{11}x_{1}+a_{12}x_{2}+\cdots +a_{1n}x_{n}=b_{1}\\
  a_{21}x_{1}+a_{22}x_{2}+\cdots +a_{2n}x_{n}=b_{2}\\
  \cdots \\
  a_{n1}x_{1}+a_{n2}x_{2}+\cdots +a_{nn}x_{n}=b_{n}\\
  \end{array}
  
Has a unique solution.

2. The coefficient matrix has no zeros on its main diagonal, namely  
:math:`a_{11},a_{22},\cdots,a_{nn}` are nonzeros.

Main idea of Jacobi:

To begin, solve the :math:`1^{st}` equation for :math:`x_{1}`, the :math:`2^{nd}` equation :math:`x_{2}` for and so on to obtain the rewritten equations:

.. math::
  \begin{array}{l}
  x_{1}=\cfrac{1}{a_{11}}(b_{1}-a_{12}x_{2}-a_{13}x_{3}-\cdots-a_{1n}x_{n})\\
  x_{2}=\cfrac{1}{a_{22}}(b_{2}-a_{21}x_{1}-a_{23}x_{3}-\cdots-a_{2n}x_{n})\\
  \cdots \\
  x_{n}=\cfrac{1}{a_{nn}}(b_{n}-a_{n1}x_{1}-a_{n2}x_{2}-\cdots-a_{n,n-1}x_{n-1})\\
  \end{array}
  
Then make an initial guess of the solution :math:`x^{(0)}=(x_{1}^{(0)},x_{2}^{(0)},x_{3}^{(0)},\cdots,x_{n}^{(0)})`.
Substitute these values into the right hand side the of the rewritten equations to obtain the first approximation, 
:math:`(x_{1}^{(1)},x_{2}^{(1)},x_{3}^{(1)},\cdots,x_{n}^{(1)})`.
This accomplishes one iteration. In the same way, the second approximation :math:`(x_{1}^{(2)},x_{2}^{(2)},x_{3}^{(2)},\cdots,x_{n}^{(2)})` is computed by substituting the first approximation’s :math:`x`-
values into the right hand side of the rewritten equations. 

By repeated iterations, we form a sequence of approximations :math:`\mathbf{x}^{(k)}=(x_{1}^{(2)},x_{2}^{(2)},x_{3}^{(2)},\cdots,x_{n}^{(2)})^{T},k=1,2,3,\cdots`

For each :math:`k \ge 1`, generate the components :math:`x_{i}^{(k)}` of :math:`\mathbf{x}^{(k)}` from :math:`\mathbf{x}^{(k-1)}` by

.. math::
  x_{i}^{(k)}=\cfrac{1}{a_{ii}} \bigg[\sum_{j=1,j\ne i}^{n}(-a_{ij}x_{j}^{(k-1)})+b_{i}\bigg]
  
for 

.. math::  
  i=1,2,\cdots,n
  
The Jacobi Method in Matrix Form  

.. math::
  \mathbf{A}=\begin{bmatrix}
  a_{11}&a_{12}&\cdots&a_{1n} \\
  a_{21}&a_{22}&\cdots&a_{2n} \\
  \vdots&\vdots  & \ddots  &\vdots\\
  a_{n1}&a_{n2}&\cdots&a_{nn} \\
  \end{bmatrix}  
  
and  
  
.. math::  
  \mathbf{b}=\begin{bmatrix}
  b_{1}\\ b_{2}\\\vdots\\ b_{n}
  \end{bmatrix}  
  
for 
  
.. math::  
  \mathbf{x}=\begin{bmatrix}
  x_{1}\\ x_{2}\\\vdots\\ x_{n}
  \end{bmatrix}  
  
We split :math:`\mathbf{A}` into   

.. math::  
  \mathbf{A}=\mathbf{L}+\mathbf{D}+\mathbf{U}
  
where  

.. math::  
  \mathbf{L}=\begin{bmatrix}
  0&0&\cdots&0 \\
  a_{21}&0&\cdots&0 \\
  \vdots&\vdots  & \ddots  &\vdots\\
  a_{n1}&a_{n2}&\cdots&0\\
  \end{bmatrix}\quad
  \mathbf{D}=\begin{bmatrix}
  a_{11}&0&\cdots&0 \\
  0&a_{22}&\cdots&0 \\
  \vdots&\vdots  & \ddots  &\vdots\\
  0&0&\cdots&a_{nn}\\
  \end{bmatrix}\quad
  \mathbf{U}=\begin{bmatrix}
  0&a_{12}&\cdots&a_{1n} \\
  0&0&\cdots&a_{2n} \\
  \vdots&\vdots  & \ddots  &\vdots\\
  0&0&\cdots&0\\
  \end{bmatrix}\quad
  
:math:`\mathbf{A}\mathbf{x}=\mathbf{b}` is transformed into 

.. math:: 
  \begin{array}{l}
  (\mathbf{L}+\mathbf{D}+\mathbf{U})\mathbf{x}=\mathbf{b}\\
  {\mathbf{D}}\mathbf{x}=\mathbf{b}-(\mathbf{L}+\mathbf{U})\mathbf{x}\\
  \mathbf{x}= {\mathbf{D}^{-1}}\mathbf{b}-{\mathbf{D}^{-1}}(\mathbf{L}+\mathbf{U})\mathbf{x}\\
  \end{array}
  
where  

.. math:: 
  \mathbf{D}^{-1}=\begin{bmatrix}
  \cfrac{1}{a_{11}} &0&\cdots&0 \\
  0&\cfrac{1}{a_{22}}&\cdots&0 \\
  \vdots&\vdots  & \ddots  &\vdots\\
  0&0&\cdots&\cfrac{1}{a_{nn}}\\
  \end{bmatrix}\quad  
  
The matrix form of Jacobi iterative method is

.. math::
  \mathbf{x}^{(k+1)}= -{\mathbf{D}^{-1}}(\mathbf{L}+\mathbf{U})\mathbf{x}^{(k)}+{\mathbf{D}^{-1}}\mathbf{b}\quad k=1,2,3\\
  
- 
 
.. math::
  \begin{array}{l}
  a_{11}x_{1}^{(k+1)}+a_{12}x_{2}^{(k)}+\cdots +a_{1n}x_{n}^{(k)}=b_{1}\\
  a_{21}x_{1}^{(k)}+a_{22}x_{2}^{(k+1)}+\cdots +a_{2n}x_{n}^{(k)}=b_{2}\\
  \cdots \\
  a_{n1}x_{1}^{(k)}+a_{n2}x_{2}^{(k)}+\cdots +a_{nn}x_{n}^{(k+1)}=b_{n}\\
  \end{array}  
  
Gauss-Seidel iterative method
----------------------------------
For each :math:`k \ge 1`, generate the components :math:`x_{i}^{(k)}` of :math:`\mathbf{x}^{(k)}` from :math:`\mathbf{x}^{(k-1)}` by

.. math::
  x_{i}^{(k)}=\cfrac{1}{a_{ii}} \bigg[-\sum_{j=1}^{i-1}(a_{ij}x_{j}^{(k)})-\sum_{j=i+1}^{n}(a_{ij}x_{j}^{(k-1)})+b_{i}\bigg]
  
Namely,

.. math::
  \begin{array}{l}
  a_{11}x_{1}^{(k+1)}+a_{12}x_{2}^{(k)}+a_{13}x_{3}^{(k)}+\cdots +a_{1n}x_{n}^{(k)}=b_{1}\\
  a_{21}x_{1}^{(k+1)}+a_{22}x_{2}^{(k+1)}+a_{23}x_{3}^{(k)}+\cdots +a_{2n}x_{n}^{(k)}=b_{2}\\
  a_{31}x_{1}^{(k+1)}+a_{32}x_{2}^{(k+1)}+a_{33}x_{3}^{(k+1)}+\cdots +a_{2n}x_{n}^{(k)}=b_{3}\\
  \cdots \\
  a_{n1}x_{1}^{(k+1)}+a_{n2}x_{2}^{(k+1)}+a_{n3}x_{3}^{(k+1)}+\cdots +a_{nn}x_{n}^{(k+1)}=b_{n}\\
  \end{array}  
  
Matrix-based formula

The solution is obtained iteratively via

.. math::
  \begin{array}{l}
  (\mathbf{L}+\mathbf{D})\mathbf{x}^{(k+1)}+\mathbf{U}\mathbf{x}^{(k)}=\mathbf{b}\\
  (\mathbf{L}+\mathbf{D})\mathbf{x}^{(k+1)}=\mathbf{b}-\mathbf{U}\mathbf{x}^{(k)}\\
  \mathbf{L}_{*}=\mathbf{L}+\mathbf{D}\\
  \mathbf{L}_{*}\mathbf{x}^{(k+1)}=\mathbf{b}-\mathbf{U}\mathbf{x}^{(k)}\\
  \end{array}
  
Why the matrix-based formula works  
----------------------------------------
The system of linear equations may be rewritten as:

.. math::
  \begin{array}{l}
  \mathbf{A}\mathbf{x}=\mathbf{b}\\
  (\mathbf{L}_{*}+\mathbf{U})\mathbf{x}=\mathbf{b}\\
  \mathbf{L}_{*}\mathbf{x}+\mathbf{U}\mathbf{x}=\mathbf{b}\\
  \mathbf{L}_{*}\mathbf{x}=\mathbf{b}-\mathbf{U}\mathbf{x}\\
  \mathbf{x}=\mathbf{L}_{*}^{-1}(\mathbf{b}-\mathbf{U}\mathbf{x})\\
  \end{array}

The Gauss–Seidel method now solves the left hand side of this expression for :math:`\mathbf{x}`, using previous value for 
:math:`\mathbf{x}`, on the right hand side. Analytically, this may be written as:

.. math::
  \mathbf{x}^{(k+1)}=\mathbf{L}_{*}^{-1}(\mathbf{b}-\mathbf{U}\mathbf{x}^{(k)})\\
  
Solving the 1D Poisson equation using finite differences
--------------------------------------------------------------
Consider the 1D Poisson equation

.. math::
  \cfrac{d^{2}u}{dx^{2}} =f(x)=-1
  
on :math:`\Omega=[0,1]` with boundary conditions
  
.. math:: 
  \begin{array}{l}
  u(0)=0\\
  u'(1)=\cfrac{du}{dx}\bigg|_{x=1}=0 
  \end{array}
  
which has analytical solution

.. math:: 
  u=x-\cfrac{1}{2}x^{2}
  
Finite differences
``````````````````````````````````

We we will use this specific example to investigate various approaches to solving partial differential equations with finite differences, in which we discretize the domain by defining :math:`N` equally
spaced points  

.. math:: 
  \begin{array}{l}
  \cfrac{u_{i+1}-2u_{i}+u_{i-1}}{\Delta x^{2}} =f_{i}\\
  {u_{i+1}-2u_{i}+u_{i-1}}={\Delta x^{2}}f_{i}\\
  \end{array}
  
The boundary conditions require special care. For :math:`x=0` we have a Dirichlet boundary condition
which allows us to fix the value :math:`u_{1}=0`. For :math:`x=1` we have a Neumann boundary condition
:math:`du/dx = 0`. This is a symmetry boundary condition, so that in this case we can imagine a ’ghost’
point :math:`u_{N+1}` which is always equal to :math:`u_{N-1}`. This leads to the expression for point :math:`x_{N}`:

.. math:: 
  u_{N}=u_{N-1}
  
-  
  
.. math::   
  \begin{array}{l}
  \cfrac{u_{i+1}-2u_{i}+u_{i-1}}{\Delta x^{2}} =f_{i}\\
  {u_{i+1}-2u_{i}+u_{i-1}}={\Delta x^{2}}f_{i}=b{i}\\
  {u_{i+1}-2u_{i}+u_{i-1}}=b_{i}\\
  u_{i-1}-2u_{i}+u_{i+1}=b_{i}\\
  \end{array}  
  
-  
  
.. math:: 
  \begin{array}{l}
  u_{2-1}-2u_{2}+u_{2+1}=b_{2}\\
  u_{3-1}-2u_{3}+u_{3+1}=b_{3}\\
  \cdots\\
  u_{i-1}-2u_{i}+u_{i+1}=b_{i}\\
  \cdots\\
  u_{N-1-1}-2u_{N-1}+u_{N-1+1}=b_{N-1}\\
  \end{array}  
  
-  
  
.. math:: 
  \begin{array}{l}
  u_{1}-2u_{2}+u_{3}=b_{2}\\
  u_{2}-2u_{3}+u_{4}=b_{3}\\
  \cdots\\
  u_{i-1}-2u_{i}+u_{i+1}=b_{i}\\
  \cdots\\
  u_{N-2}-2u_{N-1}+u_{N}=b_{N-1}\\
  \end{array}   
  
-  
  
.. math::   
  \begin{array}{l}
  -2u_{2}+u_{3}=b_{2}-u_{1}=\hat{b}_{2}\\
  u_{2}-2u_{3}+u_{4}=b_{3}=\hat{b}_{3}\\
  \cdots\\
  u_{i-1}-2u_{i}+u_{i+1}=b_{i}=\hat{b}_{i}\\
  \cdots\\
  u_{N-3}-2u_{N-2}+u_{N-1}=b_{N-2}=\hat{b}_{N-2}\\
  u_{N-2}-2u_{N-1}=b_{N-1}-u_{N}=\hat{b}_{N-1}\\
  \end{array}    
  
-  
  
.. math::  
  \begin{bmatrix}
  -2&1  &0 &\cdots&\cdots&0\\
  1& -2 & 1&\ddots&&\vdots\\
  0&1& -2 & 1&\ddots&\vdots\\
  \vdots&\ddots&\ddots&\ddots&\ddots&0  \\
  \vdots&&\ddots&1& -2 & 1\\
  0&\cdots&\cdots& 0&1&-2 \\
  \end{bmatrix}
  \begin{bmatrix}
  u_{2}\\u_{3}\\\vdots\\u_{i}\\\vdots\\u_{N-2}\\u_{N-1}
  \end{bmatrix}=
  \begin{bmatrix}
  \hat{b}_{2}\\\hat{b}_{3}\\\vdots\\\hat{b}_{i}\\\vdots\\\hat{b}_{N-2}\\\hat{b}_{N-1}
  \end{bmatrix}
  
-  
  
.. math::  
  \mathbf{A}=
  \begin{bmatrix}
  -2&1  &0 &\cdots&\cdots&0\\
  1& -2 & 1&\ddots&&\vdots\\
  0&1& -2 & 1&\ddots&\vdots\\
  \vdots&\ddots&\ddots&\ddots&\ddots&0  \\
  \vdots&&\ddots&1& -2 & 1\\
  0&\cdots&\cdots& 0&1&-2 \\
  \end{bmatrix}  
  
-  
  
.. math::  
  \begin{array}{l}
  a_{ii}=-2\\
  a_{i,i+1}=1\\
  a_{i+1,i}=1\\
  \text{else}\\
  a_{i,j}=0\\
  \end{array}    
  
-  
  
.. math:: 
  \begin{array}{l}
  -2u_{2}^{(k+1)}+u_{3}^{(k)}=\hat{b}_{2}\\
  u_{2}^{(k+1)}-2u_{3}^{(k+1)}+u_{4}^{(k)}=\hat{b}_{3}\\
  u_{3}^{(k+1)}-2u_{4}^{(k+1)}+u_{5}^{(k)}=\hat{b}_{4}\\
  \cdots\\
  u_{i-1}^{(k+1)}-2u_{i}^{(k+1)}+u_{i+1}^{(k)}=\hat{b}_{i}\\
  \cdots\\
  u_{N-3}^{(k+1)}-2u_{N-2}^{(k+1)}+u_{N-1}^{(k)}=\hat{b}_{N-2}\\
  u_{N-2}^{(k+1)}-2u_{N-1}^{(k+1)}=\hat{b}_{N-1}\\
  \end{array}    
  
  
-  
  
.. math:: 
  \begin{array}{l}
  \displaystyle \cfrac{u_{i-1}-2u_{i}+u_{i+1}}{\Delta x^{2}} =f_{i}\\
  \displaystyle \cfrac{u_{i-1}^{(k+1)}-2u_{i}^{(k+1)}+u_{i+1}^{(k)}}{\Delta x^{2}} =f_{i}\\
  \displaystyle \cfrac{u_{i-1}^{(k+1)}-2u_{i}^{(k)}+2u_{i}^{(k)}-2u_{i}^{(k+1)}+u_{i+1}^{(k)}}{\Delta x^{2}} =f_{i}\\
  \end{array}  
  
-  
  
.. math::
  \begin{array}{l}
  \displaystyle \cfrac{u_{i-1}^{(k+1)}-2u_{i}^{(k)}-2(u_{i}^{(k+1)}-u_{i}^{(k)})+u_{i+1}^{(k)}}{\Delta x^{2}} =f_{i}\\
  \displaystyle \cfrac{u_{i-1}^{(k+1)}-2u_{i}^{(k)}+u_{i+1}^{(k)}}{\Delta x^{2}}
  -\cfrac{2}{\Delta x^{2}}(u_{i}^{(k+1)}-u_{i}^{(k)})=f_{i}\\
  -\cfrac{2}{\Delta x^{2}}(u_{i}^{(k+1)}-u_{i}^{(k)})=f_{i}-\displaystyle \cfrac{u_{i-1}^{(k+1)}-2u_{i}^{(k)}+u_{i+1}^{(k)}}{\Delta x^{2}}\\
  \end{array}  
  
-  
  
.. math::
  \begin{array}{l}
  u_{i}^{(k+1)}-u_{i}^{(k)}=\cfrac{1}{-\cfrac{2}{\Delta x^{2}}} \Bigg\{f_{i}-\displaystyle \cfrac{u_{i-1}^{(k+1)}-2u_{i}^{(k)}+u_{i+1}^{(k)}}{\Delta x^{2}}\Bigg\}\\
  u_{i}^{(k+1)}=u_{i}^{(k)}+\cfrac{1}{-\cfrac{2}{\Delta x^{2}}} \Bigg\{f_{i}-\displaystyle \cfrac{u_{i-1}^{(k+1)}-2u_{i}^{(k)}+u_{i+1}^{(k)}}{\Delta x^{2}}\Bigg\}\\
  \end{array}  
  
Let  

.. math::
  r_{i}=f_{i}-\displaystyle \cfrac{u_{i-1}^{(k+1)}-2u_{i}^{(k)}+u_{i+1}^{(k)}}{\Delta x^{2}}

then

.. math::
  \begin{array}{l}
  d_{i}= -\cfrac{2}{\Delta x^{2}}\\
  u_{i}^{(k+1)}=u_{i}^{(k)}+\cfrac{1}{-\cfrac{2}{\Delta x^{2}}}r_{i}\\
  u_{i}^{(k+1)}=u_{i}^{(k)}+\cfrac{1}{d_{i}}r_{i}\\
  \end{array}
  
Solving the 2D Poisson equation using finite differences  
--------------------------------------------------------------

.. math::
  \displaystyle \cfrac{u_{i-1,j}-2u_{i,j}+u_{i+1,j}}{\Delta x^{2}} 
  +\cfrac{u_{i,j-1}-2u_{i,j}+u_{i,j+1}}{\Delta y^{2}} 
  =f_{i,j}\\
  
-
  
.. math::
  \begin{array}{l}
  \displaystyle \cfrac{u_{2-1,j}-2u_{2,j}+u_{2+1,j}}{\Delta x^{2}} 
  +\cfrac{u_{2,j-1}-2u_{2,j}+u_{2,j+1}}{\Delta y^{2}} =f_{2,j}\\
  \displaystyle \cfrac{u_{3-1,j}-2u_{3,j}+u_{3+1,j}}{\Delta x^{2}} 
  +\cfrac{u_{3,j-1}-2u_{3,j}+u_{3,j+1}}{\Delta y^{2}} =f_{3,j}\\
  \cdots \\
  \displaystyle \cfrac{u_{i-1,j}-2u_{i,j}+u_{i+1,j}}{\Delta x^{2}} 
  +\cfrac{u_{i,j-1}-2u_{i,j}+u_{i,j+1}}{\Delta y^{2}} =f_{i,j}\\
  \cdots \\
  \displaystyle \cfrac{u_{m-1-1,j}-2u_{m-1,j}+u_{m-1+1,j}}{\Delta x^{2}} 
  +\cfrac{u_{m-1,j-1}-2u_{m-1,j}+u_{m-1,j+1}}{\Delta y^{2}} =f_{m-1,j}\\
  \end{array}
  
-
  
.. math::
  \begin{array}{l}
  \displaystyle \cfrac{u_{1,j}-2u_{2,j}+u_{3,j}}{\Delta x^{2}} 
  +\cfrac{u_{2,j-1}-2u_{2,j}+u_{2,j+1}}{\Delta y^{2}} =f_{2,j}\\
  \displaystyle \cfrac{u_{2,j}-2u_{3,j}+u_{4,j}}{\Delta x^{2}} 
  +\cfrac{u_{3,j-1}-2u_{3,j}+u_{3,j+1}}{\Delta y^{2}} =f_{3,j}\\
  \cdots \\
  \displaystyle \cfrac{u_{i-1,j}-2u_{i,j}+u_{i+1,j}}{\Delta x^{2}} 
  +\cfrac{u_{i,j-1}-2u_{i,j}+u_{i,j+1}}{\Delta y^{2}} =f_{i,j}\\
  \cdots \\
  \displaystyle \cfrac{u_{m-2,j}-2u_{m-1,j}+u_{m,j}}{\Delta x^{2}} 
  +\cfrac{u_{m-1,j-1}-2u_{m-1,j}+u_{m-1,j+1}}{\Delta y^{2}} =f_{m-1,j}\\
  \end{array}
  
-
  
.. math::  
  \begin{array}{l}
  \displaystyle \cfrac{u_{1,2}-2u_{2,2}+u_{3,2}}{\Delta x^{2}} 
  +\cfrac{u_{2,1}-2u_{2,2}+u_{2,3}}{\Delta y^{2}} =f_{2,2}\\
  \displaystyle \cfrac{u_{2,2}-2u_{3,2}+u_{4,2}}{\Delta x^{2}} 
  +\cfrac{u_{3,1}-2u_{3,2}+u_{3,3}}{\Delta y^{2}} =f_{3,2}\\
  \cdots \\
  \displaystyle \cfrac{u_{i-1,2}-2u_{i,2}+u_{i+1,2}}{\Delta x^{2}} 
  +\cfrac{u_{i,1}-2u_{i,2}+u_{i,3}}{\Delta y^{2}} =f_{i,2}\\
  \cdots \\
  \displaystyle \cfrac{u_{m-2,2}-2u_{m-1,2}+u_{m,2}}{\Delta x^{2}} 
  +\cfrac{u_{m-1,1}-2u_{m-1,2}+u_{m-1,3}}{\Delta y^{2}} =f_{m-1,2}\\
  \end{array}  
  
-
  
.. math:: 
  \begin{array}{l}
  \displaystyle \cfrac{u_{1,3}-2u_{2,3}+u_{3,3}}{\Delta x^{2}} 
  +\cfrac{u_{2,2}-2u_{2,3}+u_{2,4}}{\Delta y^{2}} =f_{2,3}\\
  \displaystyle \cfrac{u_{2,3}-2u_{3,3}+u_{4,3}}{\Delta x^{2}} 
  +\cfrac{u_{3,2}-2u_{3,3}+u_{3,4}}{\Delta y^{2}} =f_{3,3}\\
  \cdots \\
  \displaystyle \cfrac{u_{i-1,3}-2u_{i,3}+u_{i+1,3}}{\Delta x^{2}} 
  +\cfrac{u_{i,2}-2u_{i,3}+u_{i,4}}{\Delta y^{2}} =f_{i,3}\\
  \cdots \\
  \displaystyle \cfrac{u_{m-2,3}-2u_{m-1,3}+u_{m,3}}{\Delta x^{2}} 
  +\cfrac{u_{m-1,2}-2u_{m-1,3}+u_{m-1,4}}{\Delta y^{2}} =f_{m-1,3}\\
  \end{array}
  
-
  
.. math::   
  \begin{array}{l}
  \displaystyle \cfrac{u_{1,n-1}-2u_{2,n-1}+u_{3,n-1}}{\Delta x^{2}} 
  +\cfrac{u_{2,n-2}-2u_{2,n-1}+u_{2,n}}{\Delta y^{2}} =f_{2,n-1}\\
  \displaystyle \cfrac{u_{2,n-1}-2u_{3,n-1}+u_{4,n-1}}{\Delta x^{2}} 
  +\cfrac{u_{3,n-2}-2u_{3,n-1}+u_{3,n}}{\Delta y^{2}} =f_{3,n-1}\\
  \cdots \\
  \displaystyle \cfrac{u_{i-1,n-1}-2u_{i,n-1}+u_{i+1,n-1}}{\Delta x^{2}} 
  +\cfrac{u_{i,n-2}-2u_{i,n-1}+u_{i,n}}{\Delta y^{2}} =f_{i,n-1}\\
  \cdots \\
  \displaystyle \cfrac{u_{m-2,n-1}-2u_{m-1,n-1}+u_{m,n-1}}{\Delta x^{2}} 
  +\cfrac{u_{m-1,n-1}-2u_{m-1,n-1}+u_{m-1,n}}{\Delta y^{2}} =f_{m-1,n-1}\\
  \end{array}
  
-
  
.. math::  
  \mathbf{u}=\begin{bmatrix}
  u_{2,2}\\u_{3,2}\\\vdots\\u_{m-1,2}\\
  u_{2,3}\\u_{3,3}\\\vdots\\u_{m-1,3}\\
  \vdots\\
  u_{2,n-1}\\u_{3,n-1}\\\vdots\\u_{m-1,n-1}\\
  \end{bmatrix}\quad 
  \mathbf{f}=\begin{bmatrix}
  f_{2,2}\\f_{3,2}\\\vdots\\f_{m-1,2}\\
  f_{2,3}\\f_{3,3}\\\vdots\\f_{m-1,3}\\
  \vdots\\
  f_{2,n-1}\\f_{3,n-1}\\\vdots\\f_{m-1,n-1}\\
  \end{bmatrix}\quad   

-

.. math::
  \displaystyle \cfrac{u_{i-1,j}-2u_{i,j}+u_{i+1,j}}{\Delta x^{2}} 
  +\cfrac{u_{i,j-1}-2u_{i,j}+u_{i,j+1}}{\Delta y^{2}} 
  =f_{i,j}\\  
  
-

.. math::
   \cfrac{1}{\Delta x^{2}}u_{i-1,j}
  -(\cfrac{2}{\Delta x^{2}}+\cfrac{2}{\Delta y^{2}})u_{i,j}
  +\cfrac{1}{\Delta x^{2}}u_{i+1,j}
  +\cfrac{1}{\Delta y^{2}}u_{i,j-1}
  +\cfrac{1}{\Delta y^{2}}u_{i,j+1}
  =f_{i,j}\\
  
-

.. math::
  \begin{array}{l}
  a=\cfrac{1}{\Delta x^{2}},
  b=-(\cfrac{2}{\Delta x^{2}}+\cfrac{2}{\Delta y^{2}}),
  c=\cfrac{1}{\Delta x^{2}},
  d=\cfrac{1}{\Delta y^{2}},
  e=\cfrac{1}{\Delta y^{2}}\\
  au_{i-1,j}+bu_{i,j}+cu_{i+1,j}+du_{i,j-1}+eu_{i,j+1}  =f_{i,j}\\
  \end{array}  
  
-

.. math::
  \begin{array}{l}
  bu_{2,2}+cu_{3,2}+eu_{2,3}=f_{2,2}-au_{1,2}-du_{2,1}=\hat{f}_{2,2}\\
  au_{2,2}+bu_{3,2}+cu_{4,2}+eu_{3,3}=f_{3,2}-du_{3,1}=\hat{f}_{3,2}\\
  au_{3,2}+bu_{4,2}+cu_{5,2}+eu_{4,3}=f_{4,2}-du_{4,1}=\hat{f}_{4,2}\\
  \cdots \\
  au_{i-1,2}+bu_{i,2}+cu_{i+1,2}+eu_{i,3}=f_{i,2}-du_{i,1}=\hat{f}_{i,2}\\
  \cdots \\
  au_{m-2,2}+bu_{m-1,2}+cu_{m,2}+eu_{m-1,3}=f_{m-1,2}-du_{m-1,1}=\hat{f}_{m-1,2}\\
  \end{array}
  
-

.. math::
  \begin{array}{l}
  au_{i-1,j}+bu_{i,j}+cu_{i+1,j}+du_{i,j-1}+eu_{i,j+1}=f_{i,j}\\
  au_{i-1,3}+bu_{i,3}+cu_{i+1,3}+du_{i,2}+eu_{i,4}=f_{i,3}\\
  bu_{2,3}+cu_{3,3}+du_{2,2}+eu_{2,4}=f_{2,3}-au_{1,3}=\hat{f}_{2,3}\\
  du_{2,2}+bu_{2,3}+cu_{3,3}+eu_{2,4}=f_{2,3}-au_{1,3}=\hat{f}_{2,3}\\
  au_{2,3}+bu_{3,3}+cu_{4,3}+du_{3,2}+eu_{3,4}=f_{3,3}=\hat{f}_{3,3}\\
  du_{3,2}+au_{2,3}+bu_{3,3}+cu_{4,3}+eu_{3,4}=f_{3,3}=\hat{f}_{3,3}\\
  au_{3,3}+bu_{4,3}+cu_{5,3}+du_{4,2}+eu_{4,4}=f_{4,3}=\hat{f}_{4,3}\\
  \cdots\\
  au_{m-2,3}+bu_{m-1,3}+cu_{m,3}+du_{m-1,2}+eu_{m-1,4}=f_{m-1,3}=\hat{f}_{m-1,3}\\
  du_{m-1,2}+au_{m-2,3}+bu_{m-1,3}+cu_{m,3}+eu_{m-1,4}=f_{m-1,3}=\hat{f}_{m-1,3}\\
  \end{array}    
  
-

.. math::
  \begin{bmatrix}
  b& c &0&\cdots &e &\cdots&\cdots&&0\\
  a& b &c&\ddots&\ &\ddots&&&\vdots\\
  0&a& b &c&\ddots &&\ddots&&\vdots\\
  \vdots&\ddots &\ddots&\ddots&&\ddots&&&e\\
  d& &\ddots&\ddots&&\ddots&&&\vdots\\
  \vdots&\ddots   & & &\ddots&a& b &c&0\\
  \vdots&  &\ddots&&&\ddots&a& b &c\\
  0&\cdots  &&d&\cdots&&0&a& b\\
  \end{bmatrix}  
  
-

.. math::
  au_{i-1,j}^{(k+1)}+bu_{i,j}^{(k+1)}+cu_{i+1,j}^{(k)}+du_{i,j-1}^{(k+1)}+eu_{i,j+1}^{(k)}  =f_{i,j}\\  
  
-

.. math::
  \displaystyle \cfrac{u_{i-1,j}^{(k+1)}-2u_{i,j}^{(k+1)}+u_{i+1,j}^{(k)}}{\Delta x^{2}} 
  +\cfrac{u_{i,j-1}^{(k+1)}-2u_{i,j}^{(k+1)}+u_{i,j+1}^{(k)}}{\Delta y^{2}} 
  =f_{i,j}\\  
  
-

.. math::
  \begin{array}{l}
  \displaystyle \cfrac{u_{i-1,j}^{(k+1)}-2u_{i,j}^{(k+1)}+u_{i+1,j}^{(k)}}{\Delta x^{2}} 
  +\cfrac{u_{i,j-1}^{(k+1)}-2u_{i,j}^{(k+1)}+u_{i,j+1}^{(k)}}{\Delta y^{2}} 
  =f_{i,j}\\
  \displaystyle \cfrac{u_{i-1,j}^{(k+1)}-2u_{i,j}^{(k)}-2(u_{i,j}^{(k+1)}-u_{i,j}^{(k)})+u_{i+1,j}^{(k)}}{\Delta x^{2}} 
  +\cfrac{u_{i,j-1}^{(k+1)}-2u_{i,j}^{(k)}-2(u_{i,j}^{(k+1)}-u_{i,j}^{(k)})+u_{i,j+1}^{(k)}}{\Delta y^{2}} 
  =f_{i,j}\\
  \displaystyle \cfrac{u_{i-1,j}^{(k+1)}-2u_{i,j}^{(k)}+u_{i+1,j}^{(k)}}{\Delta x^{2}}+\cfrac{-2}{\Delta x^{2}}(u_{i,j}^{(k+1)}-u_{i,j}^{(k)})\\ 
  +\cfrac{u_{i,j-1}^{(k+1)}-2u_{i,j}^{(k)}+u_{i,j+1}^{(k)}}{\Delta y^{2}} +\cfrac{-2}{\Delta y^{2}}(u_{i,j}^{(k+1)}-u_{i,j}^{(k)}) 
  =f_{i,j}\\
  \displaystyle \cfrac{u_{i-1,j}^{(k+1)}-2u_{i,j}^{(k)}+u_{i+1,j}^{(k)}}{\Delta x^{2}}
  +\cfrac{u_{i,j-1}^{(k+1)}-2u_{i,j}^{(k)}+u_{i,j+1}^{(k)}}{\Delta y^{2}}
  +\bigg\{\cfrac{-2}{\Delta x^{2}}+\cfrac{-2}{\Delta y^{2}}\bigg\}(u_{i,j}^{(k+1)}-u_{i,j}^{(k)}) 
  =f_{i,j}\\
  \end{array} 
  
-

.. math::
  \displaystyle \bigg\{\cfrac{-2}{\Delta x^{2}}+\cfrac{-2}{\Delta y^{2}}\bigg\}(u_{i,j}^{(k+1)}-u_{i,j}^{(k)}) =f_{i,j}-\cfrac{u_{i-1,j}^{(k+1)}-2u_{i,j}^{(k)}+u_{i+1,j}^{(k)}}{\Delta x^{2}}
  -\cfrac{u_{i,j-1}^{(k+1)}-2u_{i,j}^{(k)}+u_{i,j+1}^{(k)}}{\Delta y^{2}}\\  
  
-

.. math::
  \begin{array}{l}
  u_{i,j}^{(k+1)}-u_{i,j}^{(k)}
   =\cfrac{f_{i,j}-\cfrac{u_{i-1,j}^{(k+1)}-2u_{i,j}^{(k)}+u_{i+1,j}^{(k)}}{\Delta x^{2}}
  -\cfrac{u_{i,j-1}^{(k+1)}-2u_{i,j}^{(k)}+u_{i,j+1}^{(k)}}{\Delta y^{2}}}{\cfrac{-2}{\Delta x^{2}}+\cfrac{-2}{\Delta y^{2}}} \\
  u_{i,j}^{(k+1)}
   =u_{i,j}^{(k)}+\cfrac{f_{i,j}-\cfrac{u_{i-1,j}^{(k+1)}-2u_{i,j}^{(k)}+u_{i+1,j}^{(k)}}{\Delta x^{2}}
  -\cfrac{u_{i,j-1}^{(k+1)}-2u_{i,j}^{(k)}+u_{i,j+1}^{(k)}}{\Delta y^{2}}}{\cfrac{-2}{\Delta x^{2}}+\cfrac{-2}{\Delta y^{2}}} \\
  \end{array}  
  
Let

.. math::
  \begin{array}{l}
  r_{i,j}=f_{i,j}-\Bigg[\cfrac{u_{i-1,j}^{(k+1)}-2u_{i,j}^{(k)}+u_{i+1,j}^{(k)}}{\Delta x^{2}}
  +\cfrac{u_{i,j-1}^{(k+1)}-2u_{i,j}^{(k)}+u_{i,j+1}^{(k)}}{\Delta y^{2}}\Bigg]\\
  d_{i,j}=\cfrac{-2}{\Delta x^{2}}+\cfrac{-2}{\Delta y^{2}}
  \end{array}
  
then

.. math::
  u_{i,j}^{(k+1)}=u_{i,j}^{(k)}+\cfrac{r_{i,j}}{d_{i,j}}
  
  
Non-Stationary Methods: Conjugate Gradient Algorithm
------------------------------------------------------------
#. `共轭梯度法-苏州大学 <https://www.bilibili.com/video/BV16a4y1t76z/>`_
#. `数值分析2020春苏州大学 <https://www.bilibili.com/video/BV18741177td/>`_
#. `Linear Conjugate Gradient Algorithm <https://indrag49.github.io/Numerical-Optimization/conjugate-gradient-methods-1.html>`_
#. `Gradient Descent <https://ml-cheatsheet.readthedocs.io/en/latest/gradient_descent.html>`_
#. `如何通俗地解释梯度下降法 <https://www.bilibili.com/video/BV1a94y1S7PP/>`_
#. `One-Dimensional Gradient Descent <https://d2l.ai/chapter_optimization/gd.html>`_
#. `Easiest Way to Understand Gradient Descent Step by Step <https://www.youtube.com/watch?v=wRCczEGBi1g>`_
#. `The Concept of Conjugate Gradient Descent in Python <https://ikuz.eu/machine-learning-and-computer-science/the-concept-of-conjugate-gradient-descent-in-python/>`_
#. `Complete Step-by-step Conjugate Gradient Algorithm from Scratch <https://towardsdatascience.com/complete-step-by-step-conjugate-gradient-algorithm-from-scratch-202c07fb52a8/>`_
#. `Conjugate Gradient Descent <https://gregorygundersen.com/blog/2022/03/20/conjugate-gradient-descent/>`_
#. `Descent method — Steepest descent and conjugate gradient in Python <https://sophiamyang.github.io/DS/optimization/descentmethod2/descentmethod2.html>`_
#. `Finite difference method for 1D Poisson equation with mixed boundary conditions <https://mathematica.stackexchange.com/questions/220627/finite-difference-method-for-1d-poisson-equation-with-mixed-boundary-conditions/>`_

Non-stationary methods differ from stationary methods in that the iterative matrix changes at every iteration. These methods work by forming a basis of a sequence of matrix powers times the initial residual.
The basis is called as the Krylov subspace and mathematically given by :math:`\mathcal{K}_{n}(\mathbf{A},\mathbf{b})=\text{span}\{\mathbf{b},\mathbf{A}\mathbf{b},\mathbf{A}^{2}\mathbf{b},\cdots,\mathbf{A}^{n-1}\mathbf{b}\}`.
The approximate solution to the linear system is found by minimizing the residual over the subspace formed. In this paper, we discuss the conjugate gradient method which is one of the most effective methods for symmetric positive definite systems.

The conjugate gradient method proceeds by calculating the vector sequence of successive approximate solution, residual corresponding the approximate solution, and search direction used in updating the solution and residuals. The approximate solution
:math:`\mathbf{u}^{(k)}` is updated at every iteration by a scalar multiple :math:`\alpha_{k}` of the search direction vector
:math:`\mathbf{p}^{(k)}`:

.. math::
  \mathbf{u}^{(k+1)}=\mathbf{u}^{(k)}+\alpha_{k}\mathbf{p}^{(k)}
  
Correspondingly, the residuals :math:`\mathbf{r}^{(k)}` are updated as

.. math::
  \mathbf{r}^{(k+1)}=\mathbf{r}^{(k)}+\alpha_{k}\mathbf{q}^{(k)}  
  
where 

.. math::
  \mathbf{q}^{(k)}=\mathbf{A}\mathbf{p}^{(k)}
  
The search directions are then updated using the residuals

.. math::
  \mathbf{p}^{(k+1)}=\mathbf{r}^{(k+1)}+\beta{k}\mathbf{p}^{(k)}  
  
  
where the choice :math:`\beta_{k}=\cfrac{{\mathbf{r}^{(k)}}^{\mathbf{T}}\mathbf{r}^{(k)}}{{\mathbf{r}^{(k-1)}}^{\mathbf{T}}\mathbf{r}^{(k-1)}}` 
ensures that :math:`\mathbf{p}^{(k+1)}` and :math:`\mathbf{r}^{(k+1)}` are orthogonal to all previous :math:`\mathbf{A}\mathbf{p}^{(k)}`
and :math:`\mathbf{r}^{(k)}` respectively. 

Conjugate gradient algorithm
----------------------------------------
1. Given :math:`\mathbf{b}`
2. Given matrix operator :math:`\mathbf{A}`
3. :math:`\mathbf{u}^{(0)}=\mathbf{b}`
4. :math:`\mathbf{r}^{(0)}=\mathbf{b}-\mathbf{A}\mathbf{u}^{(0)}`
5. :math:`\mathbf{p}^{(0)}=\mathbf{r}^{(0)}`
6. :math:`k=0`
7. :math:`\rho_{0}={\mathbf{r}^{(0)}}^{\mathbf{T}}\mathbf{r}^{(0)}`
8. while tolerance met ( or :math:`k<N` ) do
9.    :math:`\mathbf{q}^{(k)}=\mathbf{A}\mathbf{p}^{(k)}`
10.   :math:`\alpha_{k}=\cfrac{\rho_{k}}{{\mathbf{p}^{(k)}}^{\mathbf{T}}\mathbf{q}^{(k)}}`
11.   :math:`\mathbf{u}^{(k+1)}=\mathbf{u}^{(k)}+\alpha_{k}\mathbf{p}^{(k)}`
12.   :math:`\mathbf{r}^{(k+1)}=\mathbf{r}^{(k)}-\alpha_{k}\mathbf{q}^{(k)}`
13.   :math:`\rho_{k+1}={\mathbf{r}^{(k+1)}}^{\mathbf{T}}\mathbf{r}^{(k+1)}`
14.   :math:`\beta_{k}=\cfrac{\rho_{k+1}}{\rho_{k}}`
15.   :math:`\mathbf{p}^{(k+1)}=\mathbf{r}^{(k+1)}+\beta_{k}\mathbf{p}^{(k)}`
16.   check convergence; continue if necessary
17.   :math:`k\gets k+1`
18. end while 

Linear Conjugate Gradient Algorithm
----------------------------------------

  
Steepest Descent Algorithm  

1. Given any initial value :math:`\mathbf{x}_{0}` calculate the residual :math:`\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}`.
2. Go along the selected direction :math:`\mathbf{p}=\mathbf{r}_{0}`, calculate

.. math::
  \alpha =\cfrac{(\mathbf{r}_{0},\mathbf{r}_{0})}{(\mathbf{A}\mathbf{r}_{0},\mathbf{r}_{0})} 
  
update

.. math::
  \mathbf{x}_{1}=\mathbf{x}_{0}+\alpha \mathbf{r}_{0}  

3. Repeat  

The Method of Steepest Descent
---------------------------------------  
Suppose we want to find the minimizer of an objective function, having the quadratic form:

.. math::
  f(\mathbf{x})=\cfrac{1}{2} \mathbf{x}^{\text{T}}\mathbf{A}\mathbf{x}-\mathbf{x}^{\text{T}} \mathbf{b}
  
where :math:`\mathbf{A}` is a :math:`n\times n` symmetric positive definite matrix. 
Since :math:`\mathbf{A}` is positive definite, the minimum value point of the quadratic function :math:`f(\mathbf{x})` must exist and be a stagnation point:

.. math::
  f(\mathbf{x}^{*})=\text{min}f(\mathbf{x})\quad\Leftrightarrow \quad
  \nabla f(\mathbf{x}^{*})=\mathbf{A}\mathbf{x}^{*}-\mathbf{b}=0
  
Calculate the gradient first

.. math::
  \begin{array}{l}
  \displaystyle f(\mathbf{x})=\cfrac{1}{2}\sum_{j=1}^{n}\sum_{k=1}^{n}a_{kj}x_{k}x_{j}-\sum_{j=1}^{n}x_{j}b_{j}\\
  \displaystyle \cfrac{\partial f}{\partial x_{l}}=\cfrac{1}{2}\sum_{k=1}^{n}a_{kl}x_{k}
  +\cfrac{1}{2}\sum_{j=1}^{n}a_{lj}x_{j}-b_{l}=\sum_{k=1}^{n}a_{kl}x_{k}-b_{l}\\
  \nabla f=\text{grad } f=\mathbf{A}\mathbf{x}-\mathbf{b}
  \end{array}
  
1. Given any initial value :math:`\mathbf{x}_{0}` calculate the residual :math:`\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}`.
2. Go along the selected direction :math:`\mathbf{p}`, for example, :math:`\mathbf{p}=\mathbf{r}_{0}` direction

.. math::
  \mathbf{x}_{1}= \mathbf{x}_{0}+\alpha \ast \mathbf{p}

-

.. math::
  \begin{array}{l}
  \nabla f=\text{grad } f=\mathbf{A}\mathbf{x}-\mathbf{b}\\
  \mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}\\
  \nabla f\bigg|_{x=x_{1}}\cdot\mathbf{p}=(\mathbf{r}_{1},\mathbf{p})=(\mathbf{b}-\mathbf{A}\mathbf{x}_{1},\mathbf{p})=0\\
  (\mathbf{b}-\mathbf{A}\mathbf{x}_{1},\mathbf{p})=
  (\mathbf{b}-\mathbf{A}\mathbf{x}_{0}-\alpha\mathbf{A}\mathbf{p},\mathbf{p})=
  (\mathbf{r}_{0}-\alpha\mathbf{A}\mathbf{p},\mathbf{p})=0\\
  \Rightarrow \alpha =\cfrac{(\mathbf{r}_{0},\mathbf{p})}{(\mathbf{A}\mathbf{p},\mathbf{p})} 
  \end{array}  
  
Gradient

.. math::
  \begin{array}{l}
  \nabla f=\text{grad } f=\mathbf{A}\mathbf{x}-\mathbf{b}\\
  -\nabla f=-\text{grad } f=\mathbf{b}-\mathbf{A}\mathbf{x}\\
  \end{array}

-
  
.. math::
  \begin{array}{l}
  f(\mathbf{x})=\cfrac{1}{2} \mathbf{x}^{\text{T}}\mathbf{A}\mathbf{x}-\mathbf{x}^{\text{T}} \mathbf{b}\\
  f(\mathbf{x}_{1})=\cfrac{1}{2} \mathbf{x}_{1}^{\text{T}}\mathbf{A}\mathbf{x}_{1}-\mathbf{x}_{1}^{\text{T}} \mathbf{b}\\
  \mathbf{x}_{1}= \mathbf{x}_{0}+\alpha \ast \mathbf{p}\\
  f(\mathbf{x}_{1})=\hat{f}(\alpha)=\cfrac{1}{2} {(\mathbf{x}_{0}+\alpha \ast \mathbf{p})}^{\text{T}}\mathbf{A}{(\mathbf{x}_{0}+\alpha \ast \mathbf{p})}-{(\mathbf{x}_{0}+\alpha \ast \mathbf{p})}^{\text{T}} \mathbf{b}\\
  \end{array}  
  
-
  
.. math::
  \begin{array}{l}
  f(\mathbf{x}_{1})=\hat{f}(\alpha)=\cfrac{1}{2} {(\mathbf{x}_{0}+\alpha \ast \mathbf{p})}^{\text{T}}\mathbf{A}{(\mathbf{x}_{0}+\alpha \ast \mathbf{p})}-{(\mathbf{x}_{0}+\alpha \ast \mathbf{p})}^{\text{T}} \mathbf{b}\\
  \cfrac{d\hat{f}}{d\alpha} =\cfrac{1}{2}\mathbf{p}^{\text{T}}\mathbf{A}{(\mathbf{x}_{0}+\alpha \ast \mathbf{p})}
  +\cfrac{1}{2} {(\mathbf{x}_{0}+\alpha \ast \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}-\mathbf{p}^{\text{T}}\mathbf{b}\\
  \cfrac{d\hat{f}}{d\alpha} =\cfrac{1}{2}\mathbf{p}^{\text{T}}\mathbf{A}\mathbf{x}_{0}
  +\cfrac{1}{2}\alpha \mathbf{p}^{\text{T}}\mathbf{A}{(\mathbf{p})}
  +\cfrac{1}{2} {(\mathbf{x}_{0})}^{\text{T}}\mathbf{A}\mathbf{p}
  +\cfrac{1}{2}\alpha  {( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}
  -\mathbf{p}^{\text{T}}\mathbf{b}\\
  \cfrac{d\hat{f}}{d\alpha} =\cfrac{1}{2}\mathbf{p}^{\text{T}}\mathbf{A}\mathbf{x}_{0}
  +\cfrac{1}{2} {(\mathbf{x}_{0})}^{\text{T}}\mathbf{A}\mathbf{p}
  +\alpha  {( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}
  -\mathbf{p}^{\text{T}}\mathbf{b}\\
  \end{array}

-
  
.. math::  
  {\displaystyle \left(\mathbf {AB} \right)^{\operatorname {T} }=\mathbf {B} ^{\operatorname {T} }\mathbf {A} ^{\operatorname {T} }.}  
  
-
  
.. math:: 
  \begin{array}{l}
  {\displaystyle \left(\mathbf {AB} \right)^{\operatorname {T} }=\mathbf {B} ^{\operatorname {T} }\mathbf {A} ^{\operatorname {T} }}  \\
  (\mathbf{p}^{\text{T}}(\mathbf{A}\mathbf{x}_{0}))^{\text{T}}=(\mathbf{A}\mathbf{x}_{0})^{\text{T}}(\mathbf{p}^{\text{T}})^{\text{T}}
  =(\mathbf{A}\mathbf{x}_{0})^{\text{T}}\mathbf{p}=\mathbf{x}_{0}^{\text{T}}\mathbf{A}^{\text{T}}\mathbf{p}
  \end{array} 

where :math:`\mathbf{A}` is a :math:`n\times n` symmetric positive definite matrix. 
then

.. math:: 
  \mathbf {A}=\mathbf {A}^{\text{T}}\\

-

.. math:: 
  \mathbf{x}_{0}^{\text{T}}\mathbf{A}^{\text{T}}\mathbf{p}=\mathbf{x}_{0}^{\text{T}}\mathbf{A}\mathbf{p}

  
Here again :math:`\mathbf{p}^{\text{T}}\mathbf{A}\mathbf{x}_{0}`, 
:math:`{(\mathbf{x}_{0})}^{\text{T}}\mathbf{A}\mathbf{p}`  and :math:`\mathbf{p}^{\text{T}}\mathbf{b}` are scalars.
then,

.. math:: 
  \begin{array}{l}
  (\mathbf{p}^{\text{T}}\mathbf{A}\mathbf{x}_{0})^{\text{T}}=\mathbf{p}^{\text{T}}\mathbf{A}\mathbf{x}_{0}\\
  ({(\mathbf{x}_{0})}^{\text{T}}\mathbf{A}\mathbf{p})^{\text{T}}={(\mathbf{x}_{0})}^{\text{T}}\mathbf{A}\mathbf{p}\\
  (\mathbf{p}^{\text{T}}\mathbf{b})^{\text{T}}=\mathbf{p}^{\text{T}}\mathbf{b}
  \end{array}
  
-
  
.. math:: 
  \begin{array}{l}
  \cfrac{d\hat{f}}{d\alpha} =\cfrac{1}{2}\mathbf{p}^{\text{T}}\mathbf{A}\mathbf{x}_{0}
  +\cfrac{1}{2} {(\mathbf{x}_{0})}^{\text{T}}\mathbf{A}\mathbf{p}
  +\alpha  {( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}
  -\mathbf{p}^{\text{T}}\mathbf{b}=0\\
  \cfrac{d\hat{f}}{d\alpha} =\mathbf{p}^{\text{T}}\mathbf{A}\mathbf{x}_{0}
  +\alpha  {( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}
  -\mathbf{p}^{\text{T}}\mathbf{b}=0\\
  \alpha  {( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}=\mathbf{p}^{\text{T}}(\mathbf{b}-\mathbf{A}\mathbf{x}_{0})\\
  \alpha=\cfrac{\mathbf{p}^{\text{T}}(\mathbf{b}-\mathbf{A}\mathbf{x}_{0})}
  {{( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}} 
  \end{array}
  
-
  
.. math::
  \begin{array}{l}
  \mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}\\
  -\nabla f(\mathbf{x}_{0})=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}=\mathbf{r}_{0}\\
  \end{array}   
  
-
  
.. math::
  \begin{array}{l}
  \alpha=\cfrac{\mathbf{p}^{\text{T}}(\mathbf{b}-\mathbf{A}\mathbf{x}_{0})}
  {{( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}} \\
  \alpha=\cfrac{\mathbf{p}^{\text{T}}(-\nabla f(\mathbf{x}_{0}))}
  {{( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}} \\
  \alpha=\cfrac{\mathbf{p}^{\text{T}}(\mathbf{r}_{0})}
  {{( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}} \\
  \end{array}  
  
Inner product
  
.. math::
  \left \langle \mathbf{x},\mathbf{y} \right \rangle= \mathbf{x}^{\text{T}}\mathbf{y}\\ 

then

.. math::
  \alpha=\cfrac{\mathbf{p}^{\text{T}}(\mathbf{r}_{0})}
  {{( \mathbf{p})}^{\text{T}}\mathbf{A}\mathbf{p}}=
  \cfrac{\left \langle \mathbf{p},\mathbf{r}_{0} \right \rangle}{\left \langle \mathbf{p},\mathbf{A}\mathbf{p} \right \rangle}
  =\cfrac{\left \langle \mathbf{r}_{0},\mathbf{p} \right \rangle}{\left \langle \mathbf{A}\mathbf{p},\mathbf{p} \right \rangle} 
  
We can also write it as

1. Given any initial value :math:`\mathbf{x}_{0}` calculate the residual :math:`-\nabla f(\mathbf{x}_{0})=\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}`.
2. Go along the selected direction :math:`\mathbf{p}_{0}=\mathbf{r}_{0}=-\nabla f(\mathbf{x}_{0})`, calculate

.. math::
  \begin{array}{l}
  \alpha_{0} =\cfrac{(\mathbf{r}_{0},\mathbf{p}_{0})}{(\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})} \\
  \mathbf{x}_{1}= \mathbf{x}_{0}+\alpha_{0} \ast \mathbf{p}_{0}\\
  \end{array}

3. Go along the selected direction :math:`\mathbf{p}_{1}=\mathbf{r}_{1}=-\nabla f(\mathbf{x}_{1})=\mathbf{b}-\mathbf{A}\mathbf{x}_{1}`, calculate

.. math::  
  \begin{array}{l}
  \alpha_{1} =\cfrac{(\mathbf{r}_{1},\mathbf{p}_{1})}{(\mathbf{A}\mathbf{p}_{1},\mathbf{p}_{1})} \\
  \mathbf{x}_{2}= \mathbf{x}_{1}+\alpha_{1} \ast \mathbf{p}_{1}\\
  \end{array} 
  
4. Go along the selected direction :math:`\mathbf{p}_{k}=\mathbf{r}_{k}=-\nabla f(\mathbf{x}_{k})=\mathbf{b}-\mathbf{A}\mathbf{x}_{k}`, calculate  

.. math::  
  \begin{array}{l}
  \alpha_{k} =\cfrac{(\mathbf{r}_{k},\mathbf{p}_{k})}{(\mathbf{A}\mathbf{p}_{k},\mathbf{p}_{k})} \\
  \mathbf{x}_{k+1}= \mathbf{x}_{k}+\alpha_{k} \ast \mathbf{p}_{k}\\
  \end{array}   

A-orthogonality
------------------------
Let :math:`\mathbf{A}` be a symmetric and positive definite matrix. Two vectors :math:`\mathbf{u}` and :math:`\mathbf{v}` are :math:`\mathbf{A}`-orthogonal if

.. math::
  \mathbf{u}^{\text{T}}\mathbf{A}\mathbf{v}=0,\quad \mathbf{u}\ne\mathbf{v}

Suppose that :math:`\{\mathbf{p}_{k}\}` is a sequence of :math:`n` mutually conjugate directions. Then the :math:`\mathbf{p}_{k}` form a basis of :math:`\mathbf{R}^{n}`, so
we can expand the solution :math:`\mathbf{x}_{*}` of :math:`\mathbf{A}\mathbf{a}=\mathbf{b}` in this basis:

.. math::
  \mathbf{x}_{*} = \sum_{i=1}^{n}\alpha_{i}\mathbf{p}_{i}
  
The coefficients are given by

.. math::
  \mathbf{b}= \mathbf{A}\mathbf{x}_{*} = \sum_{i=1}^{n}\alpha_{i}\mathbf{A}\mathbf{p}_{i}

The conjugate gradient method 
-----------------------------------

Matrix :math:`\mathbf{A}` is symmetric and positive definite, we say that two non-zero vectors :math:`\mathbf{u}` and :math:`\mathbf{v}` are conjugate (with respect to :math:`\mathbf{A}`) if

.. math::
  (\mathbf{u},\mathbf{v})_{\mathbf{A}}=(\mathbf{A}\mathbf{u},\mathbf{v})=\mathbf{v}^{\text{T}}\mathbf{A}\mathbf{u}=
  \sum_{j=1}^{n}\sum_{k=1}^{n}a_{jk}u_{j}v{k}=0
  
So, two vectors are conjugate if they are orthogonal with respect to this inner product. Being conjugate
is a symmetric relation: if :math:`\mathbf{u}` is conjugate to :math:`\mathbf{v}`, then :math:`\mathbf{v}` is conjugate to :math:`\mathbf{u}`. (Note: This notion of conjugate is
not related to the notion of complex conjugate.)  

1. Given any initial value :math:`\mathbf{x}_{0}` calculate the residual :math:`\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}`.
2. Go along the selected direction :math:`\mathbf{p}_{0}=\mathbf{r}_{0}`, calculate

.. math::
  \mathbf{x}_{1}= \mathbf{x}_{0}+\alpha_{0} \ast \mathbf{p}_{0}, 
  \quad \alpha_{0} =\cfrac{(\mathbf{r}_{0},\mathbf{p}_{0})}{(\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})} 
  
3. A new residual :math:`\mathbf{r}_{1}=\mathbf{b}-\mathbf{A}\mathbf{x}_{1}` is obtained from :math:`\mathbf{x}_{1}`, and based on :math:`\mathbf{r}_{1}`, A projection is made on :math:`\mathbf{p}_{0}` to obtain a new forward direction :math:`\mathbf{p}_{1}`  

.. math::
  \begin{array}{l}
  \mathbf{p}_{1}= \mathbf{r}_{1}-\beta_{1} \ast \mathbf{p}_{0}, 
  \quad \beta_{1} =\cfrac{(\mathbf{A}\mathbf{p}_{0},\mathbf{r}_{1})}{(\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})} 
  \end{array}

  
3. Repeat  

.. math::
  \mathbf{x}_{k+1}= \mathbf{x}_{k}+\alpha_{k} \ast \mathbf{p}_{k}\to \mathbf{r}_{k+1}\to \mathbf{p}_{k+1}\to \mathbf{x}_{k+2}
  
-
  
.. math::
  \left\{\begin{array}{ll}
  \mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}, & \mathbf{p}_{0}=\mathbf{r}_{0}\\
  \mathbf{x}_{k+1}= \mathbf{x}_{k}+\alpha_{k} \ast \mathbf{p}_{k}\to &\mathbf{r}_{k+1}=\mathbf{p}_{k+1}
  \end{array}\right.
  \Rightarrow 
  \left\{\begin{array}{ll}
  \mathbf{r}_{0},\mathbf{r}_{1},\cdots,\mathbf{r}_{n}\\
  \mathbf{p}_{0},\mathbf{p}_{1},\cdots,\mathbf{p}_{n}\\
  \end{array}\right.  
  
Suppose we want to find the minimizer of an objective function, having the quadratic form:

.. math::
  f(\mathbf{x})=\cfrac{1}{2} \mathbf{x}^{\text{T}}\mathbf{A}\mathbf{x}-\mathbf{x}^{\text{T}} \mathbf{b}
  
Gradient

.. math::
  \begin{array}{l}
  \nabla f=\text{grad } f=\mathbf{A}\mathbf{x}-\mathbf{b}\\
  -\nabla f=-\text{grad } f=\mathbf{b}-\mathbf{A}\mathbf{x}\\
  \end{array}

We can also write it as

1. Given any initial value :math:`\mathbf{x}_{0}` calculate the residual :math:`-\nabla f(\mathbf{x}_{0})=\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}`.
2. Go along the selected direction :math:`\mathbf{p}_{0}=\mathbf{r}_{0}=-\nabla f(\mathbf{x}_{0})`, calculate

.. math::
  \begin{array}{l}
  \alpha_{0} =\cfrac{(\mathbf{r}_{0},\mathbf{p}_{0})}{(\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})} \\
  \mathbf{x}_{1}= \mathbf{x}_{0}+\alpha_{0} \ast \mathbf{p}_{0}\\
  \end{array}
  
3. A new residual :math:`-\nabla f(\mathbf{x}_{1})=\mathbf{r}_{1}=\mathbf{b}-\mathbf{A}\mathbf{x}_{1}` is obtained from :math:`\mathbf{x}_{1}`, and based on :math:`\mathbf{r}_{1}`, A projection is made on :math:`\mathbf{p}_{0}` to obtain a new forward direction :math:`\mathbf{p}_{1}`    

.. math::
  \begin{array}{l}
  \mathbf{p}_{1}= \mathbf{r}_{1}-\beta_{1} \ast \mathbf{p}_{0}, 
  \quad \beta_{1} =\cfrac{(\mathbf{A}\mathbf{p}_{0},\mathbf{r}_{1})}{(\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})} 
  \end{array}
  
-
  
.. math::
  \begin{array}{l}
  \alpha_{1} =\cfrac{(\mathbf{r}_{1},\mathbf{p}_{1})}{(\mathbf{A}\mathbf{p}_{1},\mathbf{p}_{1})} \\
  \mathbf{x}_{2}= \mathbf{x}_{1}+\alpha_{1} \ast \mathbf{p}_{1}\\
  \end{array}  
  

4. A new residual :math:`-\nabla f(\mathbf{x}_{k})=\mathbf{r}_{k}=\mathbf{b}-\mathbf{A}\mathbf{x}_{k}` is obtained from :math:`\mathbf{x}_{k}`, and based on :math:`\mathbf{r}_{k}`, A projection is made on :math:`\mathbf{p}_{k-1}` to obtain a new forward direction :math:`\mathbf{p}_{k}`    

.. math::
  \begin{array}{l}
  \mathbf{p}_{k}= \mathbf{r}_{k}-\beta_{k} \ast \mathbf{p}_{k-1}, 
  \quad \beta_{k} =\cfrac{(\mathbf{A}\mathbf{p}_{k-1},\mathbf{r}_{k})}{(\mathbf{A}\mathbf{p}_{k-1},\mathbf{p}_{k-1})} 
  \end{array}
  
-
  
.. math::
  \begin{array}{l}
  \alpha_{k} =\cfrac{(\mathbf{r}_{k},\mathbf{p}_{k})}{(\mathbf{A}\mathbf{p}_{k},\mathbf{p}_{k})} \\
  \mathbf{x}_{k+1}= \mathbf{x}_{k}+\alpha_{k} \ast \mathbf{p}_{k}\\
  \end{array}
  
Example 1  
-----------------------------------------------------------------
Solve :math:`\mathbf{A}\mathbf{x}=\mathbf{b}`

where

.. math:: 
  \mathbf{A} =\begin{bmatrix}
  3&2 \\
  2&6
  \end{bmatrix}\quad
  \mathbf{x}=\begin{bmatrix}
  x_{1}\\x_{2}
  \end{bmatrix}\quad
  \mathbf{b}=\begin{bmatrix}
  2\\-8
  \end{bmatrix}
  
The quadratic form:

.. math:: 
  f(\mathbf{x})=\cfrac{1}{2}\mathbf{x}^{\text{T}}\mathbf{A}\mathbf{x}-\mathbf{b}^{\text{T}}\mathbf{x}
  
-
  
.. math:: 
  \begin{array}{l}
  f(\mathbf{x})=\cfrac{1}{2}\begin{bmatrix}
    x_{1}&x_{2}
  \end{bmatrix}
   \begin{bmatrix}
    a_{11}& a_{12}\\
    a_{21}&a_{22}
  \end{bmatrix}\begin{bmatrix}
   x_{1}\\x_{2}
  \end{bmatrix}
  -\begin{bmatrix}
    b_{1}&b_{2}
  \end{bmatrix}\begin{bmatrix}
   x_{1}\\x_{2}
  \end{bmatrix}\\
  f(\mathbf{x})=\cfrac{1}{2}\begin{bmatrix}
    x_{1}&x_{2}
  \end{bmatrix}\begin{bmatrix}
   a_{11}x_{1}+a_{12}x_{2}\\a_{21}x_{1}+a_{22}x_{2}
  \end{bmatrix}-b_{1}x_{1}-b_{2}x_{2}\\
  f(\mathbf{x})=\cfrac{1}{2}(a_{11}x_{1}^{2}+a_{12}x_{1}x_{2}+a_{21}x_{1}x_{2}+a_{22}x_{2}^{2})-b_{1}x_{1}-b_{2}x_{2}
  \end{array}
  
In general
  
.. math:: 
  f(\mathbf{x})=\cfrac{1}{2}a_{ij}x_{i}x_{j}-b_{j}x_{j}  
  
The exact solution

.. math:: 
  \mathbf{x}=\begin{bmatrix}
  x_{1}\\x_{2}
  \end{bmatrix}=
  \begin{bmatrix}
  2\\-2
  \end{bmatrix}
  
  
using steepest descent method
`````````````````````````````````````````
1. Given any initial value :math:`\mathbf{x}_{0}` calculate the residual :math:`-\nabla f(\mathbf{x}_{0})=\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}`.

.. math:: 
  \mathbf{x}_{0}=
  \begin{bmatrix}
  0\\0
  \end{bmatrix}
  
-
  
.. math::   
  -\nabla f(\mathbf{x}_{0})=\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}=\mathbf{b}=
  \begin{bmatrix}
  2\\-8
  \end{bmatrix}
  
2. Go along the selected direction :math:`\mathbf{p}_{0}=\mathbf{r}_{0}=-\nabla f(\mathbf{x}_{0})`, calculate

.. math::
  \begin{array}{l}
  \alpha_{0} =\cfrac{(\mathbf{r}_{0},\mathbf{p}_{0})}{(\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})} \\
  \mathbf{r}_{0}=[2,-8]^{\text{T}}\\
  \mathbf{A}\mathbf{p}_{0}=[-10,44]^{\text{T}}\\
  (\mathbf{r}_{0},\mathbf{p}_{0})=2*2+8*8=68\\
  (\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})=332\\
  \alpha_{0}=68/332=0.2048\\
  \mathbf{x}_{1}= \mathbf{x}_{0}+\alpha_{0} \ast \mathbf{p}_{0}=[0.4096,-1.6386]^{\text{T}}\\
  \end{array}
  
3. Go along the selected direction :math:`\mathbf{p}_{1}=\mathbf{r}_{1}=-\nabla f(\mathbf{x}_{1})=\mathbf{b}-\mathbf{A}\mathbf{x}_{1}=[4.0482,1.0120]^{\text{T}}`, calculate

.. math::  
  \begin{array}{l}
  \alpha_{1} =\cfrac{(\mathbf{r}_{1},\mathbf{p}_{1})}{(\mathbf{A}\mathbf{p}_{1},\mathbf{p}_{1})} \\
  \mathbf{x}_{2}= \mathbf{x}_{1}+\alpha_{1} \ast \mathbf{p}_{1}\\
  \end{array}   
  
using conjugate gradient method
`````````````````````````````````````````
1. Given any initial value :math:`\mathbf{x}_{0}` calculate the residual :math:`-\nabla f(\mathbf{x}_{0})=\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}`.

.. math:: 
  \mathbf{x}_{0}=
  \begin{bmatrix}
  0\\0
  \end{bmatrix}
  
-
  
.. math::   
  -\nabla f(\mathbf{x}_{0})=\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}=\mathbf{b}=
  \begin{bmatrix}
  2\\-8
  \end{bmatrix}
  
2. Go along the selected direction :math:`\mathbf{p}_{0}=\mathbf{r}_{0}=-\nabla f(\mathbf{x}_{0})=[2,-8]^{\text{T}}`, calculate

.. math::
  \begin{array}{l}
  \alpha_{0} =\cfrac{(\mathbf{r}_{0},\mathbf{p}_{0})}{(\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})} \\
  \mathbf{r}_{0}=[2,-8]^{\text{T}}\\
  \mathbf{A}\mathbf{p}_{0}=[-10,44]^{\text{T}}\\
  (\mathbf{r}_{0},\mathbf{p}_{0})=2*2+8*8=68\\
  (\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})=332\\
  \alpha_{0}=68/332=0.2048\\
  \mathbf{x}_{1}= \mathbf{x}_{0}+\alpha_{0} \ast \mathbf{p}_{0}=[0.4096,-1.6386]^{\text{T}}\\
  \end{array}
  
3. A new residual :math:`-\nabla f(\mathbf{x}_{1})=\mathbf{r}_{1}=\mathbf{b}-\mathbf{A}\mathbf{x}_{1}=[4.0482,1.0120]^{\text{T}}` is obtained from :math:`\mathbf{x}_{1}`, and based on :math:`\mathbf{r}_{1}`, A projection is made on :math:`\mathbf{p}_{0}` to obtain a new forward direction :math:`\mathbf{p}_{1}`    

.. math::
  \begin{array}{l}
  \quad \beta_{1} =\cfrac{(\mathbf{A}\mathbf{p}_{0},\mathbf{r}_{1})}{(\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})}=-0.2561 \\
  \mathbf{p}_{1}= \mathbf{r}_{1}-\beta_{1} \ast \mathbf{p}_{0} =[4.5603,-1.0364]^{\text{T}}
  \end{array}
  
-
  
.. math::
  \begin{array}{l}
  \alpha_{1} =\cfrac{(\mathbf{r}_{1},\mathbf{p}_{1})}{(\mathbf{A}\mathbf{p}_{1},\mathbf{p}_{1})}=0.3487 \\
  \mathbf{x}_{2}= \mathbf{x}_{1}+\alpha_{1} \ast \mathbf{p}_{1}=[2,-2]^{\text{T}}\\
  \end{array}    
  
 
Solving the 1D Poisson equation using conjugate gradient method
-----------------------------------------------------------------
Consider the 1D Poisson equation

.. math::
  \cfrac{d^{2}u}{dx^{2}} =f(x)=-1
  
on :math:`\Omega=[0,1]` with boundary conditions
  
.. math:: 
  \begin{array}{l}
  u(0)=0\\
  u'(1)=\cfrac{du}{dx}\bigg|_{x=1}=0 
  \end{array}
  
which has analytical solution

.. math:: 
  u=x-\cfrac{1}{2}x^{2}

We we will use this specific example to investigate various approaches to solving partial differential equations with finite differences, in which we discretize the domain by defining :math:`N` equally
spaced points  

.. math:: 
  \begin{array}{l}
  \cfrac{u_{i+1}-2u_{i}+u_{i-1}}{\Delta x^{2}} =f_{i}\\
  {u_{i+1}-2u_{i}+u_{i-1}}={\Delta x^{2}}f_{i}\\
  -u_{i-1}+2u_{i}-u_{i+1}=-{\Delta x^{2}}f_{i}={\Delta x^{2}}\\
  \end{array}  

The boundary conditions require special care. For :math:`x = 0` we have a Dirichlet boundary condition
which allows us to fix the value :math:`u_{1} = 0`. For :math:`x = 1` we have a Neumann boundary condition
:math:`du/dx = 0`. This is a symmetry boundary condition, so that in this case we can imagine a ’ghost’
point :math:`u_{N+1}` which is always equal to :math:`u_{N-1}`. This leads to the expression for point :math:`x_{N}`:  


.. math:: 
  \begin{array}{l}
  u_{N}=u_{N-1}+\cfrac{1}{2}\Delta x^{2}\\
  \end{array}
  
  
then

.. math:: 
  \begin{array}{l}
  -u_{i-1}+2u_{i}-u_{i+1}={\Delta x^{2}}\\
  -u_{1}+2u_{2}-u_{3}={\Delta x^{2}}\\
  -u_{2}+2u_{3}-u_{4}={\Delta x^{2}}\\
  \cdots\\
 -u_{N-3}+2u_{N-2}-u_{N-1}={\Delta x^{2}}\\
 -u_{N-2}+2u_{N-1}-u_{N}={\Delta x^{2}}\\
  \end{array}   
  
-
  
.. math::   
  \begin{array}{l}
  0+2u_{2}-u_{3}={\Delta x^{2}}\\
  -u_{2}+2u_{3}-u_{4}={\Delta x^{2}}\\
  \cdots\\
  -u_{N-3}+2u_{N-2}-u_{N-1}={\Delta x^{2}}\\
  -u_{N-2}+2u_{N-1}-(u_{N-1}+\cfrac{1}{2}\Delta x^{2} )={\Delta x^{2}}\\
  \end{array}  
  
-
  
.. math::   
  \begin{array}{l}
  0+2u_{2}-u_{3}={\Delta x^{2}}\\
  -u_{2}+2u_{3}-u_{4}={\Delta x^{2}}\\
  \cdots\\
 -u_{N-3}+2u_{N-2}-u_{N-1}={\Delta x^{2}}\\
 -u_{N-2}+u_{N-1}=\cfrac{3}{2}{\Delta x^{2}}\\
  \end{array}    
    
-
  
.. math::   
  \begin{array}{l}
  2u_{2}-u_{3}={\Delta x^{2}}\\
  -u_{2}+2u_{3}-u_{4}={\Delta x^{2}}\\
  \cdots\\
 -u_{N-3}+2u_{N-2}-u_{N-1}={\Delta x^{2}}\\
 -u_{N-2}+u_{N-1}=\cfrac{3}{2}{\Delta x^{2}}\\
  \end{array}    

  
This can be written as  
  
.. math::     
  \mathbf{A}\mathbf{u}=\mathbf{b}  
  
where(for :math:`N=5`)

.. math::  
  \mathbf{A}=
  \begin{bmatrix}
  2& -1 & 0& 0& 0\\
  -1& 2 & -1& 0& 0\\
  0&-1& 2 & -1& 0\\
  0&0&-1& 2 & -1\\
   0&0 &0 &-1  &1
  \end{bmatrix}  
  
-
  
.. math::   
  \mathbf{u}=\begin{bmatrix}
   u_{2}\\u_{3}\\u_{4}\\u_{5}\\u_{6}
  \end{bmatrix}=\mathbf{v}=
  \begin{bmatrix}
  v_{1}\\v_{2}\\v_{3}\\v_{4}\\v_{5}
  \end{bmatrix}  
  
-
  
.. math:: 
  \mathbf{b}=
  \begin{bmatrix}
  \Delta x^{2}\\\Delta x^{2}\\\Delta x^{2}\\\Delta x^{2}\\\cfrac{3}{2}\Delta x^{2}
  \end{bmatrix}
  
-
  
.. math:: 
  \begin{array}{l}
  1&\underbrace{\Delta x}& 2&\underbrace{\Delta x}&3&\underbrace{\Delta x}&4&\underbrace{\Delta x}&5&\underbrace{\Delta x}&6&\underbrace{\Delta x}&7\\
  &&v_{1}&&v_{2}&&v_{3}&&v_{4}&&v_{5}\\
  \end{array}
  
-
  
.. math:: 
  \begin{array}{l}
  N_{\text{Points}}=7\\
  N=N_{\text{Points}}-2=5\\
  \Delta x=\cfrac{1}{N_{\text{Points}}-1} =\cfrac{1}{7-1}=\cfrac{1}{6}\\
  \Delta x^{2}=\cfrac{1}{36}\\
  \cfrac{3}{2}\Delta x^{2}=\cfrac{3}{72}\\
  \end{array} 
  
-
  
.. math:: 
  \begin{array}{l}
  x_{1}=0\Delta x=0\\
  x_{2}=1\Delta x=1/6\\
  x_{3}=2\Delta x=2/6\\
  x_{4}=3\Delta x=3/6\\
  x_{5}=4\Delta x=4/6\\
  x_{6}=5\Delta x=5/6\\
  x_{7}=6\Delta x=1\\
  \end{array}  
  
  
   
Solve :math:`\mathbf{A}\mathbf{x}=\mathbf{b}`

where

.. math::  
  \mathbf{A}=
  \begin{bmatrix}
  2& -1 & 0& 0& 0\\
  -1& 2 & -1& 0& 0\\
  0&-1& 2 & -1& 0\\
  0&0&-1& 2 & -1\\
   0&0 &0 &-1  &1
  \end{bmatrix}  
  
-
  
.. math::   
  \mathbf{v}=
  \begin{bmatrix}
  v_{1}\\v_{2}\\v_{3}\\v_{4}\\v_{5}
  \end{bmatrix}    
  
-
  
.. math::   
  \mathbf{b}=
  \begin{bmatrix}
  \Delta x^{2}\\\Delta x^{2}\\\Delta x^{2}\\\Delta x^{2}\\\cfrac{1}{2}\Delta x^{2}
  \end{bmatrix}
  =  \begin{bmatrix}
  \cfrac{1}{36} \\\cfrac{1}{36} \\\cfrac{1}{36} \\\cfrac{1}{36} \\\cfrac{1}{72} 
  \end{bmatrix}  
  

The quadratic form:

.. math:: 
  f(\mathbf{v})=\cfrac{1}{2}\mathbf{v}^{\text{T}}\mathbf{A}\mathbf{v}-\mathbf{b}^{\text{T}}\mathbf{v}
  
The exact soltion:

.. math:: 
  u(x)=x-\cfrac{1}{2}x^{2}\\
  
-
  
.. math:: 
  \mathbf{u}=\begin{bmatrix}
  u(0)\\u(1/6)\\u(2/6)\\u(3/6)\\u(4/6)\\u(5/6)\\u(1)
  \end{bmatrix}
  =\begin{bmatrix}
  0\\11/72\\20/72\\27/72\\32/72\\35/72\\36/72
  \end{bmatrix}
  =\begin{bmatrix}
  0\\0.1528\\0.2778\\0.3750\\0.4444\\0.4861\\0.5000
  \end{bmatrix}
  

another form
-----------------------

.. math:: 
  \cfrac{u_{i-1}-2u_{i}+u_{i+1}}{\Delta x^{2}} =f_{i}=-1

-
  
.. math::   
  \begin{array}{l}
  \cfrac{u_{1}-2u_{2}+u_{3}}{\Delta x^{2}} =f_{2}=-1\\
  \cfrac{u_{2}-2u_{3}+u_{4}}{\Delta x^{2}} =f_{3}=-1\\
  \cdots\\
  \cfrac{u_{i-1}-2u_{i}+u_{i+1}}{\Delta x^{2}} =f_{i}=-1\\
  \cdots\\
  \cfrac{u_{N-3}-2u_{N-2}+u_{N-1}}{\Delta x^{2}} =f_{N-2}=-1\\
  \cfrac{u_{N-2}-2u_{N-1}+u_{N}}{\Delta x^{2}} =f_{N-1}=-1\\
  \end{array}
  
-
  
.. math:: 
  \begin{array}{l}
  u_{1}=0\\
  u_{N}=u_{N-1}+\cfrac{1}{2}\Delta x^{2} \\
  \end{array}  
  
  
-
  
.. math:: 
  \begin{array}{l}
  u_{1}=0\\
  \cfrac{u_{1}-2u_{2}+u_{3}}{\Delta x^{2}} =f_{2}=-1\\
  \cfrac{-2u_{2}+u_{3}}{\Delta x^{2}} =f_{2}-\cfrac{u_{1}}{\Delta x^{2}}=-1-\cfrac{u_{1}}{\Delta x^{2}}=-1\\
  \cfrac{-2u_{2}+u_{3}}{\Delta x^{2}} =\hat{f}_{2}=f_{2}-\cfrac{u_{1}}{\Delta x^{2}}=-1-\cfrac{u_{1}}{\Delta x^{2}}=-1\\
  \end{array}  

-
  
.. math:: 
  \begin{array}{l}
  u_{N}=u_{N-1}+\cfrac{1}{2}\Delta x^{2} \\
  \cfrac{u_{N-2}-2u_{N-1}+u_{N}}{\Delta x^{2}} =f_{N-1}=-1\\
  \cfrac{u_{N-2}-2u_{N-1}+u_{N-1}+\cfrac{1}{2}\Delta x^{2}}{\Delta x^{2}} =f_{N-1}=-1\\
  \cfrac{u_{N-2}-u_{N-1}}{\Delta x^{2}} =\hat{f}_{N-1}=f_{N-1}-\cfrac{1}{2}=-\cfrac{3}{2}\\
  \end{array}  
  
Let

.. math::
  \mathbf{u}=\begin{bmatrix}
  u_{2}\\u_{3}\\\vdots\\u_{N-2}\\u_{N-1}
  \end{bmatrix}
  
then

.. math::
  \mathbf{A}\mathbf{u}=\mathbf{f}
  
where

.. math::
  \mathbf{f}=\begin{bmatrix}
  \hat{f}_{2}\\f_{3}\\\vdots\\f_{N-2}\\\hat{f}_{N-1}
  \end{bmatrix}=\begin{bmatrix}
  -1\\-1\\\vdots\\-1\\-3/2
  \end{bmatrix}

-
  
.. math::
  \mathbf{A}\mathbf{u}=\begin{bmatrix}
  \cfrac{-2u_{2}+u_{3}}{\Delta x^{2}}\\
  \cfrac{u_{2}-2u_{3}+u_{4}}{\Delta x^{2}}\\
  \vdots\\
  \cfrac{u_{N-3}-2u_{N-2}+u_{N-1}}{\Delta x^{2}}\\
  \cfrac{u_{N-2}-u_{N-1}}{\Delta x^{2}}
  \end{bmatrix}  
  
using the Method of Steepest Descent

Let :math:`\mathbf{b}=\mathbf{f}`, then
  
.. math::
  g(\mathbf{u})=\cfrac{1}{2} \mathbf{u}^{\text{T}}\mathbf{A}\mathbf{u}-\mathbf{b}^{\text{T}}\mathbf{u}

1. Given any initial value :math:`\mathbf{u}_{0}` calculate the residual :math:`-\nabla g(\mathbf{u}_{0})=\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{u}_{0}`.
2. Go along the selected direction :math:`\mathbf{p}_{0}=\mathbf{r}_{0}=-\nabla g(\mathbf{u}_{0})`, calculate

.. math::
  \begin{array}{l}
  \alpha_{0} =\cfrac{(\mathbf{r}_{0},\mathbf{p}_{0})}{(\mathbf{A}\mathbf{p}_{0},\mathbf{p}_{0})} \\
  \mathbf{u}_{1}= \mathbf{u}_{0}+\alpha_{0} \ast \mathbf{p}_{0}\\
  \end{array}
  
where 

.. math::
  \mathbf{A}\mathbf{p}_{0}=\begin{bmatrix}
  \cfrac{-2p_{2}^{(0)}+p_{3}^{(0)}}{\Delta x^{2}}\\
  \cfrac{p_{2}^{(0)}-2p_{3}^{(0)}+p_{4}^{(0)}}{\Delta x^{2}}\\
  \vdots\\
  \cfrac{p_{N-3}^{(0)}-2p_{N-2}^{(0)}+p_{N-1}^{(0)}}{\Delta x^{2}}\\
  \cfrac{p_{N-2}^{(0)}-p_{N-1}^{(0)}}{\Delta x^{2}}
  \end{bmatrix}  

3. Go along the selected direction :math:`\mathbf{p}_{1}=\mathbf{r}_{1}=-\nabla g(\mathbf{u}_{1})=\mathbf{b}-\mathbf{A}\mathbf{u}_{1}`, calculate

.. math::  
  \begin{array}{l}
  \alpha_{1} =\cfrac{(\mathbf{r}_{1},\mathbf{p}_{1})}{(\mathbf{A}\mathbf{p}_{1},\mathbf{p}_{1})} \\
  \mathbf{u}_{2}= \mathbf{u}_{1}+\alpha_{1} \ast \mathbf{p}_{1}\\
  \end{array} 
  
where   
  
.. math::
  \mathbf{A}\mathbf{p}_{1}=\begin{bmatrix}
  \cfrac{-2p_{2}^{(1)}+p_{3}^{(1)}}{\Delta x^{2}}\\
  \cfrac{p_{2}^{(1)}-2p_{3}^{(1)}+p_{4}^{(1)}}{\Delta x^{2}}\\
  \vdots\\
  \cfrac{p_{N-3}^{(1)}-2p_{N-2}^{(1)}+p_{N-1}^{(1)}}{\Delta x^{2}}\\
  \cfrac{p_{N-2}^{(1)}-p_{N-1}^{(1)}}{\Delta x^{2}}
  \end{bmatrix}    
  
4. Go along the selected direction :math:`\mathbf{p}_{k}=\mathbf{r}_{k}=-\nabla g(\mathbf{u}_{k})=\mathbf{b}-\mathbf{A}\mathbf{u}_{k}`, calculate  

.. math::  
  \begin{array}{l}
  \alpha_{k} =\cfrac{(\mathbf{r}_{k},\mathbf{p}_{k})}{(\mathbf{A}\mathbf{p}_{k},\mathbf{p}_{k})} \\
  \mathbf{u}_{k+1}= \mathbf{u}_{k}+\alpha_{k} \ast \mathbf{p}_{k}\\
  \end{array}   
  
where   
  
.. math::
  \mathbf{A}\mathbf{p}_{k}=\begin{bmatrix}
  \cfrac{-2p_{2}^{(k)}+p_{3}^{(k)}}{\Delta x^{2}}\\
  \cfrac{p_{2}^{(k)}-2p_{3}^{(k)}+p_{4}^{(k)}}{\Delta x^{2}}\\
  \vdots\\
  \cfrac{p_{N-3}^{(k)}-2p_{N-2}^{(k)}+p_{N-1}^{(k)}}{\Delta x^{2}}\\
  \cfrac{p_{N-2}^{(k)}-p_{N-1}^{(k)}}{\Delta x^{2}}
  \end{bmatrix}     
  
The conjugate gradient method from code
----------------------------------------------
1. Given any initial value :math:`\mathbf{x}_{0}` calculate the residual :math:`-\nabla f(\mathbf{x}_{0})=\mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}`.
2. Go along the selected direction :math:`\mathbf{p}_{0}=\mathbf{r}_{0}=-\nabla f(\mathbf{x}_{0})`, calculate

.. math::
  \begin{array}{l}
  \mathbf{r}_{0}=\mathbf{b}-\mathbf{A}\mathbf{x}_{0}\\
  \mathbf{p}_{0}=\mathbf{r}_{0}\\
  \delta_{new}=(\mathbf{r}_{0},\mathbf{r}_{0})\\
  \text{while_loop:}\\
  \alpha_{k} =\cfrac{(\mathbf{r}_{k},\mathbf{r}_{k})}{(\mathbf{A}\mathbf{p}_{k},\mathbf{p}_{k})} \\
  \mathbf{x}_{k+1}= \mathbf{x}_{k}+\alpha_{k} \ast \mathbf{p}_{k}\\
  \mathbf{r}_{k+1}=\mathbf{b}-\mathbf{A}\mathbf{x}_{k+1}\\
  \delta_{old}=\delta_{new}\\
  \delta_{new}=(\mathbf{r}_{k+1},\mathbf{r}_{k+1})\\
  \beta_{k}=\cfrac{\delta_{new}}{\delta_{old}}=\cfrac{(\mathbf{r}_{k+1},\mathbf{r}_{k+1})}{(\mathbf{r}_{k},\mathbf{r}_{k})}\\
  \mathbf{p}_{k+1}=\mathbf{r}_{k+1}+\beta_{k}\mathbf{p}_{k}
  \end{array}
  
-
  
.. math::  
  \begin{array}{l}
  \mathbf{r}_{k}=\mathbf{b}-\mathbf{A}\mathbf{x}_{k}\\
  \mathbf{r}_{k+1}=\mathbf{b}-\mathbf{A}\mathbf{x}_{k+1}\\
  \mathbf{r}_{k+1}-\mathbf{r}_{k}=-\mathbf{A}(\mathbf{x}_{k+1}-\mathbf{x}_{k})\\
  \mathbf{x}_{k+1}= \mathbf{x}_{k}+\alpha_{k} \ast \mathbf{p}_{k}\\
  \mathbf{r}_{k+1}-\mathbf{r}_{k}=-\mathbf{A}(\mathbf{x}_{k+1}-\mathbf{x}_{k})=-\alpha_{k} \mathbf{A} \mathbf{p}_{k}\\
  \mathbf{r}_{k+1}-\mathbf{r}_{k}=-\alpha_{k} \mathbf{A} \mathbf{p}_{k}\\
  \mathbf{r}_{k+1}=\mathbf{r}_{k}-\alpha_{k} \mathbf{A} \mathbf{p}_{k}\\
  \end{array}  
  
Solving the 2D Poisson equation
-------------------------------------
Consider the 2D Poisson equation

.. math::
  \cfrac{\partial^{2}u}{\partial x^{2}}+\cfrac{\partial^{2}u}{\partial x^{2}} =f(x,y)=-1
  
on :math:`\Omega=[0,1]\times[0,1]` with boundary conditions
  
.. math:: 
  \begin{array}{l}
  u(0,j)=0\\
  u'(1,j)=\cfrac{\partial u}{\partial x}\bigg|_{x=1}=0 
  \end{array}
  
which has analytical solution

.. math:: 
  u(x,y)=x-\cfrac{1}{2}x^{2}

We we will use this specific example to investigate various approaches to solving partial differential equations with finite differences, in which we discretize the domain by defining :math:`N` equally
spaced points  

.. math:: 
  \begin{array}{l}
  \cfrac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{\Delta x^{2}} 
  + \cfrac{u_{i,j+1}-2u_{i,j}+u_{i,j-1}}{\Delta y^{2}} =f_{i,j}\\
  \cfrac{u_{i-1,j}-2u_{i,j}+u_{i+1,j}}{\Delta x^{2}} 
  + \cfrac{u_{i,j-1}-2u_{i,j}+u_{i,j+1}}{\Delta y^{2}} =f_{i,j}\\
  \end{array}  

The boundary conditions require special care. For :math:`x = 0` we have a Dirichlet boundary condition
which allows us to fix the value :math:`u_{1} = 0`. For :math:`x = 1` we have a Neumann boundary condition
:math:`du/dx = 0`. This is a symmetry boundary condition, so that in this case we can imagine a ’ghost’
point :math:`u_{M+1,j}` which is always equal to :math:`u_{M-1,j}`. This leads to the expression for point :math:`x_{M,j}`  


.. math:: 
  \begin{array}{l}
  u_{M,j}=u_{M-1,j}+\cfrac{1}{2}\Delta x^{2}\\
  \end{array}
  
  
then

.. math:: 
  \begin{array}{l}
  j=2\\
  \cfrac{u_{i-1,2}-2u_{i,2}+u_{i+1,2}}{\Delta x^{2}} 
  + \cfrac{u_{i,2-1}-2u_{i,2}+u_{i,2+1}}{\Delta y^{2}} =f_{i,2}\\
  \cfrac{u_{i-1,2}-2u_{i,2}+u_{i+1,2}}{\Delta x^{2}} 
  + \cfrac{u_{i,1}-2u_{i,2}+u_{i,3}}{\Delta y^{2}} =f_{i,2}\\
  \cfrac{u_{i-1,2}-2u_{i,2}+u_{i+1,2}}{\Delta x^{2}} 
  + \cfrac{-u_{i,2}+u_{i,3}}{\Delta y^{2}} =f_{i,2}\\
  \end{array}  

-
  
.. math:: 
  \begin{array}{l}
  j=N-1\\
  \cfrac{u_{i-1,N-1}-2u_{i,N-1}+u_{i+1,N-1}}{\Delta x^{2}} 
  + \cfrac{u_{i,N-1-1}-2u_{i,N-1}+u_{i,N-1+1}}{\Delta y^{2}} =f_{N-1,2}\\
  \cfrac{u_{i-1,N-1}-2u_{i,N-1}+u_{i+1,N-1}}{\Delta x^{2}} 
  + \cfrac{u_{i,N-2}-2u_{i,N-1}+u_{i,N}}{\Delta y^{2}} =f_{N-1,2}\\
  \cfrac{u_{i-1,N-1}-2u_{i,N-1}+u_{i+1,N-1}}{\Delta x^{2}} 
  + \cfrac{u_{i,N-2}-u_{i,N-1}}{\Delta y^{2}} =f_{N-1,2}\\
  \end{array}
  
-
  
.. math:: 
  \begin{array}{l}
  i=2\\
  \cfrac{u_{2-1,j}-2u_{2,j}+u_{2+1,j}}{\Delta x^{2}} 
  + \cfrac{u_{2,j-1}-2u_{2,j}+u_{2,j+1}}{\Delta y^{2}} =f_{2,j}\\
  \cfrac{u_{1,j}-2u_{2,j}+u_{3,j}}{\Delta x^{2}} 
  + \cfrac{u_{2,j-1}-2u_{2,j}+u_{2,j+1}}{\Delta y^{2}} =f_{2,j}\\
  \cfrac{-2u_{2,j}+u_{3,j}}{\Delta x^{2}} 
  + \cfrac{u_{2,j-1}-2u_{2,j}+u_{2,j+1}}{\Delta y^{2}} =f_{2,j}\\
  \end{array}  
  
-
  
.. math:: 
  \begin{array}{l}
  i=M-1\\
  \cfrac{u_{M-1-1,j}-2u_{M-1,j}+u_{M-1+1,j}}{\Delta x^{2}} 
  + \cfrac{u_{M-1,j-1}-2u_{M-1,j}+u_{M-1,j+1}}{\Delta y^{2}} =f_{M-1,j}\\
  \cfrac{u_{M-2,j}-2u_{M-1,j}+u_{M,j}}{\Delta x^{2}} 
  + \cfrac{u_{M-1,j-1}-2u_{M-1,j}+u_{M-1,j+1}}{\Delta y^{2}} =f_{M-1,j}\\
  \cfrac{u_{M-2,j}-2u_{M-1,j}+u_{M-1,j}+\cfrac{1}{2}\Delta x^{2}}{\Delta x^{2}} 
  + \cfrac{u_{M-1,j-1}-2u_{M-1,j}+u_{M-1,j+1}}{\Delta y^{2}} =f_{M-1,j}\\
  \cfrac{u_{M-2,j}-u_{M-1,j}}{\Delta x^{2}} 
  + \cfrac{u_{M-1,j-1}-2u_{M-1,j}+u_{M-1,j+1}}{\Delta y^{2}} =f_{M-1,j}-\cfrac{1}{2}\\
  \end{array} 