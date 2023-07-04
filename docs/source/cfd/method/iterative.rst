Iterative Solver 
==================================

#. `Scientific Computing || 02 Week 7 19 1 Introduction to spectral methods 10 46 <https://www.youtube.com/watch?v=ymsY8IFbOwY/>`_
#. `Solving the Discrete Poisson Equation using Jacobi, SOR, Conjugate Gradients, and the FFT <https://people.eecs.berkeley.edu/~demmel/cs267/lecture24/lecture24.html>`_
#. `Iteration methods <https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html>`_
#. `Numerical Methods for Elliptic Equations <http://www.fem.unicamp.br/~phoenics/SITE_PHOENICS/Apostilas/CFD-1_U%20Michigan_Hong/Lecture11.pdf>`_
#. `Iterative Methods: Part 2 <https://crunchingnumbers.live/2017/07/09/iterative-methods-part-2/>`_
#. `FFTW Tutorial <https://github.com/jonathanschilling/fftw_tutorial/>`_


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
  