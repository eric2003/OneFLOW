Implicit Compact Pade (ICP) Scheme
==================================

#. `CFD Julia: A Learning Module Structuring an Introductory Course on Computational Fluid Dynamics <https://www.mdpi.com/2311-5521/4/3/159/>`_
#. `Compact finite difference <https://en.wikipedia.org/wiki/Compact_finite_difference>`_
#. `AN INTRODUCTION TO COMPACT FINITE DIFFERENCE <https://people.bath.ac.uk/ensdasr/COMPACT/dasr.compact.pdf>`_
#. `Two-Dimensional Compact-Finite-Difference Schemes for Solving the bi-Laplacian Operator with Homogeneous Wall-Normal Derivatives <https://www.mdpi.com/2227-7390/9/19/2508/>`_

Approximation of First Derivative
----------------------------------------
Given the values of a function on a set of nodes the finite
difference approximation to the derivative of the function is
expressed as a linear combination of the given function
values. For simplicity consider a uniformly spaced mesh
where the nodes are indexed by :math:`i`. The independent variable
at the nodes is :math:`x_{i}=h(i-1)` for :math:`1\le i\le N` and the function
values at the nodes :math:`f_{i}=f(x_{i})` are given. The finite difference
approximation :math:`{f}'_{i}` to the first derivative :math:`(df/dx)(x_{i})` at the node 
:math:`i` depends on the function values at nodes near :math:`i`.
For second- and fourth-order central differences the
approximation :math:`{f}'_{i}` depends on the sets :math:`(f_{i-1},f_{i+1})` and
:math:`(f_{i-2},f_{i-1},f_{i+1},f_{i+2})`, respectively. In the spectral
methods, however, the value off :math:`{f}'_{i}` depends on all the nodal
values. The Pade or compact finite difference schemes mimic this global dependence. The schemes
presented here are generalizations of the Pade scheme.

These generalizations are derived by writing approximations of the form:

.. math::
  \begin{align}
  &\beta{f}'_{i-2}+\alpha{f}'_{i-1}+{f}'_{i}+\alpha{f}'_{i+1}+\beta{f}'_{i+2}\\
    = &c\cfrac{{f}_{i+3}-{f}_{i-3}}{6h}+b\cfrac{{f}_{i+2}-{f}_{i-2}}{4h}+a\cfrac{{f}_{i+1}-{f}_{i-1}}{2h}
  \end{align}
  
Approximation of Second Derivative
----------------------------------------
The derivation of compact approximations for the second
derivative proceeds exactly analogous to the first derivative.
Once again the starting point is a relation of the form

.. math::
  \begin{align}
  &\beta{f}''_{i-2}+\alpha{f}''_{i-1}+{f}''_{i}+\alpha{f}''_{i+1}+\beta{f}''_{i+2}\\
  = &c\cfrac{{f}_{i+3}-2{f}_{i}+{f}_{i-3}}{9h^{2}}+b\cfrac{{f}_{i+2}-2{f}_{i}+{f}_{i-2}}{4h^{2}}+a\cfrac{{f}_{i+1}-2{f}_{i}+{f}_{i-1}}{h^{2}}
  \end{align}
  
where :math:`{f}''_{i}` represents the finite difference approximation to
the second derivative at node :math:`i`. The relations between the
coefficients :math:`a,b,c` and :math:`\alpha,\beta` are derived by matching the
Taylor series coefficients of various orders. The first unmatched coefficient determines the formal truncation
error of the approximation. These constraints are:

.. math::
  \begin{align}
  a+b+c & = 1+2\alpha+2\beta\quad\text{(second order)}\\
  a+2^{2}b+3^{2}c & = \cfrac{4!}{2!}(\alpha+2^{2}\beta)\quad\text{(fourth order)}\\
  a+2^{4}b+3^{4}c & = \cfrac{6!}{4!}(\alpha+2^{4}\beta)\quad\text{(sixth   order)}\\
  a+2^{6}b+3^{6}c & = \cfrac{8!}{6!}(\alpha+2^{6}\beta)\quad\text{(eighth  order)}\\
  a+2^{8}b+3^{8}c & = \cfrac{10!}{8!}(\alpha+2^{8}\beta)\quad\text{(tenth  order)}\\
  \end{align}
  
By choosing :math:`\beta=0` and :math:`c = 0` a one-parameter family of
fourth-order schemes is generated. This family has

.. math::
  \beta=0,\quad c=0,\quad a=\cfrac{4}{3}(1-\alpha),\quad b=\cfrac{1}{3}(-1+10\alpha)

-
  
.. math::
  \begin{array}{c}
  \beta=0,c=0\\
  a+2^{2}b+3^{2}c  = \cfrac{4!}{2!}(\alpha+2^{2}\beta)\quad\text{(fourth order)}\\
  \Rightarrow a+2^{2}b  = \cfrac{4!}{2!}(\alpha)\quad\text{(fourth order)}\\
  \Rightarrow a+4b  = 4\times 3(\alpha)\quad\text{(fourth order)}\\
  \cfrac{4}{3}(1-\alpha)+4\times \cfrac{1}{3}(-1+10\alpha)
  =\cfrac{4}{3}-\cfrac{4}{3}\alpha-\cfrac{4}{3}+\cfrac{4}{3}(10\alpha)=\cfrac{4}{3}(9\alpha)=4\times (3\alpha)
  \end{array}  
  
For :math:`\alpha=\cfrac{1}{10}`  the classical Pade scheme is recovered.

then

.. math::
   \alpha=\cfrac{1}{10},\quad \beta=0,\quad b=0, \quad c=0,\quad a=\cfrac{4}{3}(1-\cfrac{1}{10})=\cfrac{12}{10}

-

.. math::   
  \cfrac{1}{10}{f}''_{i-1}+{f}''_{i}+\cfrac{1}{10}{f}''_{i+1}=\cfrac{12}{10}\left(\cfrac{{f}_{i+1}-2{f}_{i}+{f}_{i-1}}{h^{2}}\right)   


-

.. math::   
  \cfrac{1}{12}{f}''_{i-1}+\cfrac{10}{12}{f}''_{i}+\cfrac{1}{12}{f}''_{i+1}=\cfrac{{f}_{i+1}-2{f}_{i}+{f}_{i-1}}{h^{2}}
  
-

.. math::   
  \cfrac{\partial u}{\partial t}=\alpha\cfrac{\partial ^{2} u}{\partial x^{2}}\Rightarrow 
  \cfrac{\partial ^{2} u}{\partial x^{2}}=\cfrac{1}{\alpha}\cfrac{\partial u}{\partial t}
  
We usually use the order of accuracy as an indicator of the ability of finite difference schemes to approximate the exact solution as it tells us how the discretization error will reduce with mesh refinement. Another way of measuring the order of accuracy is the modified wavenumber approach. In this approach, we see how much the modified wave number is different from the true wave number. The solution of a nonlinear partial differential equation usually contains several frequencies and the modified wavenumber approach provides a way to assess how well the different components of the solution are represented.
The modified wavenumber varies for every finite difference scheme and can be found using Fourier analysis of the differencing scheme.

Compact finite difference schemes have very good resolution characteristics and can be used for capturing high-frequency waves.
In compact formulation, we express the derivative of a function as a linear combination of values of function defined on a set of nodes.
The compact method tries to mimic global dependence and hence has good resolution characteristics. Lele investigated the behavior of compact finite difference schemes for approximating a first and second derivative. The fourth-order accurate approximation to the second derivative is given by

.. math::
  \left.\frac{1}{12} \frac{\partial^{2} u}{\partial x^{2}}\right|_{i-1}+\left.\frac{10}{12} \frac{\partial^{2} u}{\partial x^{2}}\right|_{i}+\left.\frac{1}{12} \frac{\partial^{2} u}{\partial x^{2}}\right|_{i+1}=\frac{u_{i-1}-2 u_{i}+u_{i+1}}{\Delta x^{2}}
  
Taking the linear combination of heat equation and using the same coefficients as the above equation, we get
  
.. math::
  \left.\frac{1}{12} \frac{\partial u}{\partial t}\right|_{i-1}+\left.\frac{10}{12} \frac{\partial u}{\partial t}\right|_{i}+\left.\frac{1}{12} \frac{\partial u}{\partial t}\right|_{i+1}=\alpha\frac{u_{i-1}-2 u_{i}+u_{i+1}}{\Delta x^{2}}  
  
We then use Crank-Nicolson numerical scheme for time discretization and we arrive to below equation

.. math::
  \begin{align}
  &\frac{1}{12} \frac{u^{n+1}_{i-1}-u^{n}_{i-1}}{\Delta t}+\frac{10}{12} \frac{u^{n+1}_{i}-u^{n}_{i}}{\Delta t}+\frac{1}{12} \frac{u^{n+1}_{i+1}-u^{n}_{i+1}}{\Delta t}\\
   = &\frac{\alpha}{2\Delta x^{2}}  (u^{n+1}_{i-1}-2 u^{n+1}_{i}+u^{n+1}_{i+1}+u^{n}_{i-1}-2 u^{n}_{i}+u^{n}_{i+1})
  \end{align}
  
Simplifying above equation, we get the implicit compact Pade (ICP) scheme for heat equation  

.. math::
  a_{i}u^{n+1}_{i-1}+b_{i}u^{n+1}_{i}+c_{i}u^{n+1}_{i+1}=d_{i}
  
where

.. math::
  \begin{array}{l}
  a_{i}=\cfrac{12}{\Delta x^{2}}-\cfrac{2}{\alpha\Delta t}\\
  b_{i}=\cfrac{-24}{\Delta x^{2}}-\cfrac{20}{\alpha\Delta t}\\
  c_{i}=\cfrac{12}{\Delta x^{2}}-\cfrac{2}{\alpha\Delta t}\\
  d_{i}=\cfrac{-2}{\alpha\Delta t}(u^{n+1}_{i-1}+10 u^{n+1}_{i}+u^{n+1}_{i+1})
  -\cfrac{12}{\Delta x^{2}}(u^{n}_{i-1}-2 u^{n}_{i}+u^{n}_{i+1})\\
  \end{array}  
  
proof:

.. math::
  \begin{align}
  &\frac{1}{12\Delta t}u^{n+1}_{i-1}- \frac{1}{12\Delta t}u^{n}_{i-1}+\frac{10}{12\Delta t}u^{n+1}_{i} -\frac{10}{12\Delta t}u^{n}_{i}+\frac{1}{12\Delta t}u^{n+1}_{i+1}-\frac{1}{12\Delta t}u^{n}_{i+1}\\
  = &\frac{\alpha}{2\Delta x^{2}}  (u^{n+1}_{i-1}-2 u^{n+1}_{i}+u^{n+1}_{i+1}+u^{n}_{i-1}-2 u^{n}_{i}+u^{n}_{i+1})
  \end{align}
  
-
  
.. math::
  \begin{align}
  &(\frac{1}{12\Delta t}-\frac{\alpha}{2\Delta x^{2}})u^{n+1}_{i-1}- \frac{1}{12\Delta t}u^{n}_{i-1}\\
  +&(\frac{10}{12\Delta t}+\frac{2\alpha}{2\Delta x^{2}})u^{n+1}_{i} -\frac{10}{12\Delta t}u^{n}_{i}\\
  +&(\frac{1}{12\Delta t}-\frac{\alpha}{2\Delta x^{2}})u^{n+1}_{i+1}-\frac{1}{12\Delta t}u^{n}_{i+1}\\
  =&\frac{\alpha}{2\Delta x^{2}}  (u^{n}_{i-1}-2 u^{n}_{i}+u^{n}_{i+1})
  \end{align} 
  
-
  
.. math::
  \begin{align}
  &(\frac{1}{12\Delta t}-\frac{\alpha}{2\Delta x^{2}})u^{n+1}_{i-1}
  +(\frac{10}{12\Delta t}+\frac{2\alpha}{2\Delta x^{2}})u^{n+1}_{i}
  +(\frac{1}{12\Delta t}-\frac{\alpha}{2\Delta x^{2}})u^{n+1}_{i+1}\\
  =&\frac{\alpha}{2\Delta x^{2}}  (u^{n}_{i-1}-2 u^{n}_{i}+u^{n}_{i+1})
   + \frac{1}{12\Delta t}u^{n}_{i-1}
   + \frac{10}{12\Delta t}u^{n}_{i}
   + \frac{1}{12\Delta t}u^{n}_{i+1}
  \end{align}

-
  
.. math::
  \begin{align}
  &(\frac{1}{12}-\frac{\alpha\Delta t}{2\Delta x^{2}})u^{n+1}_{i-1}
  +(\frac{10}{12}+\frac{2\alpha\Delta t}{2\Delta x^{2}})u^{n+1}_{i}
  +(\frac{1}{12}-\frac{\alpha\Delta t}{2\Delta x^{2}})u^{n+1}_{i+1}\\
  =&\frac{\alpha\Delta t}{2\Delta x^{2}}  (u^{n}_{i-1}-2 u^{n}_{i}+u^{n}_{i+1})
   + \frac{1}{12}u^{n}_{i-1}
   + \frac{10}{12}u^{n}_{i}
   + \frac{1}{12}u^{n}_{i+1}
  \end{align}  
  
Simplifying above equation, we get the implicit compact Pade (ICP) scheme for heat equation  

.. math::
  a_{i}u^{n+1}_{i-1}+b_{i}u^{n+1}_{i}+c_{i}u^{n+1}_{i+1}=d_{i}  
  
where

.. math::
  \begin{array}{l}
  a_{i}=\cfrac{1}{12}-\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  b_{i}=\cfrac{10}{12}+2\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  c_{i}=\cfrac{1}{12}-\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  d_{i}=\cfrac{\alpha\Delta t}{2\Delta x^{2}}  (u^{n}_{i-1}-2 u^{n}_{i}+u^{n}_{i+1})
   + \cfrac{1}{12}u^{n}_{i-1}
   + \cfrac{10}{12}u^{n}_{i}
   + \cfrac{1}{12}u^{n}_{i+1}\\
  \end{array} 

-

.. math::
  \begin{array}{l}
  a_{i}=\cfrac{1}{12}-\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  b_{i}=\cfrac{10}{12}+2\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  c_{i}=\cfrac{1}{12}-\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  d_{i}=(\cfrac{1}{12}+\cfrac{\alpha\Delta t}{2\Delta x^{2}})u^{n}_{i-1}+(\cfrac{10}{12}-2\cfrac{\alpha\Delta t}{2\Delta x^{2}})u^{n}_{i}+(\cfrac{1}{12}+\cfrac{\alpha\Delta t}{2\Delta x^{2}}) u^{n}_{i+1}
  \end{array} 
  
-

.. math::
  \begin{array}{l}
  a_{i}=\cfrac{1}{12}-\cfrac{\alpha\Delta t}{2\Delta x^{2}}\quad \hat{a}_{i}=\cfrac{1}{12}+\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  b_{i}=\cfrac{10}{12}+2\cfrac{\alpha\Delta t}{2\Delta x^{2}}\quad \hat{b}_{i}=\cfrac{10}{12}-2\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  c_{i}=\cfrac{1}{12}-\cfrac{\alpha\Delta t}{2\Delta x^{2}}\quad \hat{c}_{i}=\cfrac{1}{12}+\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  d_{i}=\hat{a}_{i}u^{n}_{i-1}+ \hat{b}_{i}u^{n}_{i}+\hat{c}_{i}u^{n}_{i+1}
  \end{array}    
  
version 2

.. math::
  \begin{array}{l}
  a_{i}=12\cfrac{1}{12}-12\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  b_{i}=12\cfrac{10}{12}+12\times2\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  c_{i}=12\cfrac{1}{12}-12\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  d_{i}=(12\cfrac{1}{12}+12\cfrac{\alpha\Delta t}{2\Delta x^{2}})u^{n}_{i-1}+(12\cfrac{10}{12}-12\times2\cfrac{\alpha\Delta t}{2\Delta x^{2}})u^{n}_{i}+(12\cfrac{1}{12}+12\cfrac{\alpha\Delta t}{2\Delta x^{2}}) u^{n}_{i+1}
  \end{array} 
  
-

.. math::
  \begin{array}{l}
  a_{i}=\cfrac{1}{1}-12\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  b_{i}=\cfrac{10}{1}+12\times2\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  c_{i}=\cfrac{1}{1}-12\cfrac{\alpha\Delta t}{2\Delta x^{2}}\\
  d_{i}=(\cfrac{1}{1}+12\cfrac{\alpha\Delta t}{2\Delta x^{2}})u^{n}_{i-1}+(\cfrac{10}{1}-12\times2\cfrac{\alpha\Delta t}{2\Delta x^{2}})u^{n}_{i}+(\cfrac{1}{1}+12\cfrac{\alpha\Delta t}{2\Delta x^{2}}) u^{n}_{i+1}
  \end{array}  
  
-

.. math::
  \begin{array}{l}
  a_{i}=\cfrac{1}{\alpha\Delta t}-\cfrac{12}{2\Delta x^{2}}\\
  b_{i}=\cfrac{10}{\alpha\Delta t}+24\cfrac{1}{2\Delta x^{2}}\\
  c_{i}=\cfrac{1}{\alpha\Delta t}-\cfrac{12}{2\Delta x^{2}}\\
  d_{i}=(\cfrac{1}{\alpha\Delta t}+\cfrac{12}{2\Delta x^{2}})u^{n}_{i-1}+(\cfrac{10}{\alpha\Delta t}-\cfrac{24}{2\Delta x^{2}})u^{n}_{i}+(\cfrac{1}{\alpha\Delta t}+\cfrac{12}{2\Delta x^{2}}) u^{n}_{i+1}
  \end{array}   
  
-

.. math::
  \begin{array}{l}
  a_{i}=\cfrac{2}{\alpha\Delta t}-\cfrac{12}{\Delta x^{2}}\\
  b_{i}=\cfrac{20}{\alpha\Delta t}+\cfrac{24}{\Delta x^{2}}\\
  c_{i}=\cfrac{2}{\alpha\Delta t}-\cfrac{12}{\Delta x^{2}}\\
  d_{i}=(\cfrac{2}{\alpha\Delta t}+\cfrac{12}{\Delta x^{2}})u^{n}_{i-1}+(\cfrac{20}{\alpha\Delta t}-\cfrac{24}{\Delta x^{2}})u^{n}_{i}+(\cfrac{2}{\alpha\Delta t}+\cfrac{12}{\Delta x^{2}}) u^{n}_{i+1}
  \end{array} 
  
-

.. math::
  \begin{array}{l}
  a_{i}=\cfrac{2}{\alpha\Delta t}-\cfrac{12}{\Delta x^{2}}\\
  b_{i}=\cfrac{20}{\alpha\Delta t}+\cfrac{24}{\Delta x^{2}}\\
  c_{i}=\cfrac{2}{\alpha\Delta t}-\cfrac{12}{\Delta x^{2}}\\
  d_{i}=\cfrac{2}{\alpha\Delta t}(u^{n}_{i-1}+10u^{n}_{i}+u^{n}_{i+1})
  +\cfrac{12}{\Delta x^{2}}(u^{n}_{i-1}-2u^{n}_{i}+u^{n}_{i+1})\\
  \end{array}   
  
-

.. math::
  \begin{array}{l}
  a_{i}=\cfrac{12}{\Delta x^{2}}-\cfrac{2}{\alpha\Delta t}\\
  b_{i}=\cfrac{-24}{\Delta x^{2}}+\cfrac{-20}{\alpha\Delta t}\\
  c_{i}=\cfrac{12}{\Delta x^{2}}-\cfrac{2}{\alpha\Delta t}\\
  d_{i}=\cfrac{-2}{\alpha\Delta t}(u^{n}_{i-1}+10u^{n}_{i}+u^{n}_{i+1})
 +\cfrac{-12}{\Delta x^{2}}(u^{n}_{i-1}-2u^{n}_{i}+u^{n}_{i+1})\\
  \end{array}   