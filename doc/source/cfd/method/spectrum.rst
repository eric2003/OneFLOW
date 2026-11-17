Spectral Methods 
==================================

#. `Scientific Computing || 02 Week 7 19 1 Introduction to spectral methods 10 46 <https://www.youtube.com/watch?v=ymsY8IFbOwY/>`_
#. `Introduction to Spectral Methods for Partial Differential Equations <https://www.youtube.com/watch?v=rsdqvrldgHM/>`_
#. `Fourier Series <https://mathworld.wolfram.com/FourierSeries.html>`_
#. `Spectral Numerical Method <https://www.youtube.com/watch?v=wbk6QaTzqOY/>`_


Examples of spectral methods
---------------------------------
Here we presume an understanding of basic multivariate calculus and Fourier series. If :math:`g(x,y)` is a known, complex-valued
function of two real variables, and :math:`g` is periodic in x and y ( that is, :math:`g(x,y)=g(x+2\pi,y)=g(x,y+2\pi)` ) then wei are
interested in  finding a function :math:`f(x,y)` so that

.. math::
  \left(\cfrac{\partial }{\partial x^{2}}+\cfrac{\partial }{\partial y^{2}} \right)f(x,y)=g(x,y)\quad \text{for all }x,y
  
where the expression on the left denotes the second partial derivatives of :math:`f` in :math:`x` and :math:`y` respectively.
This is the Poisson equation, and can be physically interpreted as some sort of heat conduction problem, or a problem in potential theory, among other possibilities.
If we write :math:`f` and :math:`g` in Fourier series:

.. math::
  \begin{array}{c}
  f=\sum a_{j,k}e^{\mathbf{i}(jx+ky)}\\
  g=\sum b_{j,k}e^{\mathbf{i}(jx+ky)}
  \end{array}
  
and substitute into the differential equation, we obtain this equation:

.. math::
  \begin{array}{l}
  \cfrac{\partial f}{\partial x}=\sum j\mathbf{i}a_{j,k}e^{\mathbf{i}(jx+ky)}\\
  \cfrac{\partial^{2} f}{\partial x^{2}}=-\sum j^{2}a_{j,k}e^{\mathbf{i}(jx+ky)}\\
  \cfrac{\partial f}{\partial y}=\sum k\mathbf{i}a_{j,k}e^{\mathbf{i}(jx+ky)}\\
  \cfrac{\partial ^{2}f}{\partial y^{2}}=-\sum k^{2}a_{j,k}e^{\mathbf{i}(jx+ky)}\\
  \end{array}

-

.. math:: 
  \begin{array}{l}
  -\sum (j^{2}+k^{2})a_{j,k}e^{\mathbf{i}(jx+ky)}=\sum b_{j,k}e^{\mathbf{i}(jx+ky)}\\
  -(j^{2}+k^{2})a_{j,k}=b_{j,k}
  \end{array}  
  
We have exchanged partial differentiation with an infinite sum, which is legitimate if we assume for instance that :math:`f` has a continuous second derivative. By the uniqueness theorem for Fourier expansions, we must then equate the Fourier coefficients term by term, giving  

.. math:: 
  \begin{array}{l}
  a_{j,k}=-\cfrac{b_{j,k}}{j^{2}+k^{2}}
  \end{array}  
  
Example 1 of spectral methods
---------------------------------  
Consider the following PDE with constant :math:`a`

.. math::
   \cfrac{\partial u}{\partial t}=a\cfrac{\partial u}{\partial x}

Represent the solution :math:`u(x,t)` as

.. math::
  u(x,t)=\sum_{n}\hat{u}_{n}(t)\phi_{n}(x)
  
where :math:`\phi_{n}` are some basis functions and :math:`\hat{u}_{n}` are corresponding weighting factors.

In a spectral method we take :math:`\phi_{n}` as :math:`e^{\mathbf{i}2\pi n\cfrac{x}{L}}=\text{cos}(2\pi n\cfrac{x}{L})+\mathbf{i}\text{sin}(2\pi n\cfrac{x}{L})`,
where :math:`\mathbf{i}` is the imaginary number and :math:`L` is the domian size.

- Here, :math:`\phi_{n}` are periodic functions and we assume periodic boundary conditions.
- This gives :math:`u` as a discrete inverse Fourier transform, and we take a(truncated) summation from :math:`n=-N/2+1` to :math:`n=N/2`.

This results in 

.. math::
  u(x,t)=\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}(t)e^{\mathbf{i}{2\pi}\cfrac{nx}{L}}
  
In particular, consider a grid of :math:`N` points with uniform spacing :math:`\Delta x`. If the points are indexed by :math:`j`, starting at :math:`j=0`,
we have :math:`\Delta x=\cfrac{L}{N}`, and :math:`x_{j}=j\Delta x = j\cfrac{L}{N}`. Let :math:`u_{j}=u(x_{j})`,

Then the inverse discrete Fourier transform(IDFT) of :math:`\hat{u}_{n}` evaluated at grid points :math:`x_{j}` and denoted 
:math:`F^{-1}`, is given by 

.. math::
  \begin{array}{l}
  \displaystyle u(x_{j},t)=\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}(t)e^{\mathbf{i}{2\pi}\cfrac{nx_{j}}{L}}\\
  \displaystyle u(x_{j},t)=u_{j}(t)=\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}(t)e^{\mathbf{i}{2\pi}\cfrac{nj\Delta x}{L}}\\
  \displaystyle u(x_{j},t)=u_{j}(t)=\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}(t)e^{\mathbf{i}{2\pi}\cfrac{nj}{L}\cfrac{L}{N}}\\
  \displaystyle u(x_{j},t)=u_{j}(t)=\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}(t)e^{\mathbf{i}{2\pi}\cfrac{nj}{N}}\\
  \end{array}
  
-
  
.. math::
  u(x_{j},t)=u_{j}(t)=F^{-1}(\hat{u}_{n}(t))=\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}(t)e^{\mathbf{i}{2\pi}\cfrac{nj}{N}},j=0,\cdots,n-1.\\  
  
or   

.. math::
  u_{j}=F^{-1}(\hat{u}_{n})=\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}e^{\mathbf{i}{2\pi}\cfrac{nj}{N}},j=0,\cdots,n-1.\\
  
The corresponding discrete Fourier transform(DFT) is given by

.. math::
  \hat{u}_{n}=F(u_{j})=\cfrac{1}{N}\sum_{j=0}^{N-1}u_{j}e^{-\mathbf{i}{2\pi}\cfrac{nj}{N}},n=-\cfrac{N}{2}+1,\cdots,\cfrac{N}{2}
  
Consider the PDE  

.. math::
  \begin{array}{l}
  \displaystyle u(x,t)=\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}(t)e^{\mathbf{i}{2\pi}\cfrac{nx}{L}}\\
  \displaystyle \cfrac{\partial u(x,t)}{\partial t}=\sum_{n=-N/2+1}^{N/2}\cfrac{\partial\hat{u}_{n}(t)}{\partial t} e^{\mathbf{i}{2\pi}\cfrac{nx}{L}}\\
  \displaystyle \cfrac{\partial u(x,t)}{\partial x}={\mathbf{i}{2\pi}\cfrac{n}{L}}\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}(t) e^{\mathbf{i}{2\pi}\cfrac{nx}{L}}\\
  \end{array}  
  
then

.. math::
  \displaystyle \sum_{n=-N/2+1}^{N/2}\cfrac{\partial\hat{u}_{n}(t)}{\partial t} e^{\mathbf{i}{2\pi}\cfrac{nx}{L}}
  =a{\mathbf{i}{2\pi}\cfrac{n}{L}}\sum_{n=-N/2+1}^{N/2}\hat{u}_{n}(t) e^{\mathbf{i}{2\pi}\cfrac{nx}{L}}

-

.. math::
  \cfrac{\partial\hat{u}_{n}(t)}{\partial t}=a{\mathbf{i}{2\pi}\cfrac{n}{L}}\hat{u}_{n}
  
Given some initial condition for :math:`u_{j}`, we can use the DFT to compute the initial :math:`\hat{u}_{n}`. This ODE can ten be
solved for :math:`\hat{u}_{n}(t)`. Then :math:`{u}_{j}(t)` can be computed form :math:`\hat{u}_{n}(t)` using the IDFT.
Note, for this particular problem, the ODE has a simple analytic solution:

.. math::
  \hat{u}_{n}(t)=\hat{u}_{n}(0)e^{\left(a\mathbf{i}{2\pi}\cfrac{n}{L}\right)t}
  
Example 2 of spectral methods  
-----------------------------------
Consider the following PDE:

.. math::
  \cfrac{\partial ^{2}u}{\partial x^{2}} +\cfrac{\partial ^{2}u}{\partial y^{2}}=f(x,y)
  
Represent the solution :math:`u(x,y)` as  

.. math::
  u(x,y)=\sum_{m,n}\hat{u}_{m,n}\phi_{m}(x)\phi_{n}(y)
  
where :math:`\phi_{m}` and :math:`\phi_{n}` are some basis functions and :math:`\hat{u}_{m,n}` are corresponding weighting factors.

In a spectral method we take

.. math::
  \phi_{m}(x)=e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}=\text{cos}(2\pi m\cfrac{x}{L_{x}})+\mathbf{i}\text{sin}(2\pi m\cfrac{x}{L_{x}})
  
-  

.. math::
  \phi_{n}(y)=e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}=\text{cos}(2\pi n\cfrac{y}{L_{y}})+\mathbf{i}\text{sin}(2\pi n\cfrac{y}{L_{y}})
  
where :math:`\mathbf{i}` is the imaginary number. :math:`L_{x}` and :math:`L_{y}` are :math:`x` and :math:`y` direction Length.
then

.. math::
  u(x,y)=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}
  
-
  
.. math::  
  u(x,y)=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi(m\cfrac{x}{L_{x}}+n\cfrac{y}{L_{y}}) }  
  
Let

.. math::  
  \begin{array}{c}
  \Delta x = \cfrac{L_{x}}{M}， x_{i}=i\Delta x= i\cfrac{L_{x}}{M} \\
  \Delta y = \cfrac{L_{y}}{N}， y_{j}=j\Delta y= j\cfrac{L_{y}}{N} \\
  \end{array}
  
Then the inverse discrete Fourier transform(IDFT):

.. math::  
  \begin{array}{l}
  \displaystyle u(x,y)=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi(m\cfrac{x}{L_{x}}+n\cfrac{y}{L_{y}}) }\\
  \displaystyle u(x_{i},y_{j})=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi(m\cfrac{x_{i}}{L_{x}}+n\cfrac{y_{j}}{L_{y}}) }\\
  \displaystyle u(x_{i},y_{j})=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi(m\cfrac{i\Delta x}{L_{x}}+n\cfrac{j\Delta y}{L_{y}}) }\\
  \displaystyle u(x_{i},y_{j})=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi(m\cfrac{iL_{x}/M}{L_{x}}+n\cfrac{jL_{y}/N}{L_{y}}) }\\
  \displaystyle u(x_{i},y_{j})=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi(\cfrac{im}{M}+\cfrac{jn}{N}) }\\
  \end{array}
  
or

.. math::  
  u_{i,j}=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi(\cfrac{im}{M}+\cfrac{jn}{N}) }\\  
  
and   

.. math::    
  i=0,1,\cdots,M-1; \quad j=0,1,\cdots ,N-1
  
The corresponding discrete Fourier transform(DFT) is given by  

.. math::
  \hat{u}_{m,n}=\cfrac{1}{MN} \sum_{i=0}^{M-1}\sum_{j=0}^{N-1}u_{i,j}e^{-\mathbf{i}2\pi(\cfrac{im}{M}+\cfrac{jn}{N}) }\\
  
where

.. math::
  \begin{array}{l}
  m=-\cfrac{M}{2}+1,\cdots,\cfrac{M}{2}\\
  n=-\cfrac{N}{2}+1,\cdots,\cfrac{N}{2}\\
  \end{array}
  
Consider the PDE:  

.. math::
  \begin{array}{l}
  \displaystyle u(x,y)=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}\\
  \displaystyle \cfrac{\partial u(x,y)}{\partial x} ={\mathbf{i}2\pi m\cfrac{1}{L_{x}}}\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}\\
  \displaystyle \cfrac{\partial ^{2}u(x,y)}{\partial x^{2}} =\left\{\mathbf{i}2\pi m\cfrac{1}{L_{x}}\right\}^{2}\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}\\
  \displaystyle \cfrac{\partial u(x,y)}{\partial y} ={\mathbf{i}2\pi n\cfrac{1}{L_{y}}}\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}\\
  \displaystyle \cfrac{\partial ^{2}u(x,y)}{\partial y^{2}} =\left\{\mathbf{i}2\pi n\cfrac{1}{L_{y}}\right\}^{2}\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{u}_{m,n}e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}\\
  \end{array}
  
and 

.. math::
  f(x,y)=\sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}\hat{f}_{m,n}e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}\\
  
then

.. math::
  \begin{array}{l}
  \cfrac{\partial ^{2}u}{\partial x^{2}} +\cfrac{\partial ^{2}u}{\partial y^{2}}=f(x,y)\\
  \displaystyle \Rightarrow \sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}\hat{u}_{m,n}\left\{\left(\mathbf{i}2\pi m\cfrac{1}{L_{x}}\right )^2+\left(\mathbf{i}2\pi n\cfrac{1}{L_{y}}\right )^2\right\}\\
  \displaystyle = \sum_{m=-M/2+1}^{M/2}\sum_{n=-N/2+1}^{N/2}e^{\mathbf{i}2\pi m\cfrac{x}{L_{x}}}e^{\mathbf{i}2\pi n\cfrac{y}{L_{y}}}\hat{f}_{m,n}
  \end{array}    

-
  
.. math::
  \hat{u}_{m,n}\left\{\left(\mathbf{i}2\pi m\cfrac{1}{L_{x}}\right )^2+\left(\mathbf{i}2\pi n\cfrac{1}{L_{y}}\right )^2\right\}=\hat{f}_{m,n}  
  
-
  
.. math::
  \begin{array}{l}
  -\hat{u}_{m,n}\left\{\left(\cfrac{2\pi m}{L_{x}}\right )^2+\left(\cfrac{2\pi n}{L_{y}}\right )^2\right\}=\hat{f}_{m,n}\\
  \hat{u}_{m,n}
  =-\cfrac{\hat{f}_{m,n}}{\left(\cfrac{2\pi m}{L_{x}}\right )^2+\left(\cfrac{2\pi n}{L_{y}}\right )^2}\\
  \end{array}  