Fourier Series
==================================

#. `The Fast Fourier Transform (FFT): Most Ingenious Algorithm Ever? <https://www.youtube.com/watch?v=h7apO7q16V0/>`_

Wave mathematical representation
----------------------------------
Traveling sinusoidal waves are often represented mathematically in terms of their velocity :math:`v` (in the :math:`x` direction), frequency :math:`f` and wavelength :math:`\lambda` as:

.. math::
  y(x,t)=A\text{cos}(2\pi(\cfrac{x}{\lambda}-ft))=A\text{cos}(\cfrac{2\pi}{\lambda} (x-vt))
  
where :math:`y` is the value of the wave at any position :math:`x` and time :math:`t`, and :math:`A` is the amplitude of the wave.
They are also commonly expressed in terms of wavenumber :math:`k` (:math:`2\pi` times the reciprocal of wavelength) and angular frequency :math:`\omega` (:math:`2\pi` times the frequency) as:

.. math::
  y(x,t)=A\text{cos}(kx-\omega t)=A\text{cos}(k(x-vt))\\
  
in which wavelength and wavenumber are related to velocity and frequency as:

.. math::
  k=\cfrac{2\pi}{\lambda}= \cfrac{2\pi f}{v}=\cfrac{\omega}{v}\\  

-

.. math::
  f=\cfrac{1}{T},T=\cfrac{1}{f}
  
- 
 
.. math::
  v=\cfrac{\lambda }{T}=f\lambda,\lambda=vT,T=\cfrac{\lambda}{v}
  
where :math:`T` is period.

angular frequency

.. math::
  \omega=\cfrac{2\pi}{T}=2\pi f=2\pi \cfrac{v}{\lambda}=kv \\
  
Orthogonal Function System and Series Expansions of Function  
--------------------------------------------------------------------------------
In scientific experiments and engineering, we often face a
periodic motion in which the periodic motion is seen most,
such as spring shock, ball pendulum phenomenon, and
spring rebound phenomenon in high school physics
phenomena, etc., they can be described with a sine function

The Fourier series of the period :math:`T=2L`
-----------------------------------------------------------------------------

.. math::
  \begin{array}{l}
  \displaystyle f(x)=\cfrac{a_{0}}{2}+\sum_{n=1}^{\infty}(a_{n}\text{cos}\cfrac{n\pi x}{L}+b_{n}\text{sin}\cfrac{n\pi x}{L} )\\
  \displaystyle a_{n}=\cfrac{1}{L}\int_{-L}^{L} f(x)\text{cos}\cfrac{n\pi x}{L}dx\quad(n=0,1,2,\dots)\\
  \displaystyle b_{n}=\cfrac{1}{L}\int_{-L}^{L} f(x)\text{sin}\cfrac{n\pi x}{L}dx\quad(n=1,2,3,\dots)
  \end{array}
  
The Fourier series of the period :math:`T=2\pi`
-----------------------------------------------------------------------------

.. math::
  \begin{array}{l}
  \displaystyle f(x)=\cfrac{a_{0}}{2}+\sum_{n=1}^{\infty}(a_{n}\text{cos}\cfrac{n\pi x}{\pi}+b_{n}\text{sin}\cfrac{n\pi x}{\pi} )\\
  \displaystyle a_{n}=\cfrac{1}{\pi}\int_{-\pi}^{\pi} f(x)\text{cos}\cfrac{n\pi x}{\pi}dx\quad(n=0,1,2,\dots)\\
  \displaystyle b_{n}=\cfrac{1}{\pi}\int_{-\pi}^{\pi} f(x)\text{sin}\cfrac{n\pi x}{\pi}dx\quad(n=1,2,3,\dots)
  \end{array} 

-
  
.. math::
  \begin{array}{l}
  \displaystyle f(x)=\cfrac{a_{0}}{2}+\sum_{n=1}^{\infty}(a_{n}\text{cos}(nx)+b_{n}\text{sin}(nx) )\\
  \displaystyle a_{n}=\cfrac{1}{\pi}\int_{-\pi}^{\pi} f(x)\text{cos}(nx) dx\quad(n=0,1,2,\dots)\\
  \displaystyle b_{n}=\cfrac{1}{\pi}\int_{-\pi}^{\pi} f(x)\text{sin}(nx)dx\quad(n=1,2,3,\dots)
  \end{array}
  
Fourier series of odd and even functions :math:`T=2\pi`
---------------------------------------------------------- 
Since the integral of an odd function over a symmetric interval is zero, the integral of an even function over a symmetric interval is twice the integral over a half interval

odd functions:

.. math::
  \begin{array}{l}
  \displaystyle f(x)=\sum_{n=1}^{\infty}(b_{n}\text{sin}(nx) )\\
  \displaystyle a_{n}=0\quad(n=0,1,2,\dots)\\
  \displaystyle b_{n}=\cfrac{2}{\pi}\int_{0}^{\pi} f(x)\text{sin}(nx)dx\quad(n=1,2,3,\dots)
  \end{array}
  
even functions:

.. math::
  \begin{array}{l}
  \displaystyle f(x)=\cfrac{a_{0}}{2}+\sum_{n=1}^{\infty}(a_{n}\text{cos}(nx))\\
  \displaystyle a_{n}=\cfrac{2}{\pi}\int_{0}^{\pi} f(x)\text{cos}(nx) dx\quad(n=0,1,2,\dots)\\
  \displaystyle b_{n}=0\quad(n=1,2,3,\dots)
  \end{array}
  
Even and Odd Extensions
------------------------------
- `Even and Odd Extensions <https://math24.net/even-odd-extensions.html>`_

Suppose that a function :math:`f(x)` is piecewise continuous and defined on the interval :math:`[0,\pi]`. To find its Fourier series, we first extend this function to the interval :math:`[-\pi,\pi]`. This can be done in two ways:  

- We can construct the even extension of :math:`f(x)`:

.. math::
  f_{\text{even}}(x)=\left\{\begin{array}{l}
  f(-x),-\pi\le x\lt0&\\
  f(x)\quad,0\le x\le\pi\\
  \end{array}\right.
  
-  or the odd extension of :math:`f(x)`:

.. math::
  f_{\text{odd}}(x)=\left\{\begin{array}{l}
  -f(-x),-\pi\le x\lt0&\\
  f(x)\quad,0\le x\le\pi\\
  \end{array}\right.
  
For the even function, the Fourier series is called the Fourier Cosine series and is given by

.. math::
  f_{\text{even}}(x)=\cfrac{a_{0}}{2}+\sum_{n=1}^{\infty}a_{n}\text{cos }nx\\

where

.. math::
  a_{n}=\cfrac{2}{\pi}\int_{0}^{\pi} f(x)\text{ cos }nx dx\quad(n=0,1,2,3,\dots)  
  
Respectively, for the odd function, the Fourier series is called the Fourier Sine series and is given by

.. math::
  f_{\text{odd}}(x)=\sum_{n=1}^{\infty}b_{n}\text{sin }nx
  
where the Fourier coefficients are  

.. math::
  b_{n}=\cfrac{2}{\pi}\int_{0}^{\pi} f(x)\text{sin }nxdx\quad(n=1,2,3,\dots)

We can also define the Fourier Sine and Cosine series for a function with an arbitrary period :math:`2L`.
Let :math:`f(x)` be defined on the interval :math:`[0,L]`. Using even extension of the function to interval :math:`[-L,L]`, we obtain
  
.. math::
  f_{\text{even}}(x)=\cfrac{a_{0}}{2}+\sum_{n=1}^{\infty}a_{n}\text{cos }\cfrac{n\pi x}{L}\\

where

.. math::
  a_{n}=\cfrac{2}{L}\int_{0}^{L} f(x)\text{ cos }\cfrac{n\pi x}{L} dx\quad(n=0,1,2,3,\dots)  
  
Respectively, for the odd function, the Fourier series is called the Fourier Sine series and is given by

.. math::
  f_{\text{odd}}(x)=\sum_{n=1}^{\infty}b_{n}\text{sin }\cfrac{n\pi x}{L}
  
where the Fourier coefficients are  

.. math::
  b_{n}=\cfrac{2}{L}\int_{0}^{L} f(x)\text{sin }\cfrac{n\pi x}{L}dx\quad(n=1,2,3,\dots)
  
Period and Angular Frequency
---------------------------------

.. math::  
  \begin{array}{l}
  T=2\pi,\quad\omega=\cfrac{2\pi}{T} =\cfrac{2\pi}{2\pi}=1\\
  T=2L,\quad\omega=\cfrac{2\pi}{T} =\cfrac{2\pi}{L}\\
  \end{array}
  
Complex Form of Fourier Series  
---------------------------------
- `Complex Form of Fourier Series <https://math24.net/complex-form-fourier-series.html>`_

Let the function :math:`f(x)` be defined on the interval :math:`[-\pi,\pi]`.
we can write the Fourier series of the function in complex form:

.. math::
  f(x)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}nx}
  
Here we have used the following notations:

.. math::
  c_{0}=\cfrac{a_{0}}{2},\quad c_{n}=\cfrac{a_{n}-ib_{n}}{2},\quad c_{-n}=\cfrac{a_{n}+ib_{n}}{2}

-

.. math::
  \begin{array}{l}
  \displaystyle a_{n}=\cfrac{1}{\pi}\int_{-\pi}^{\pi} f(x)\text{cos}(nx) dx\quad(n=0,1,2,\dots)\\
  \displaystyle b_{n}=\cfrac{1}{\pi}\int_{-\pi}^{\pi} f(x)\text{sin}(nx)dx\quad(n=1,2,3,\dots)\\
  \end{array}
  
The coefficients :math:`c_{n}` are called complex Fourier coefficients. They are defined by the formulas

.. math::
  c_{n}=\cfrac{1}{2\pi}\int_{-\pi}^{\pi} f(x)e^{-\mathbf{i}nx}dx,\quad(n=0,\pm1,\pm2,\dots)
  
If necessary to expand a function :math:`f(x)` of period :math:`2L`, we can use the following expressions:  

.. math::
  f(x)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{n\pi x}{L}}
  
where

.. math::  
  c_{n}=\cfrac{1}{2L}\int_{-L}^{L} f(x)e^{-\mathbf{i}\cfrac{n\pi x}{L}}dx,\quad(n=0,\pm1,\pm2,\dots)
  
Expand a function :math:`f(x)` of period :math:`T`:  

.. math::
  f(x)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi nx}{T}}
  
where

.. math::  
  c_{n}=\cfrac{1}{T}\int_{-T/2}^{T/2} f(x)e^{-\mathbf{i}\cfrac{2\pi nx}{T}}dx,\quad(n=0,\pm1,\pm2,\dots)  
  
Applications of Fourier Series to Differential Equations
--------------------------------------------------------------
- `Applications of Fourier Series to Differential Equations <https://math24.net/fourier-series-applications-differential-equations.html>`_  

Discrete Fourier Transform (DFT) and its Relation to Fourier Series
-------------------------------------------------------------------------
- `Discrete Fourier Transform (DFT) and its Relation to Fourier Series <https://www.youtube.com/watch?v=6ncRdqf2nYs/>`_  

Fourier Series of a function :math:`f(t)` of period :math:`T`:  

.. math::
  f(t)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi nt}{T}}
  
where

.. math::  
  c_{n}=\cfrac{1}{T}\int_{-T/2}^{T/2} f(t)e^{-\mathbf{i}\cfrac{2\pi nt}{T}}dt,\quad(n=0,\pm1,\pm2,\dots)  
  
Discrete Fourier Transform:

Let :math:`\Delta t=\cfrac{T}{N}`, then
  
.. math::  
   t_{i}=i\Delta t,\quad (i=0,1,2,\dots,N-1)
   
-
   
.. math::  
  \begin{array}{l}
  t_{0}=0\\
  t_{1}=\Delta t\\
  t_{2}=2\Delta t\\
  \dots\\
  t_{i}=i\Delta t\\
  \dots\\
  t_{N-1}=(N-1)\Delta t\\
  \end{array}
  
-
   
.. math::  
  \begin{array}{l}
  \displaystyle f(t)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi nt}{T}}\\
  \displaystyle f_{0}=f(t_{0})=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi nt_{0}}{T}}\\
  \displaystyle f_{1}=f(t_{1})=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi nt_{1}}{T}}\\
  \displaystyle f_{2}=f(t_{2})=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi nt_{2}}{T}}\\
  \dots\\
  \displaystyle f_{i}=f(t_{i})=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi nt_{i}}{T}}\\
  \dots\\
  \displaystyle f_{N-1}=f(t_{N-1})=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi nt_{N-1}}{T}}\\
  \end{array}
  
-
   
.. math:: 
  \begin{array}{l}
  \displaystyle f(t)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi nt}{T}}\\
  \displaystyle f_{0}=f(0)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi n}{T}0}\\
  \displaystyle f_{1}=f(\Delta t)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi n}{T}\Delta t}\\
  \displaystyle f_{2}=f(2\Delta t)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi n}{T}2\Delta t}\\
  \dots\\
  \displaystyle f_{i}=f(i\Delta t)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi n}{T}i\Delta t}\\
  \dots\\
  \displaystyle f_{N-1}=f((N-1)\Delta t)=\sum_{n=-\infty }^{\infty}c_{n}e^{\mathbf{i}\cfrac{2\pi n}{T}(N-1)\Delta t}\\
  \end{array}
  
The coefficients :math:`c_{n}`

.. math::  
  c_{n}=\cfrac{1}{T}\int_{-T/2}^{T/2} f(t)e^{-\mathbf{i}\cfrac{2\pi nt}{T}}dt,\quad(n=0,\pm1,\pm2,\dots)   

-

.. math::  
  c_{n}=\cfrac{1}{T}\int_{0}^{T} f(t)e^{-\mathbf{i}\cfrac{2\pi nt}{T}}dt,\quad(n=0,\pm1,\pm2,\dots)   

-

.. math::  
  c_{n}=\cfrac{1}{T}\int_{0}^{T} f(t)e^{-\mathbf{i}\cfrac{2\pi nt}{T}}dt\approx\cfrac{1}{T}\sum_{i=0}^{N-1}f(i\Delta t)e^{-\mathbf{i}\cfrac{2\pi n}{T}(i\Delta t)}\Delta t,\quad(n=0,\pm1,\pm2,\dots)   
  
-

.. math::  
  c_{n}=\cfrac{1}{N\Delta t}\sum_{i=0}^{N-1}f(i\Delta t)e^{-\mathbf{i}\cfrac{2\pi n}{N\Delta t}(i\Delta t)}\Delta t,\quad(n=0,\pm1,\pm2,\dots)     
  
-

.. math::  
  c_{n}=\cfrac{1}{N}\sum_{i=0}^{N-1}f(i\Delta t)e^{-\mathbf{i}\cfrac{2\pi n}{N}i},\quad(n=0,\pm1,\pm2,\dots)     
  
-

.. math::  
  c_{n}=\cfrac{1}{N}\sum_{i=0}^{N-1}f_{i}e^{-\mathbf{i}\cfrac{2\pi n}{N}i},\quad(n=0,\pm1,\pm2,\dots)     
  
or
  
.. math::  
  c_{n}=\cfrac{1}{N}\sum_{k=0}^{N-1}f_{k}e^{-\mathbf{i}\cfrac{2\pi n}{N}k},\quad(n=0,\pm1,\pm2,\dots)     


Fourier Series
--------------------
#. `Fourier Series Coefficients (In 10 Minutes) <https://www.youtube.com/watch?v=OyNOCs13t-g/>`_

.. math::
  f(x)=a_{0}+\sum_{n=1}^{\infty }\left(a_{n}\text{cos}\left(\cfrac{\pi nx}{L} \right)+b_{n}\text{sin}\left(\cfrac{\pi nx}{L} \right)\right)
  
Note-the argument of the trig functions has been tweaked to account for the fact :math:`f(x)` has period :math:`2L`.

.. math::
  \text{cos}(nx)\longrightarrow \text{cos}(\cfrac{\pi nx}{L})
  
-
  
.. math::  
  \begin{align}
  \int_{-L}^{L}f(x)dx & = \int_{-L}^{L}a_{0}dx+\int_{-L}^{L}\sum_{n = 1}^{\infty }\left(a_{n}\text{cos}\bigg(\cfrac{\pi nx}{L} \bigg)+b_{n}\text{sin}\bigg(\cfrac{\pi nx}{L} \bigg)\right)dx\\
  &=\int_{-L}^{L}a_{0}dx
  +\int_{-L}^{L}\sum_{n = 1}^{\infty }\left(a_{n}\text{cos}\bigg( \cfrac{\pi nx}{L} \bigg)\right)dx
  +\int_{-L}^{L}\sum_{n = 1}^{\infty }\left(b_{n}\text{sin}\bigg(\cfrac{\pi nx}{L} \bigg)\right)dx\\
  &=\int_{-L}^{L}a_{0}dx
  +\sum_{n = 1}^{\infty }a_{n}\int_{-L}^{L}\text{cos}\bigg( \cfrac{\pi nx}{L} \bigg)dx
  +\sum_{n = 1}^{\infty }b_{n}\int_{-L}^{L}\text{sin}\bigg( \cfrac{\pi nx}{L} \bigg)dx
  \end{align}
  
Note: you can't always switch the order of integration and summation with infinite sums and improper integrals. In this case, we can.

.. math::  
  \int_{-L}^{L}\text{cos}\bigg( \cfrac{\pi nx}{L} \bigg)dx=\cfrac{L}{\pi n} \int_{-\pi n}^{\pi n}\text{cos}(u)du
  
-
  
.. math:: 
  \cfrac{L}{\pi n} \int_{-\pi n}^{\pi n}\text{cos}(u)du=\cfrac{L}{\pi n} \big[\text{sin}(u)\big]_{-\pi n}^{\pi n}=0
  
-

.. math::  
  \int_{-L}^{L}\text{sin}\bigg( \cfrac{\pi nx}{L} \bigg)dx=\cfrac{L}{\pi n} \int_{-\pi n}^{\pi n}\text{sin}(u)du
  
-
  
.. math:: 
  \cfrac{L}{\pi n} \int_{-\pi n}^{\pi n}\text{sin}(u)du=-\cfrac{L}{\pi n} \big[\text{cos}(u)\big]_{-\pi n}^{\pi n}=0  
  
-
  
.. math::
  \text{sin}(\pi n)=0=\text{sin}(-\pi n)

-
  
.. math::
  \text{cos}(\pi n)=\text{cos}(-\pi n)
  
so
  
.. math::
  \begin{align}
  \int_{-L}^{L}f(x)dx & = \int_{-L}^{L}(a_{0})dx\\ 
  & = a_{0} \int_{-L}^{L}dx=a_{0}[x]_{-L}^{L}=a_{0}2L
  \end{align}  
  
- 
 
.. math::
  a_{0}=\bigg(\cfrac{1}{2L}\bigg)\int_{-L}^{L}f(x)dx 
  
- 
 
.. math::
  a_{n}=\cfrac{1}{L} \int_{-L}^{L}f(x)\text{cos}\bigg(\cfrac{\pi nx}{L} \bigg)dx    
  
- 
 
.. math::
  b_{n}=\cfrac{1}{L} \int_{-L}^{L}f(x)\text{sin}\bigg(\cfrac{\pi nx}{L} \bigg)dx  
  
Definition of the Fourier Transform and its inverse
-------------------------------------------------------
The so-called Fourier Transform pair is defined this way,

.. math::
  \begin{array}{c}
  \displaystyle \mathcal{F}[f(x)]=\int_{-\infty}^{\infty}f(x)e^{-\mathbf{i}\omega x}dx=F(\omega)\\
  \displaystyle \mathcal{F}^{-1}[F(\omega)]=\cfrac{1}{2\pi} \int_{-\infty}^{\infty}F(\omega)e^{\mathbf{i}\omega x}d\omega=f(x)
  \end{array}
  
or

.. math::
  \begin{array}{c}
  \displaystyle F(\omega)=\int_{-\infty}^{\infty}f(x)e^{-\mathbf{i}\omega x}dx\\
  \displaystyle f(x)=\cfrac{1}{2\pi} \int_{-\infty}^{\infty}F(\omega)e^ {\mathbf{i}\omega x}d\omega\\
  \end{array}
  
Integration by parts
------------------------------
The integration by parts formula states:

.. math::
  \int_{a}^{b}u(x)v'(x)dx =[u(x)v(x)]_{a}^{b}-\int_{a}^{b}u'(x)v(x)dx

Fourier Transforms of derivatives
-------------------------------------------
We shall find the Fourier Transform of :math:`f'(x)`:  

.. math::
  \mathcal{F}[f'(x)]=\int_{-\infty}^{\infty}f'(x)e^{-\mathbf{i}\omega x}dx
  
Let :math:`u(x)=e^{-\mathbf{i}\omega x},v(x)=f(x)`, then

.. math::
  \mathcal{F}[f'(x)]=\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f'(x)dx=[e^{-\mathbf{i}\omega x}f(x)]_{-\infty}^{\infty}
  -(-\mathbf{i}\omega)\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f(x)dx
  
for  

.. math::
  f\to 0\quad \text{as }|x|\to\infty  
  
then

.. math::
  [e^{-\mathbf{i}\omega x}f(x)]_{-\infty}^{\infty}=0
  
-
  
.. math::  
  \mathcal{F}[f'(x)]=(\mathbf{i}\omega)\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f(x)dx=(\mathbf{i}\omega)F(\omega)  
  
The most important aspect of this analysis was the assumption that :math:`f` tends to zero as :math:`x` tends to
:math:`\pm \infty`. This was one of the sufficient conditions that were stated 

Subject to having :math:`f\to 0` as :math:`|x|\to\infty`, an :math:`x`-derivative in the spatial domain is equivalent to
multiplication by :math:`j\omega` in the frequency domain.

The Fourier Transform of a second derivative  
--------------------------------------------------
We shall follow the same idea but will integrate by parts twice.

We shall find the Fourier Transform of :math:`f''(x)`:  

.. math::
  \mathcal{F}[f''(x)]=\int_{-\infty}^{\infty}f''(x)e^{-\mathbf{i}\omega x}dx
  
The integration by parts formula states:

.. math::
  \int_{a}^{b}u(x)v'(x)dx =[u(x)v(x)]_{a}^{b}-\int_{a}^{b}u'(x)v(x)dx  
  
Let :math:`u(x)=e^{-\mathbf{i}\omega x},v(x)=f'(x)`, then

.. math::
  u'(x)=(-\mathbf{i}\omega)e^{-\mathbf{i}\omega x},v'(x)=f''(x)
  
-

.. math::
  \mathcal{F}[f''(x)]=\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f''(x)dx=[e^{-\mathbf{i}\omega x}f'(x)]_{-\infty}^{\infty}
  -(-\mathbf{i}\omega)\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f'(x)dx
  
due to

.. math::
  \mathcal{F}[f'(x)]=\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f'(x)dx=[e^{-\mathbf{i}\omega x}f(x)]_{-\infty}^{\infty}
  -(-\mathbf{i}\omega)\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f(x)dx  
  
and  

.. math::
  f\to 0\quad \text{as }|x|\to\infty    

-
  
.. math::
  f'\to 0\quad \text{as }|x|\to\infty    
  
then 

.. math::
  [e^{-\mathbf{i}\omega x}f(x)]_{-\infty}^{\infty}=0

-

.. math::
  [e^{-\mathbf{i}\omega x}f'(x)]_{-\infty}^{\infty}=0

-

.. math::  
  \mathcal{F}[f'(x)]=\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f'(x)dx=
  -(-\mathbf{i}\omega)\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f(x)dx=(\mathbf{i}\omega)F(\omega)
  
-

.. math::  
  \mathcal{F}[f''(x)]=\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f''(x)dx=
  (\mathbf{i}\omega)\int_{-\infty}^{\infty}e^{-\mathbf{i}\omega x}f'(x)dx=(\mathbf{i}\omega)^{2}F(\omega)
  
Thus we need both :math:`f` and :math:`f'` to tend to zero as :math:`x\to\pm\infty` for this result to be valid. The manner in
which this analysis proceeds tells us that there is a simple rule for yet higher derivatives, namely that  

.. math::  
  \mathcal{F}\bigg[\cfrac{d^{n}f}{dx^{n}}\bigg]=(\mathbf{i}\omega)^{n}F(\omega)
  
provided that :math:`f` and its first :math:`n-1` derivatives tend to zero when :math:`|x|` is large.

Fourier Transform solutions of PDEs
------------------------------------------
Solutions of Fourier’s equation
``````````````````````````````````````  

We will solve Fourier’s equation

.. math::
  \cfrac{\partial T}{\partial t}=\alpha\cfrac{\partial ^{2}T}{\partial x^{2}}   

subject to the initial condition that

.. math::
  T=f(x)\quad \text { at } t=0 

We have to assume that the initial temperature profile, :math:`f(x)`, must decay to zero as :math:`x\to \pm \infty` in
order to be able to use Fourier Transforms. The solution procedure follows three steps:

(i) take the Fourier Transform of the given equation,
(ii) solve the transformed equation,
(iii) take the inverse Fourier Transform to find the desired solution.

Note: In many cases the final answer has to be written in terms of an integral, but sometimes one
may be able to find an explicit expression for the solution.

First we let :math:`\hat{T}(\omega,t)` be the Fourier Transform of :math:`T(x,t)` with respect to :math:`x`. Thus,

.. math::
  \hat{T}(\omega,t)=\mathcal{F}[T(x,t)]=\int_{-\infty}^{\infty}T(x,t)e^{-\mathbf{i}\omega x}dx
  
Now we’ll take the Fourier Transform of the second derivative:

.. math::
  \mathcal{F}\bigg[\cfrac{\partial ^{2}T(x,t)}{\partial x^{2}} \bigg]=\int_{-\infty}^{\infty}\cfrac{\partial ^{2}T(x,t)}{\partial x^{2}}e^{-\mathbf{i}\omega x}dx=-\omega^{2}F(\omega)=-\omega^{2}\hat{T}(\omega)
  
Note that we have used the facts that both :math:`T` and its :math:`x`-derivative tend to zero as :math:`x\to\pm\infty`.
The Fourier Transform of the time-derivative term is relatively straightforward:  

.. math::
  \mathcal{F}\bigg[\cfrac{\partial T(x,t)}{\partial t} \bigg]=\int_{-\infty}^{\infty}\cfrac{\partial T(x,t)}{\partial t}e^{-\mathbf{i}\omega x}dx
  =\cfrac{\partial }{\partial t} \int_{-\infty}^{\infty}T(x,t)e^{-\mathbf{i}\omega x}dx
  =\cfrac{\partial\hat{T}(\omega) }{\partial t}
  
Note that we have effectively swapped the order of differentiation with respect to :math:`t` and integration
with respect to :math:`x` when deriving this result.  

Now Fourier’s equation, is transformed to

.. math::
  \cfrac{\partial\hat{T}(\omega) }{\partial t}=-\alpha\omega^{2}\hat{T}(\omega) 
  
This ordinary differential equation has the solution,  

.. math::
  \hat{T}(\omega)={A}(\omega)e^{-\alpha\omega^{2}t} 
  
where the ‘constant of integration’ may be treated as a function of :math:`\omega`. While this treatment might
seem strange, the ordinary differential equation has :math:`t` as its independent variable, and therefore :math:`\omega` is
simply a passive parameter. Therefore it makes sense to ensure that :math:`A` is as general as possible by
allowing it to vary with :math:`\omega`.  

The constant of integration may be found by applying an appropriate initial condition, but the only
one which is available at present is :

.. math::
  T=f(x)\quad \text { at } t=0 

However, we may take the Fourier Transform of this to obtain

.. math::
  \hat{T}(\omega)=\mathcal{F}[f(x)]=F(\omega)\quad \text { at } t=0
  
Therefore the setting of :math:`t = 0` yields  

.. math::
  A=F
  
and so the Fourier Transform of the desired solution is,

.. math::
  \hat{T}(\omega)={F}(\omega)e^{-\alpha\omega^{2}t} 
  
We may now apply the formula for the inverse Fourier Transform to obtain the final
solution,
  
.. math::
  T(x,t)=\cfrac{1}{2\pi } \int_{-\infty}^{\infty}F(\omega )e^{-\alpha \omega ^{2}t}e^{\mathbf{i}\omega x}d\omega
  
It is difficult to simplify this final integral in any meaningful way, and therefore we have to leave it
as it is. Clearly, for any chosen pair of values of :math:`x` and :math:`t` the integral may be evaluated numerically
with ease and to any desired accuracy. If a sufficient number of pairs of :math:`x` and :math:`t` are used it becomes
possible to obtain either a contour plot of the evolution of :math:`T` in space and time, or else to show how
the temperature profile varies with time.  

Fourier Sine and Cosine Transforms
--------------------------------------
When solving partial differential equations in an infinite domain (i.e. in :math:`-\infty <x<\infty` ), it is often
necessary to use the Fourier Transform. However, not all problems are defined on an infinite domain,
but some are defined on a semi–infinite domain (:math:`0 \le x<\infty`). This where we need to use either
the Fourier Sine Transform (FST) or the Fourier Cosine Transform (FCT).

The FCT and the FST are intimately related to the Fourier Transform and they and their inverses
may be derived from it and its inverse. The Fourier Cosine Transform pair is given by

.. math::
  \begin{array}{c}
  \displaystyle F_{c}(\omega)=\mathcal{F}_{c}[f(x)]=\int_{0}^{\infty}f(x)\text{cos}(\omega x)dx\\
  \displaystyle f(x)=\mathcal{F}_{c}^{-1}[F_{c}(\omega)]=\cfrac{2}{\pi} \int_{-\infty}^{\infty}F_{c}(\omega)\text{cos}(\omega x)d\omega
  \end{array}
  
The Fourier Sine Transform pair is given by  

.. math::
  \begin{array}{c}
  \displaystyle F_{s}(\omega)=\mathcal{F}_{s}[f(x)]=\int_{0}^{\infty}f(x)\text{sin}(\omega x)dx\\
  \displaystyle f(x)=\mathcal{F}_{s}^{-1}[F_{s}(\omega)]=\cfrac{2}{\pi} \int_{-\infty}^{\infty}F_{s}(\omega)\text{sin}(\omega x)d\omega
  \end{array}
  
If we are solving an equation for :math:`u(x, y)`, say, where transforms are being taken in the :math:`x`-direction,
then the FST is used when :math:`u(x, y)` is given on the boundary :math:`x = 0`, and that the FCT is used when
the first :math:`x`-derivative of :math:`u(x, y)` is given on :math:`x = 0`. The reasons for this are technical, and they arise
naturally during the integration by parts process for evaluating the transforms of derivatives.  

Some Fourier Sine Transform examples.
-----------------------------------------
Example 1
``````````````````````````
We will consider the following unsteady one-dimensional heat transfer problem. A semi-infinite solid,
which occupies the region, :math:`0 \le x<\infty`, has the temperature profile, :math:`T = f(x)`, at :math:`t = 0`. However,
the :math:`x = 0` end of this region is maintained at the temperature, :math:`T = 0`. Determine the evolution of
the temperature profile

Let :math:`\hat{T}_{s}(\omega,t)` be the Fourier Sine Transform of :math:`T(x, t)` with respect to :math:`x`, i.e. that

.. math::
  \hat{T}_{s}(\omega,t) =\mathcal{F}_{s}[T(x,t)]=\int_{0}^{\infty}T(x,t)\text{sin}(\omega x)dx

Let :math:`\hat{T}_{c}(\omega,t)` be the Fourier Cosine Transform of :math:`T(x, t)` with respect to :math:`x`, i.e. that
  
.. math::
  \hat{T}_{c}(\omega,t) =\mathcal{F}_{c}[T(x,t)]=\int_{0}^{\infty}T(x,t)\text{cos}(\omega x)dx  
  
We will solve Fourier’s equation

.. math::
  \cfrac{\partial T}{\partial t}=\alpha\cfrac{\partial ^{2}T}{\partial x^{2}}   

subject to the initial condition that

.. math::
  T=f(x)\quad \text { at } t=0 
  
Therefore we will need to take the Fourier Sine Transform of this equation.
Beginning with the time-derivative term, we get  

.. math::
  \begin{array}{l}
  \displaystyle \mathcal{F}_{s}[T(x,t)]=\int_{0}^{\infty}T(x,t)\text{sin}(\omega x)dx=\hat{T}_{s}(\omega)\\
  \displaystyle \mathcal{F}_{s}\bigg[\cfrac{\partial T(x,t)}{\partial t} \bigg]=\int_{0}^{\infty}\cfrac{\partial T(x,t)}{\partial t} \text{sin}(\omega x)dx
  =\cfrac{\partial }{\partial t}\int_{0}^{\infty}T(x,t)\text{sin}(\omega x)dx
  =\cfrac{\partial \hat{T}_{s}(\omega) }{\partial t}
  \end{array}

-
  
.. math::
  \mathcal{F}_{s}\bigg[\cfrac{\partial T(x,t)}{\partial x} \bigg]=\int_{0}^{\infty}\cfrac{\partial T(x,t)}{\partial x} \text{sin}(\omega x)dx  


The integration by parts formula states:

.. math::
  \int_{a}^{b}u(x)v'(x)dx =[u(x)v(x)]_{a}^{b}-\int_{a}^{b}u'(x)v(x)dx
  
Let :math:`u(x)=\text{sin}(\omega x),v(x)= T(x,t)`, then  

.. math::
  u'(x)=\omega\text{cos}(\omega x)
  
-
  
.. math:: 
  \int_{0}^{\infty}\cfrac{\partial T(x,t)}{\partial x} \text{sin}(\omega x)dx
  =\bigg[\text{sin}(\omega x)T(x,t)\bigg]_{0}^{\infty}
  -\omega\int_{0}^{\infty} T(x,t)\text{cos}(\omega x)dx 

-

.. math:: 
  \begin{array}{l}
  \text{sin}(\omega 0)=\text{sin}(0)=0\\
  T(x,t)\to 0\quad \text{as }x\to\infty\\
  \bigg[\text{sin}(\omega x)T(x,t)\bigg]_{0}^{\infty}=0
  \end{array}
  
-

.. math:: 
  \mathcal{F}_{s}\bigg[\cfrac{\partial T(x,t)}{\partial x} \bigg]=\int_{0}^{\infty}\cfrac{\partial T(x,t)}{\partial x} \text{sin}(\omega x)dx 
  =-\omega\int_{0}^{\infty} T(x,t)\text{cos}(\omega x)dx
  
we will need to take the Fourier Cosine Transform of this equation.
Beginning with the time-derivative term, we get  

.. math::
  \begin{array}{l}
  \displaystyle \mathcal{F}_{c}[T(x,t)]=\int_{0}^{\infty}T(x,t)\text{cos}(\omega x)dx=\hat{T}_{c}(\omega)\\
  \displaystyle \mathcal{F}_{c}\bigg[\cfrac{\partial T(x,t)}{\partial x} \bigg]=\int_{0}^{\infty}\cfrac{\partial T(x,t)}{\partial x} \text{cos}(\omega x)dx
  \end{array}
  
Let :math:`u(x)=\text{cos}(\omega x),v(x)= T(x,t)`, then  

.. math::
  u'(x)=-\omega\text{sin}(\omega x)  

-
  
.. math::  
  \begin{align}
  \displaystyle\int_{0}^{\infty}\cfrac{\partial T(x,t)}{\partial x} \text{cos}(\omega x)dx
  & = \bigg[\text{cos}(\omega x)T(x,t)\bigg]_{0}^{\infty}
  +\omega\int_{0}^{\infty} T(x,t)\text{sin}(\omega x)dx \\ 
  & = \bigg[\text{cos}(\omega x)T(x,t)\bigg]_{0}^{\infty}
  +\omega \hat{T}_{s}(\omega)
  \end{align} 
  
The Sine Transform of a second derivative

.. math:: 
  \mathcal{F}_{s}\bigg[\cfrac{\partial ^{2}T(x,t)}{\partial x^{2}} \bigg]=\int_{0}^{\infty}\cfrac{\partial ^{2}T(x,t)}{\partial x^{2}} \text{sin}(\omega x)dx\\
  
Let :math:`u(x)=\text{sin}(\omega x),v(x)= \cfrac{\partial T(x,t)}{\partial x}`, then  

.. math::
  u'(x)=\omega\text{cos}(\omega x)  

-
  
.. math::
  \int_{0}^{\infty}\cfrac{\partial T(x,t)}{\partial x} \text{cos}(\omega x)dx  
  
-
  
.. math::
  \begin{align}
  \int_{0}^{\infty}\cfrac{\partial ^{2}T(x,t)}{\partial x^{2}} \text{sin}(\omega x)dx & = \bigg[\text{sin}(\omega x)\cfrac{\partial T(x,t)}{\partial x} \bigg]_{0}^{\infty}
  -\omega\int_{0}^{\infty}\cfrac{\partial T(x,t)}{\partial x} \text{cos}(\omega x)dx\\
  &=\bigg[\text{sin}(\omega x)\cfrac{\partial T(x,t)}{\partial x} \bigg]_{0}^{\infty}
  -\omega\bigg\{\bigg[\text{cos}(\omega x)T(x,t)\bigg]_{0}^{\infty}
    +\omega \hat{T}_{s}(\omega)\bigg\}\\
  &=\bigg[\text{sin}(\omega x)\cfrac{\partial T(x,t)}{\partial x} \bigg]_{0}^{\infty}
  -\omega\bigg[\text{cos}(\omega x)T(x,t)\bigg]_{0}^{\infty}-\omega^{2} \hat{T}_{s}(\omega)
  \end{align}  
  
Note that we had to use the :math:`T = 0 \text{ at }x = 0` boundary condition.
Note also that, if the boundary condition at :math:`x = 0` had been
a Neumann condition (the gradient of :math:`T` is specified) then we would not be able to proceed further
because the value of :math:`T` at :math:`x = 0` is needed there. Thus the Fourier Sine Transform must be used
with Dirichlet boundary conditions.

.. math::
  \mathcal{F}_{s}\bigg[\cfrac{\partial ^{2}T(x,t)}{\partial x^{2}}\bigg] =\int_{0}^{\infty}\cfrac{\partial ^{2}T(x,t)}{\partial x^{2}} \text{sin}(\omega x)dx 
  =-\omega^{2} \hat{T}_{s}(\omega)
  
We also assumed that :math:`T(x,t)` and its derivatives decay to zero as :math:`x\to\infty`. Therefore Fourier’s equation
transforms to,  

.. math::
  \cfrac{\partial \hat{T}_{s}(\omega)}{\partial t}=-\alpha\omega^{2} \hat{T}_{s}(\omega)
  
and the solution is,  

.. math::
  \hat{T}_{s}(\omega)=A(\omega)e^{-\alpha\omega^{2}t}
  
where A(ω) is the arbitrary constant. The given initial condition is
that :math:`T = f(x)` when :math:`t = 0`; on taking the Fourier Sine Transform of this, we obtain the transformed
version, namely,  

.. math::
  \hat{T}_{s}(\omega)=\mathcal{F}_{s}[f(x)]={F}_{s}(\omega)\quad\text{when }t=0  
  
Substitution yields :math:`A = {F}_s`, and hence  

.. math::
  \hat{T}_{s}(\omega)={F}_s(\omega)e^{-\alpha\omega^{2}t}
  
The final solution is obtained by taking the inverse Fourier Sine Transform of this; we get, 
 
.. math::
  T(x,t)=\mathcal{F}_{s}^{-1}[\hat{T}_{s}(\omega)]=\cfrac{2}{\pi}\int_{0}^{\infty}F_{s}(\omega) e^{-\alpha \omega ^{2}t}\text{sin}(\omega x)d\omega
  
The Discrete Fourier Transform
---------------------------------------
The Discrete Fourier Transform (DFT) is the equivalent of the continuous Fourier
Transform for signals known only at :math:`N` instants separated by sample times :math:`T` (i.e.
a finite sequence of data).

Let :math:`f(t)` be the continuous signal which is the source of the data. Let :math:`N` samples
be denoted :math:`f[0],f[1],f[2],\dots,f[k],\dots,f[N-1]`.

The Fourier Transform of the original signal, :math:`f(t)`, would be

.. math::
  F(\omega)=\int_{-\infty}^{\infty}f(t)e^{\mathbf{i}\omega t}dt

We could regard each sample :math:`f[k]` as an impulse having area :math:`f[k]`. Then, since the integrand exists only at the sample points:

.. math::
  F(\omega)=\int_{-\infty}^{\infty}f(t)e^{\mathbf{i}\omega t}dt
  
-
  
.. math::  
  \begin{align}
  F(\omega) & = \int_{-\infty}^{\infty}f(t)e^{-\mathbf{i}\omega t}dt  = \int_{0}^{(N-1)T}f(t)e^{-\mathbf{i}\omega t}dt\\
  & = f[0]e^{-\mathbf{i}0}+f[1]e^{-\mathbf{i}\omega T}+f[2]e^{-\mathbf{i}\omega 2T}
  +\dots+f[k]e^{-\mathbf{i}\omega kT}+\dots+f[N-1]e^{-\mathbf{i}\omega (N-1)T}
  \end{align}  
  
-
  
.. math::   
  F(\omega) =\sum_{k=0}^{N-1}f[k]e^{-\mathbf{i}\omega kT}
  
We could in principle evaluate this for any :math:`\omega`, but with only :math:`N` data points to start
with, only :math:`N` final outputs will be significant.

You may remember that the continuous Fourier transform could be evaluated
over a finite interval (usually the fundamental period :math:`T_{0}` ) rather than from :math:`-\infty` to
:math:`+\infty` if the waveform was periodic. Similarly, since there are only a finite number
of input data points, the DFT treats the data as if it were periodic (i.e. :math:`f(N)` to :math:`f(2N-1)` is the same as
:math:`f(0)` to :math:`f(N-1)` ).

Since the operation treats the data as if it were periodic, we evaluate the
DFT equation for the fundamental frequency (one cycle per sequence, :math:`\cfrac{1}{NT}` Hz, :math:`\cfrac{2\pi}{NT}` rad/sec.)
and its harmonics (not forgetting the d.c. component (or average) at :math:`\omega=0`).
i.e. set

.. math::  
  \omega=0,\cfrac{2\pi}{NT},\cfrac{2\pi\cdot2}{NT},\dots,\cfrac{2\pi\cdot n}{NT},\dots,\cfrac{2\pi\cdot (N-1)}{NT}
  
or, in general

.. math::  
  \begin{array}{l}
  \displaystyle F(\omega) =\sum_{k=0}^{N-1}f[k]e^{-\mathbf{i}\omega kT}\\
  \displaystyle F(\cfrac{2\pi\cdot n}{NT}) =\sum_{k=0}^{N-1}f[k]e^{-\mathbf{i}\cfrac{2\pi\cdot n}{NT} kT}
  =\sum_{k=0}^{N-1}f[k]e^{-\mathbf{i}\cfrac{2\pi}{N} nk}\\
  \displaystyle \tilde F(n) =\sum_{k=0}^{N-1}f[k]e^{-\mathbf{i}\cfrac{2\pi}{N} nk}\\
  \end{array}