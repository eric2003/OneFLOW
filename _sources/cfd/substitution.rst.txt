Substitution
==================================


Substitution for Integrals
---------------------------------------------------------------------------
- `Substitution for Definite Integrals <https://web.ma.utexas.edu/users/m408n/m408c/CurrentWeb/LM5-5-7.php/>`_
- `Substitution for Definite Integrals-math.libretexts.org <https://math.libretexts.org/Courses/Mount_Royal_University/MATH_1200%3A_Calculus_for_Scientists_I/4%3A_Integral_Calculus/4.7%3A_Definite_integrals_by_substitution.>`_
- `Lesson Explainer: Integration by Substitution: Definite Integrals <https://www.nagwa.com/en/explainers/816185750769/>`_
- `Substitution Rule for Indefinite Integrals <https://www.sfu.ca/math-coursenotes/Math%20158%20Course%20Notes/sec_SubRule.html>`_


Substitution Rule for Indefinite Integrals
---------------------------------------------------------------------------

If :math:`u=g(x)` is a differential function whose range is an interval :math:`I`, and :math:`f` is continuous on :math:`I`, then

.. math::
  \int f(g(x))\cdot {g}'(x)dx = \int f(u)du
  
Proof

By the Chain Rule, :math:`F(g(x))` is an antiderivative of :math:`f(g(x))\cdot {g}'(x)` whenever :math:`F` is an antiderivative of :math:`f`, beacuse

.. math::
  \begin{align}
  \cfrac{d}{dx}F(g(x)) & = {F}'(g(x))\cdot {g}'(x)\\ & = f(g(x)){g}'(x)
  \end{align}
  
If we make the substitution :math:`u=g(x)`, then

.. math::
  \begin{align}
  \int f(g(x))\cdot {g}'(x)dx &= \int \cfrac{d}{dx}F(g(x))dx\\
  &=F(g(x))+C\\
  &=F(u)+C\\
  &=\int {F}'(u)du\\
  &=\int f(u)du\\
  \end{align}

Substitution In Definite Integrals
---------------------------------------------------------------------------
If :math:`{g}'` is continuous on the interval :math:`[a,b]` and :math:`f` is continuous on the range of 
:math:`g(x)=u`, then

.. math::
  \int_{a}^{b}f(g(x))\cdot {g}'(x)dx=\int_{g(a)}^{g(b)}f(u)du\\
  
Proof 

Let :math:`F` denote any antiderivative of :math:`f`. Then,  

.. math::
  \begin{aligned}
  \int_{a}^{b} f(g(x)) \cdot g^{\prime}(x) d x & =F(g(x))\Bigg ]_{x=a}^{x=b} \\
  & =F(g(b))-F(g(a)) \\
  & =F(u)\Bigg]_{u=g(a)}^{u=g(b)} \\
  & =\int_{g(a)}^{g(b)} f(u) d u .
  \end{aligned}
  
Substitution Rule for Indefinite Integrals( version 2)
---------------------------------------------------------
  
For convenience, swap :math:`x` and :math:`u`  

If :math:`x=g(u)` is a differential function whose range is an interval :math:`I`, and :math:`f` is continuous on :math:`I`, then

.. math::
  \int f(x)dx = \int f(g(u))\cdot {g}'(u)du
  
Proof

By the Chain Rule, :math:`F(g(u))` is an antiderivative of :math:`f(g(u))\cdot {g}'(u)` whenever :math:`F` is an antiderivative of :math:`f`, beacuse

.. math::
  \begin{align}
  \cfrac{d}{du}F(g(u)) & = {F}'(g(u))\cdot {g}'(u)\\ & = f(g(u)){g}'(u)
  \end{align}
  
If we make the substitution :math:`x=g(u)`, then

.. math::
  \begin{align}
  \int f(g(u))\cdot {g}'(u)du &= \int \cfrac{d}{du}F(g(u))du\\
  &=F(g(u))+C\\
  &=F(x)+C\\
  &=\int {F}'(x)dx\\
  &=\int f(x)dx\\
  \end{align} 

Substitution In Definite Integrals( version 2)
---------------------------------------------------------------------------
If :math:`{g}'` is continuous on the interval :math:`[a,b]` and :math:`f` is continuous on the range of 
:math:`g(u)=x`, then

.. math::
  \int_{g(a)}^{g(b)}f(x)dx=\int_{a}^{b}f(g(u))\cdot {g}'(u)du
  
Proof 

Let :math:`F` denote any antiderivative of :math:`f`. Then,  

.. math::
  \begin{aligned}
  \int_{a}^{b} f(g(x)) \cdot g^{\prime}(x) d x & =F(g(x))\Bigg ]_{x=a}^{x=b} \\
  & =F(g(b))-F(g(a)) \\
  & =F(u)\Bigg]_{u=g(a)}^{u=g(b)} \\
  & =\int_{g(a)}^{g(b)} f(u) d u .
  \end{aligned}  
  
Some Example  
--------------------
.. math::
  \int_{g(a)}^{g(b)}f(x)dx=\int_{a}^{b}f(g(u)){g}'(u) du

-
  
.. math::
  x=g(u)\\

-
  
.. math::
  \left\{\begin{matrix}
  u=a,& x=g(a)\\
  u=b,& x=g(b)\\
  \end{matrix}\right.
  
-
  
.. math::  
  \int_{g(a)}^{g(b)}f(x)dx=\int_{a}^{b}f(g(u))J(u) du\\ 
  
-
  
.. math::
  J(u)=\text{det}\begin{bmatrix}
  \cfrac{\partial x}{\partial u}
  \end{bmatrix}=\text{det}\begin{bmatrix}
  \cfrac{\partial g}{\partial u}
  \end{bmatrix}=\text{det}\begin{bmatrix}
  \cfrac{\text{d} g}{\text{d} u}
  \end{bmatrix}=\text{det}\begin{bmatrix}
   {g}'(u)
  \end{bmatrix}={g}'(u)\\
  
Example 1:

.. math::
  x=g(u,t)=(1-u)t+\frac{1}{2}ut^2+u \\
  
.. math::
  \left\{\begin{matrix}
  u=0,t=1,& x=g(u,t)=g(0,1)=1\\
  u=1,t=1,& x=g(u,t)=g(1,1)=\frac{1}{2}+1\\
  \end{matrix}\right.  
  
  
Let 

.. math::
  f(x,t)\equiv 1
  
then there is  

.. math::
  \int_{g(a,t)}^{g(b,t)}f(x,t)dx=\int_{1}^{1+\frac{1}{2} }1dx=\frac{1}{2}  
  
-
  
.. math:: 
  \cfrac{\partial g(u,t)}{\partial u}=\cfrac{\partial [(1-u)t+\frac{1}{2}ut^2+u]}{\partial u}=-t+\frac{1}{2}t^2+1  
  
-
  
.. math::   
  \cfrac{\partial g(u,t)}{\partial u}\Bigg|_{t=1}=\cfrac{\partial [(1-u)t+\frac{1}{2}ut^2+u]}{\partial u}\Bigg|_{t=1}=(-t+\frac{1}{2}t^2+1)\Bigg|_{t=1}=\frac{1}{2}  
  
-
  
.. math:: 
  \int_{g(a,t)}^{g(b,t)}f(x,t)dx=\int_{a(t)}^{b(t)}f(g(u,t),t)\cfrac{\partial g(u,t)}{\partial u} du\\
  
 
-
  
.. math:: 
  \int_{a(t)}^{b(t)}f(g(u,t),t)\cfrac{\partial g(u,t)}{\partial u} du=\int_{0}^{1}1\times\frac{1}{2} du=\frac{1}{2}  
  
Let  

.. math::
  f(x,t)\equiv x
  
then there is  

.. math::
  \int_{g(a,t)}^{g(b,t)}f(x,t)dx=\int_{1}^{1+\frac{1}{2} }xdx=\frac{1}{2}x^2\Bigg|_{1}^{1+\frac{1}{2}} 
  =\frac{1}{2}((\frac{3}{2})^2-1^{2})=\frac{5}{8}  

-
  
.. math::  
  \cfrac{\partial g(u,t)}{\partial u}\Bigg|_{t=1}=\cfrac{\partial [(1-u)t+\frac{1}{2}ut^2+u]}{\partial u}\Bigg|_{t=1}=(-t+\frac{1}{2}t^2+1)\Bigg|_{t=1}=\frac{1}{2}  
  
-
  
.. math::  
  \int_{g(a,t)}^{g(b,t)}f(x,t)dx=\int_{a(t)}^{b(t)}f(g(u,t),t)\cfrac{\partial g(u,t)}{\partial u} du  

-
  
.. math::  
  x=g(u,t)=(1-u)t+\frac{1}{2}ut^2+u \\

-
  
.. math::     
  \begin{align}
  \displaystyle \int_{a(t)}^{b(t)}f(g(u,t),t)\cfrac{\partial g(u,t)}{\partial u} du & = \int_{0}^{1}((1-u)1+\frac{1}{2}u1^2+u)\times\frac{1}{2} du\\
  \displaystyle & = \int_{0}^{1}(1+\frac{1}{2}u)\times\frac{1}{2} du\\
  \displaystyle& = \frac{1}{2}(u+\frac{1}{4}u^2)\Bigg|_{0}^{1}\\
  \displaystyle& = \frac{1}{2}(1+\frac{1}{4})\\
  & = \frac{5}{8}
  \end{align}
  
Example 2:

.. math:: 
  \begin{array}{l}
  x=g(u,t)\\
  f(x,t)=f(g(u,t),t)=\hat{f}(u,t)\\
  \cfrac{\text{d}\hat{f}(u,t)}{\text{d}t}\equiv \cfrac{\partial \hat{f}(u,t)}{\partial t}
  =\cfrac{\partial {f}(x,t)}{\partial t}+\cfrac{\partial {f}(x,t)}{\partial x}\cfrac{\partial {g}(u,t)}{\partial t}\\
  \cfrac{\text{d}\hat{f}}{\text{d}t}\equiv \cfrac{\partial \hat{f}}{\partial t}
  =\cfrac{\partial {f}}{\partial t}+\cfrac{\partial {f}}{\partial x}\cfrac{\partial {g}}{\partial t}
  \end{array}  
  
- 
  
.. math::
  x=g(u,t)=(1-u)t+\frac{1}{2}ut^2+u \\  
  
- 
  
.. math::
  \begin{array}{l}
  f(x,t)= xt \\  
  f(x,t)= xt=((1-u)t+\frac{1}{2}ut^2+u)t=\hat{f}(u,t) \\  
  \hat{f}(u,t)=((1-u)t+\frac{1}{2}ut^2+u)t=((1-u)t^2+\frac{1}{2}ut^3+ut)
  \end{array} 
  
- 
  
.. math::
  \cfrac{\text{d}\hat{f}(u,t)}{\partial t}\equiv\cfrac{\partial \hat{f}(u,t)}{\partial t}
  =\cfrac{\partial ((1-u)t^2+\cfrac{1}{2}ut^3+ut)}{\partial t}
  =(2(1-u)t+3\cfrac{1}{2}ut^2+u)
  
- 
  
.. math::  
  \cfrac{\partial f(x,t)}{\partial t}=x,\quad\cfrac{\partial f(x,t)}{\partial x}=t\\  
  
- 
  
.. math::
  \begin{align}
  \cfrac{\partial x}{\partial t} & = \cfrac{\partial g(u,t)}{\partial t}\\
   & = \cfrac{\partial ((1-u)t+\frac{1}{2}ut^2+u)}{\partial t}\\
   & = (1-u)+ut
  \end{align} 
  
- 
  
.. math:: 
  \begin{align}
  \cfrac{\text{d}f(x,t)}{\text{d} t} & = \cfrac{\partial {f}(x,t)}{\partial t}+\cfrac{\partial {f}(x,t)}{\partial x}\cfrac{\partial {g}(u,t)}{\partial t}\\
  & = x+t((1-u)+ut)\\
  & = x+((1-u)t+ut^2)\\
  & = ((1-u)t+\frac{1}{2}ut^2+u)+((1-u)t+ut^2)\\
  & = (2(1-u)t+\frac{3}{2}ut^2+u)\\
 \end{align}
 
A better understanding of this formula
------------------------------------------------ 

.. math:: 
  \begin{array}{c}
  x=g(u,t)=(1-u)t+\cfrac{1}{2}ut^2+u\\
  {g}'(u,t)=\cfrac{\partial g(u,t)}{\partial u}=-t+\cfrac{1}{2}t^2+1\\
  \displaystyle \int_{g(a)}^{g(b)}\cfrac{d  f(x,t)}{dt} dx=\int_{a}^{b}\cfrac{\partial \hat{f}(u,t)}{\partial t}{g}'(u,t) du
  \end{array}
  
Proof:

.. math::
  \begin{array}{c}
  \displaystyle h(x,t)\equiv \cfrac{d  f(x,t)}{dt}\\
  \displaystyle \hat{h}(u,t)\equiv \cfrac{\partial \hat{f}(u,t)}{\partial t}\\
  \end{array}
  
-
  
.. math::  
  h(x,t) = \hat{h}(u,t)

By the Chain Rule, :math:`F(g(u))` is an antiderivative of :math:`h(g(u))\cdot {g}'(u)` whenever :math:`F`
is an antiderivative of :math:`h`, beacuse

.. math::
  \begin{align}
  \cfrac{d}{du}F(g(u)) & = {F}'(g(u))\cdot {g}'(u)\\ & = h(g(u)){g}'(u)
  \end{align}
  
If we make the substitution :math:`x=g(u,t)`, then

.. math::
  \begin{align}
  \int h(g(u))\cdot {g}'(u)du &= \int \cfrac{d}{du}F(g(u))du\\
  &=F(g(u))+C\\
  &=F(x)+C\\
  &=\int {F}'(x)dx\\
  &=\int h(x)dx\\
  \end{align}  
  
-

.. math::
  \int h(x)dx = \int h(g(u))\cdot {g}'(u)du  
  
Continue to prove

If :math:`{g}'` is continuous on the interval :math:`[a,b]` and :math:`h` is continuous on the range of 
:math:`g(u,t)=x`, then

.. math::
  \int_{g(a)}^{g(b)}h(x)dx=\int_{a}^{b}h(g(u))\cdot {g}'(u)du
  
Proof 

Let :math:`F` denote any antiderivative of :math:`h`. Then,  

.. math::
  \begin{aligned}
  \int_{a}^{b} h(g(x)) \cdot g^{\prime}(x) d x & =F(g(x))\Bigg ]_{x=a}^{x=b} \\
  & =F(g(b))-F(g(a)) \\
  & =F(u)\Bigg]_{u=g(a)}^{u=g(b)} \\
  & =\int_{g(a)}^{g(b)} h(u) d u .
  \end{aligned} 


