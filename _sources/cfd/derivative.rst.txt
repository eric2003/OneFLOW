Fluid Derivative
==================================

The derivative (rate of change) of a field variable :math:`\phi(t,\mathbf{x}(t))`, which may be a scalar
or a vector quantity representing density, velocity, temperature, etc., with respect to
a fixed position in space is called the Eulerian derivative :math:`\cfrac{\partial \phi}{\partial t}` while the
derivative following a moving fluid parcel is called the Lagrangian, substantial, or
material derivative and is denoted by :math:`\cfrac{D \phi}{D t}`. The substantial derivative of variable :math:`\phi`, which can be derived through application of the chain rule to account for
changes induced by all independent variables along the path, is given by

.. math::
  \begin{align}
  \cfrac{D \phi}{D t} & = \cfrac{\partial \phi}{\partial t} \cfrac{d t}{d t}+\cfrac{\partial \phi}{\partial x} \underbrace{\cfrac{d x}{d t}}_{u}+\cfrac{\partial \phi}{\partial y} \underbrace{\cfrac{d y}{d t}}_{v}+\cfrac{\partial \phi}{\partial z} \underbrace{\cfrac{d z}{d t}}_{w} \\ & = \cfrac{\partial \phi}{\partial t}+u \cfrac{\partial \phi}{\partial x}+v \cfrac{\partial \phi}{\partial y}+w \cfrac{\partial \phi}{\partial z} \\ & = \cfrac{\partial \phi}{\underbrace{\partial t}_{\text {local rate }}}+\underbrace{\mathbf{v} \cdot \nabla \phi}_{\begin{array}{c}
  \text { convective rate } \\
  \text { of change }
  \end{array}} \\
  \end{align}
  
Eulerian Description
--------------------------
This method focuses on a fixed point in the flow field space, and it expresses some physical quantity of the fluid mass flowing through any fixed point in the space at each moment
as a function of the position :math:`\mathbf{r}` and time :math:`t` of that point.

.. math::
  \phi=\phi(\mathbf{r},t)
  
Lagrangian Description
--------------------------
This method focuses on a certain fluid particle, rather than a fixed point in space,
and it specifies the physical quantity of any particle in the process of motion as a function of
the vector variable marking the fluid particle :math:`\boldsymbol{\xi}` and time :math:`t`

.. math::
  \phi=\phi(\boldsymbol{\xi},t)
  
It is not difficult to understand that any definite fluid particle can be identified
by three corresponding numbers :math:`(\xi_{1},\xi_{2},\xi_{3})`. For example, these three numbers can be taken as the vector-radius components corresponding to
the position  :math:`P` of the particle at an initial moment, ie.

.. math::
  \boldsymbol{\xi}(P)=\mathbf{r}(P,t_{0})  
  
Interconversion of the two descriptions
---------------------------------------------------- 

.. math::
  \mathbf{r} = \mathbf{r}(\boldsymbol{\xi},t)

-

.. math::  
  \begin{align}
  x_{1}=\varphi_{1}(\xi_{1},\xi_{2},t)\\
  x_{2}=\varphi_{2}(\xi_{1},\xi_{2},t)\\
  \end{align}  
  
-

.. math::  
  \begin{align}
  \xi_{1}=\varphi^{-1}_{1}(x_{1},x_{2},t)\\
  \xi_{2}=\varphi^{-1}_{2}(x_{1},x_{2},t)\\
  \end{align}  

-

.. math::  
  \begin{array}{l}
  (\xi_{1},\xi_{2}) = const\\
  x_{1} = \varphi_{1}(\xi_{1},\xi_{2},t)=\hat{x}_{1}(t)\\
  x_{2} = \varphi_{2}(\xi_{1},\xi_{2},t)=\hat{x}_{2}(t)\\
  f(\xi_{1},\xi_{2},t) = f(\xi_{1}(x_{1},x_{2},t),\xi_{2}(x_{1},x_{2},t),t)=\hat{f}(x_{1},x_{2},t)={f}^{*}(\hat{x}_{1},\hat{x}_{2},t)\\
  \end{array}

-

.. math::
  \begin{align}
  \cfrac{\partial f(\xi_{1},\xi_{2},t)}{\partial t} & = \cfrac{\text{d} {f}^{*}(\hat{x}_{1}(t),\hat{x}_{2}(t),t)}{\text{d} t}\\
  & = \cfrac{\partial {f}^{*}(\hat{x}_{1}(t),\hat{x}_{2}(t),t)}{\partial t}\cfrac{\text{d} t}{\text{d} t}\\
  &+\cfrac{\partial {f}^{*}(\hat{x}_{1}(t),\hat{x}_{2}(t),t)}{\partial \hat{x}_{1}}\cfrac{\text{d} \hat{x}_{1}(t)}{\text{d} t}\\
  &+\cfrac{\partial {f}^{*}(\hat{x}_{1}(t),\hat{x}_{2}(t),t)}{\partial \hat{x}_{2}}\cfrac{\text{d} \hat{x}_{2}(t)}{\text{d} t}\\
  \end{align}
  
-

.. math::  
  \begin{align}
  \cfrac{\partial f(\xi_{1},\xi_{2},t)}{\partial t} & = \cfrac{\text{d} \hat{f}(x_{1},x_{2},t)}{\text{d} t}\\
  & = \cfrac{\partial \hat{f}(x_{1},x_{2},t)}{\partial t}\cfrac{\text{d} t}{\text{d} t}\\
  &+\cfrac{\partial \hat{f}(x_{1},x_{2},t)}{\partial x_{1}}\cfrac{\partial {x}_{1}(\xi_{1},\xi_{2},t)}{\partial t}\\
  &+\cfrac{\partial \hat{f}(x_{1},x_{2},t)}{\partial x_{2}}\cfrac{\partial {x}_{2}(\xi_{1},\xi_{2},t)}{\partial t}\\
  \end{align}  
  
ALE Description  
----------------------------------------------------   

.. math::  
  \begin{align}
  \chi_{1}=\psi^{-1}_{1}(\xi_{1},\xi_{2},t)\\
  \chi_{2}=\psi^{-1}_{2}(\xi_{1},\xi_{2},t)\\
  \end{align}  
  
-

.. math::  
  \begin{align}
  \xi_{1}=\psi_{1}(\chi_{1},\chi_{2},t)\\
  \xi_{2}=\psi_{2}(\chi_{1},\chi_{2},t)\\
  \end{align}  
  
-

.. math::  
  \begin{array}{l}
  (\xi_{1},\xi_{2})=const\\
  \chi_{1}=\psi^{-1}_{1}(\xi_{1},\xi_{2},t)=\hat{\psi}^{-1}_{1}(t)\\
  \chi_{2}=\psi^{-1}_{2}(\xi_{1},\xi_{2},t)=\hat{\psi}^{-1}_{2}(t)\\
  f(\xi_{1},\xi_{2},t) = f(\psi_{1}(\chi_{1},\chi_{2},t),\psi_{2}(\chi_{1},\chi_{2},t),t)=\overline{f}(\chi_{1},\chi_{2},t)=\overline{f}^{*}(\hat{\psi}^{-1}_{1}(t),\hat{\psi}^{-1}_{2}(t),t)\\
  \end{array}  
  
-

.. math::  
  \begin{align}
  \cfrac{\partial f(\xi_{1},\xi_{2},t)}{\partial t} & = \cfrac{\text{d} \overline{f}^{*}(\hat{\psi}^{-1}_{1}(t),\hat{\psi}^{-1}_{2}(t),t)}{\text{d} t}\\
  & = \cfrac{\partial \overline{f}^{*}(\hat{\psi}^{-1}_{1}(t),\hat{\psi}^{-1}_{2}(t),t)}{\partial t}\cfrac{\text{d} t}{\text{d} t}\\
  &+\cfrac{\partial \overline{f}^{*}(\hat{\psi}^{-1}_{1}(t),\hat{\psi}^{-1}_{2}(t),t)}{\partial \chi_{1}}\cfrac{\text{d} \hat{\psi}^{-1}_{1}(t)}{\text{d} t}\\
  &+\cfrac{\partial \overline{f}^{*}(\hat{\psi}^{-1}_{1}(t),\hat{\psi}^{-1}_{2}(t),t)}{\partial \chi_{2}}\cfrac{\text{d} \hat{\psi}^{-1}_{2}(t)}{\text{d} t}\\
  \end{align}  
  
-

.. math::  
  \begin{align}
  \cfrac{\partial f(\xi_{1},\xi_{2},t)}{\partial t} &= \cfrac{\text{d} \overline{f}(\chi_{1},\chi_{2},t)}{\text{d} t}\\
   &= \cfrac{\partial \overline{f}(\chi_{1},\chi_{2},t)}{\partial t}\cfrac{\text{d} t}{\text{d} t}\\
  &+\cfrac{\partial\overline{f}(\chi_{1},\chi_{2},t)}{\partial \chi_{1}}\cfrac{\partial \psi^{-1}_{1}(\xi_{1},\xi_{2},t)}{\partial t}\\
  &+\cfrac{\partial\overline{f}(\chi_{1},\chi_{2},t)}{\partial \chi_{2}}\cfrac{\partial \psi^{-1}_{2}(\xi_{1},\xi_{2},t)}{\partial t}\\
  \end{align}  

-

.. math::  
  \begin{align}
  x_{1} & = \varphi_{1}(\xi_{1},\xi_{2},t) = \phi_{1}(\chi_{1},\chi_{2},t)=\phi_{1}(\psi^{-1}_{1}(\xi_{1},\xi_{2},t),\psi^{-1}_{2}(\xi_{1},\xi_{2},t),t)\\
  x_{2} & = \varphi_{2}(\xi_{1},\xi_{2},t) = \phi_{2}(\chi_{1},\chi_{2},t)=\phi_{2}(\psi^{-1}_{1}(\xi_{1},\xi_{2},t),\psi^{-1}_{2}(\xi_{1},\xi_{2},t),t)\\
  \end{align}  
  
-

.. math:: 
  \begin{align}
  v_{1}&=v_{1}(\xi_{1},\xi_{2},t)=\cfrac{\partial \varphi_{1}(\xi_{1},\xi_{2},t)}{\partial t}\\
  &=\cfrac{\partial \phi_{1}(\chi_{1},\chi_{2},t)}{\partial t}
  +\cfrac{\partial \phi_{1}(\chi_{1},\chi_{2},t)}{\partial \chi_{1}}\cdot 
   \cfrac{\partial \psi^{-1}_{1}(\xi_{1},\xi_{2},t)}{\partial t}
  +\cfrac{\partial \phi_{1}(\chi_{1},\chi_{2},t)}{\partial \chi_{2}}\cdot 
   \cfrac{\partial \psi^{-1}_{2}(\xi_{1},\xi_{2},t)}{\partial t}\\
  v_{2}&=v_{2}(\xi_{1},\xi_{2},t)=\cfrac{\partial \varphi_{2}(\xi_{1},\xi_{2},t)}{\partial t}\\
  &=\cfrac{\partial \phi_{2}(\chi_{1},\chi_{2},t)}{\partial t}
  +\cfrac{\partial \phi_{2}(\chi_{1},\chi_{2},t)}{\partial \chi_{1}}\cdot 
   \cfrac{\partial \psi^{-1}_{1}(\xi_{1},\xi_{2},t)}{\partial t}
  +\cfrac{\partial \phi_{2}(\chi_{1},\chi_{2},t)}{\partial \chi_{2}}\cdot 
   \cfrac{\partial \psi^{-1}_{2}(\xi_{1},\xi_{2},t)}{\partial t}\\
  \end{align}  
  
-

.. math::   
  \begin{align}
  v_{j} & = v_{j}(\xi_{1},\xi_{2},t)  = \cfrac{\partial \varphi_{j}(\xi_{1},\xi_{2},t)}{\partial t}\\
   & = \cfrac{\partial \phi_{j}(\chi_{1},\chi_{2},t)}{\partial t}
  +\cfrac{\partial \phi_{j}(\chi_{1},\chi_{2},t)}{\partial \chi_{1}}\cdot 
   \cfrac{\partial \psi^{-1}_{1}(\xi_{1},\xi_{2},t)}{\partial t}
  +\cfrac{\partial \phi_{j}(\chi_{1},\chi_{2},t)}{\partial \chi_{2}}\cdot 
   \cfrac{\partial \psi^{-1}_{2}(\xi_{1},\xi_{2},t)}{\partial t}\\
  \end{align} 
  
-

.. math:: 
  \begin{align}
  v_{j} & = v_{j}(\xi_{1},\xi_{2},t)  = \cfrac{\partial \varphi_{j}(\xi_{1},\xi_{2},t)}{\partial t}\\
   & = \cfrac{\partial \phi_{j}(\chi_{1},\chi_{2},t)}{\partial t}
  +\cfrac{\partial \phi_{j}(\chi_{1},\chi_{2},t)}{\partial \chi_{k}}\cdot 
   \cfrac{\partial \psi^{-1}_{k}(\xi_{1},\xi_{2},t)}{\partial t}\\
  \end{align} 
  
-

.. math::
  \begin{align}
  v_{j} & = v_{j}(\boldsymbol\xi,t)  =
  \hat v_{j}(\boldsymbol\chi,t)
  +\cfrac{\partial x_{j}(\boldsymbol\chi,t)}{\partial \chi_{k}(\boldsymbol\xi,t)}\cdot 
  \cfrac{\partial \chi_{k}(\boldsymbol\xi,t)}{\partial t}\\
  \end{align}  

-

.. math::
  w_{k}=w_{k}(\xi_{1},\xi_{2},t)=w_{k}(\boldsymbol\xi,t)=\cfrac{\partial \chi_{k}(\xi_{1},\xi_{2},t)}{\partial t}=\cfrac{\partial \chi_{k}(\boldsymbol\xi,t)}{\partial t}  

-

.. math::
  \begin{align}
  c_{j}=v_{j}(\boldsymbol\xi,t)-\hat v_{j}(\boldsymbol\chi,t)&=\cfrac{\partial x_{j}(\boldsymbol\chi,t)}{\partial \chi_{k}(\boldsymbol\xi,t)}\cdot 
   \cfrac{\partial \chi_{k}(\boldsymbol\xi,t)}{\partial t}\\
  &=\cfrac{\partial x_{j}(\boldsymbol\chi,t)}{\partial \chi_{k}(\boldsymbol\xi,t)}\cdot w_{k}(\boldsymbol\xi,t)\\
  &=\cfrac{\partial x_{j}}{\partial \chi_{k}}\cdot w_{k}\\
  \end{align}  
  
-

.. math::
  \underbrace{\mathbf{c}}_{\text{convective velocity}}=\underbrace{\mathbf{v}}_{\text{fluid velocity}}-\underbrace{\hat{\mathbf{v}}}_{\text{mesh velocity}}=\cfrac{\partial x_{j}}{\partial \chi_{k}}\cdot \underbrace{w_{k}}_{\text{referential velocity}}  
  
Material, spatial, and referential time derivatives 
------------------------------------------------------

In order to relate the time derivative in the material, spatial,
and referential domains, let a scalar physical quantity be
described by :math:`f(\mathbf{x},t)`, :math:`f^{*}(\boldsymbol{\chi},t)`, and :math:`f^{**}(\boldsymbol{\xi},t)` in the spatial,
referential, and material domains respectively. Stars are
employed to emphasize that the functional forms are, in
general, different.

.. math::
  f^{**}(\boldsymbol{\xi},t)=f(\mathbf{x},t)=f(\boldsymbol{\varphi}(\boldsymbol{\xi},t),t)
  
-
  
.. math::
  \cfrac{\text{d}F(t)}{\text{d}t}=\cfrac{\partial f^{**}(\boldsymbol{\xi},t)}{\partial t}
  =\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\cfrac{\partial f(\mathbf{x},t)}{\partial x_{j}(\boldsymbol{\xi},t)}\cdot\cfrac{\partial x_{j}(\boldsymbol{\xi},t)}{\partial t} \\
  
-
  
.. math::
  \text{grad }f=\nabla f =\begin{bmatrix}
  \cfrac{\partial f(\mathbf{x},t)}{\partial \mathbf{x}}
  \end{bmatrix}^{\text{T}}
  =\begin{bmatrix}
   \cfrac{\partial f(\mathbf{x},t)}{\partial {x}_{1}}\\
   \cfrac{\partial f(\mathbf{x},t)}{\partial {x}_{2}}
  \end{bmatrix} 
  
-
  
.. math::
  \mathbf{v}=\begin{bmatrix}
  \cfrac{\partial x_{1}(\boldsymbol{\xi},t)}{\partial t}\\
  \cfrac{\partial x_{2}(\boldsymbol{\xi},t)}{\partial t}
  \end{bmatrix}
  =\begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial t}\Bigg|_{\boldsymbol{\xi}}\\
  \cfrac{\partial x_{2}}{\partial t}\Bigg|_{\boldsymbol{\xi}}
  \end{bmatrix}
  
-
  
.. math::
  \cfrac{\partial f^{**}}{\partial t}
  =\cfrac{\partial f}{\partial t}+\cfrac{\partial f}{\partial\mathbf{x}}\cdot\mathbf{v}  
  
Note that this is the well-known equation that relates the
material and the spatial time derivatives. Dropping the stars
to ease the notation, this relation is finally cast as

.. math::
  \begin{align}
  \cfrac{\text{d} f}{\text{d} t}=\cfrac{\partial f}{\partial t}\Bigg|_{\boldsymbol{\xi}}
  =\cfrac{\partial f}{\partial t}\Bigg|_{\mathbf{x}}+\mathbf{v}\cdot\nabla f\\
  \end{align}
  
-
  
.. math:: 
  \begin{align}
  \cfrac{\text{d} f}{\text{d} t}=\cfrac{\partial f}{\partial t}\Bigg|_{\boldsymbol{\xi}}
  =\cfrac{\partial f}{\partial t}\Bigg|_{\mathbf{x}}+\mathbf{v}\cdot\nabla f\\
  \end{align}
  
-
  
.. math:: 
  \cfrac{\text{d} f}{\text{d} t}=\cfrac{\partial f}{\partial t}+\mathbf{v}\cdot\nabla f  
  
With the help of mapping :math:`\boldsymbol{\psi}`, the transformation from
the referential description :math:`f^{*}(\boldsymbol{\chi},t)` of the scalar physical
quantity to the material description :math:`f^{**}(\boldsymbol{\xi},t)` can be written
as
  
.. math:: 
  f^{**}(\boldsymbol{\xi},t)=f^{*}(\boldsymbol{\chi},t)=f^{*}(\boldsymbol{\psi}^{-1}(\boldsymbol{\xi},t),t)  

-
  
.. math:: 
  \cfrac{\text{d}F(t)}{\text{d}t}=\cfrac{\partial f^{**}(\boldsymbol{\xi},t)}{\partial t}
  =\cfrac{\partial f^{*}(\boldsymbol{\chi},t)}{\partial t}+\cfrac{\partial f^{*}(\boldsymbol{\chi},t)}{\partial \chi_{j}(\boldsymbol{\xi},t)}\cdot\cfrac{\partial \chi_{j}(\boldsymbol{\xi},t)}{\partial t} \\

-
  
.. math:: 
  \cfrac{\partial f^{**}(\boldsymbol{\xi},t)}{\partial t}
  =\cfrac{\partial f^{*}(\boldsymbol{\chi},t)}{\partial t}+\cfrac{\partial f^{*}(\boldsymbol{\chi},t)}{\partial\boldsymbol{\chi}(\boldsymbol{\xi},t)}\cdot\mathbf{w}  
  
-
 
.. math:: 
  \mathbf{w}=\begin{bmatrix}
  \cfrac{\partial \chi_{1}(\boldsymbol{\xi},t)}{\partial t}\\
  \cfrac{\partial \chi_{2}(\boldsymbol{\xi},t)}{\partial t}\\
  \end{bmatrix}=
  \begin{bmatrix}
  \cfrac{\partial \chi_{1}}{\partial t}\Bigg|_{\boldsymbol{\xi}}\\
  \cfrac{\partial \chi_{2}}{\partial t}\Bigg|_{\boldsymbol{\xi}}\\
  \end{bmatrix}
  
-
 
.. math::
  \cfrac{\partial f^{**}}{\partial t}
  =\cfrac{\partial f^{*}}{\partial t}+\cfrac{\partial f^{*}}{\partial\boldsymbol{\chi}}\cdot\mathbf{w}  
  
Note that this equation relates the material and the referential time derivatives. However, it also requires the
evaluation of the gradient of the considered quantity in
the referential domain. This can be done, but in computational mechanics it is usually easier to work in the spatial
(or material) domain.
  
.. math::
  f^{**}(\boldsymbol{\xi},t)=f^{*}(\boldsymbol{\chi},t)=f(\mathbf{x},t) =f(\boldsymbol{\phi}(\boldsymbol{\chi},t),t)  
  
-
 
.. math::
  f^{*}(\boldsymbol{\chi},t)=f(\boldsymbol{\phi}(\boldsymbol{\chi},t),t) 

-
 
.. math::
  \cfrac{\partial f^{*}(\boldsymbol{\chi},t)}{\partial {\chi}_{j}}=
  \cfrac{\partial f(\boldsymbol{\phi}(\boldsymbol{\chi},t),t)}{\partial {x}_{k}}\cdot 
  \cfrac{\partial {\phi}_{k}(\boldsymbol{\chi},t)}{\partial {\chi}_{j}}
  
-
 
.. math::
  \cfrac{\partial f^{*}(\boldsymbol{\chi},t)}{\partial {\chi}_{j}}=
  \cfrac{\partial f(\mathbf{x},t)}{\partial {x}_{k}}\cdot 
  \cfrac{\partial {x}_{k}(\boldsymbol{\chi},t)}{\partial {\chi}_{j}}\\ 
  
-
 
.. math::
  \cfrac{\text{d}F(t)}{\text{d}t}=\cfrac{\partial f^{**}(\boldsymbol{\xi},t)}{\partial t}
  =\cfrac{\partial f^{*}(\boldsymbol{\chi},t)}{\partial t}+\cfrac{\partial f^{*}(\boldsymbol{\chi},t)}{\partial \chi_{j}(\boldsymbol{\xi},t)}\cdot\cfrac{\partial \chi_{j}(\boldsymbol{\xi},t)}{\partial t} \\
  =\cfrac{\partial f^{*}(\boldsymbol{\chi},t)}{\partial t}+\cfrac{\partial f(\mathbf{x},t)}{\partial {x}_{k}}\cdot 
  \cfrac{\partial {x}_{k}(\boldsymbol{\chi},t)}{\partial {\chi}_{j}}\cdot\cfrac{\partial \chi_{j}(\boldsymbol{\xi},t)}{\partial t} \\  

Thus, using the definition of :math:`w`, the previous equation may be rearranged into

.. math::
  \begin{align}
  c_{j}=v_{j}(\boldsymbol\xi,t)-\hat v_{j}(\boldsymbol\chi,t)&=\cfrac{\partial x_{j}(\boldsymbol\chi,t)}{\partial \chi_{k}(\boldsymbol\xi,t)}\cdot 
   \cfrac{\partial \chi_{k}(\boldsymbol\xi,t)}{\partial t}\\
  &=\cfrac{\partial x_{j}(\boldsymbol\chi,t)}{\partial \chi_{k}(\boldsymbol\xi,t)}\cdot w_{k}(\boldsymbol\xi,t)\\
  &=\cfrac{\partial x_{j}}{\partial \chi_{k}}\cdot w_{k}\\
  \end{align} 
  
-
  
.. math::
  \begin{align}
  c_{k}=v_{k}(\boldsymbol\xi,t)-\hat v_{k}(\boldsymbol\chi,t)
  &=\cfrac{\partial x_{k}(\boldsymbol\chi,t)}{\partial \chi_{j}(\boldsymbol\xi,t)}\cdot 
   \cfrac{\partial \chi_{j}(\boldsymbol\xi,t)}{\partial t}\\
  &=\cfrac{\partial x_{k}(\boldsymbol\chi,t)}{\partial \chi_{j}(\boldsymbol\xi,t)}\cdot w_{j}(\boldsymbol\xi,t)\\
  &=\cfrac{\partial x_{k}}{\partial \chi_{j}}\cdot w_{j}\\
  \end{align}
  
-
  
.. math::
  \cfrac{\text{d}F(t)}{\text{d}t}=\cfrac{\partial f^{**}}{\partial t}
  =\cfrac{\partial f^{*}}{\partial t}+\cfrac{\partial f}{\partial \mathbf{x}}\cdot \mathbf{c}
  
The fundamental ALE relation between material time
derivatives, referential time derivatives, and spatial gradient
is finally cast as (stars dropped) 

.. math::
  \cfrac{\text{d} f}{\text{d} t}=\cfrac{\partial f}{\partial t}\Bigg|_{\boldsymbol{\xi}}
  =\cfrac{\partial f}{\partial t}\Bigg|_{\boldsymbol{\chi}}+\cfrac{\partial f}{\partial \mathbf{x}}\cdot \mathbf{c}
  =\cfrac{\partial f}{\partial t}\Bigg|_{\boldsymbol{\chi}}+\mathbf{c} \cdot \nabla {f}\\ 
  
  
and shows that the time derivative of the physical quantity
f for a given particle :math:`\boldsymbol{\xi}`, that is, its material derivative,
is its local derivative (with the reference coordinate :math:`\boldsymbol{\chi}` held
fixed) plus a convective term taking into account the relative
velocity :math:`\mathbf{c}` between the material and the reference system.  
  