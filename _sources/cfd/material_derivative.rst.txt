Material Time Derivative
==================================

The material time derivative of an integral is the rate of change of an integral on a material
domain. A material domain moves with the material, so that the material points on the
boundary remain on the boundary and no mass flux occurs across the boundaries. A material
domain is analogous to a Lagrangian mesh; a Lagrangian element or group of Lagrangian elements is a nice example of a material domain. The various forms for material time derivatives
of integrals are called Reynolds’ transport theorem.

The material time derivative of an integral is defined by

.. math::
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}f\text{d}\Omega=\lim_{\Delta t \to 0} \cfrac{1}{\Delta t}
  \left(\int\limits_{\Omega_{\tau+\Delta t}}f(\mathbf{x},{\tau+\Delta t})\text{d}\Omega-\int\limits_{\Omega_{\tau}}f(\mathbf{x},{\tau})\text{d}\Omega\right)
  
where :math:`\Omega_{\tau}` is the spatial domain at time :math:`t` and :math:`\Omega_{\tau+\Delta t}` is the spatial domain occupied by the same
material points at time :math:`\Omega_{\tau+\Delta t}`. The notation on the left-hand side is a little confusing because it
appears to refer to a single spatial domain. However, in this notation, which is standard, the material derivative on the integral implies that the domain Ω is a material domain. We now
transform both integrals on the RHS to the reference domain using

.. math::
  \int\limits_{\Omega}f(\mathbf{x}, t)\text{d}\Omega=\int\limits_{\Omega_{0}}f(\boldsymbol\phi(\mathbf{X},t), t)J\text{d}\Omega_{0}
  
-
  
.. math::
  \mathbf{F}=\cfrac{\partial (x_{1},x_{2},x_{3})}{\partial (X_{1},X_{2},X_{3})}=
  \begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial X_{1}}& \cfrac{\partial x_{1}}{\partial X_{2}} & \cfrac{\partial x_{1}}{\partial X_{3}}\\
  \cfrac{\partial x_{2}}{\partial X_{1}}& \cfrac{\partial x_{2}}{\partial X_{2}} & \cfrac{\partial x_{2}}{\partial X_{3}}\\
  \cfrac{\partial x_{3}}{\partial X_{1}}& \cfrac{\partial x_{3}}{\partial X_{2}} & \cfrac{\partial x_{3}}{\partial X_{3}}\\
  \end{bmatrix}
  =\begin{bmatrix}
  \cfrac{\partial x}{\partial X}& \cfrac{\partial x}{\partial Y} & \cfrac{\partial x}{\partial Z}\\
  \cfrac{\partial y}{\partial X}& \cfrac{\partial y}{\partial Y} & \cfrac{\partial y}{\partial Z}\\
  \cfrac{\partial z}{\partial X}& \cfrac{\partial z}{\partial Y} & \cfrac{\partial z}{\partial Z}\\
  \end{bmatrix}  

-
  
.. math::
  J=J(X_{1},X_{2},X_{3})=\text{det}(\mathbf{F})=\left | \cfrac{ \partial(x_{1},x_{2},x_{3})}{\partial (X_{1},X_{2},X_{3})} \right |=
  \begin{vmatrix}
  \cfrac{\partial x_{1}}{\partial X_{1}}& \cfrac{\partial x_{1}}{\partial X_{2}} & \cfrac{\partial x_{1}}{\partial X_{3}}\\
  \cfrac{\partial x_{2}}{\partial X_{1}}& \cfrac{\partial x_{2}}{\partial X_{2}} & \cfrac{\partial x_{2}}{\partial X_{3}}\\
  \cfrac{\partial x_{3}}{\partial X_{1}}& \cfrac{\partial x_{3}}{\partial X_{2}} & \cfrac{\partial x_{3}}{\partial X_{3}}\\
  \end{vmatrix}
  =\begin{vmatrix}
  \cfrac{\partial x}{\partial X}& \cfrac{\partial x}{\partial Y} & \cfrac{\partial x}{\partial Z}\\
  \cfrac{\partial y}{\partial X}& \cfrac{\partial y}{\partial Y} & \cfrac{\partial y}{\partial Z}\\
  \cfrac{\partial z}{\partial X}& \cfrac{\partial z}{\partial Y} & \cfrac{\partial z}{\partial Z}\\
  \end{vmatrix}
  
-
  
.. math::
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}f\text{d}\Omega=\lim_{\Delta t \to 0} \cfrac{1}{\Delta t}
  \left(\int\limits_{\Omega_{0}}f(\mathbf{X},{\tau+\Delta t})J(\mathbf{X},{\tau+\Delta t})\text{d}\Omega_{0}-\int\limits_{\Omega_{0}}f(\mathbf{X},{\tau})J(\mathbf{X},{\tau})\text{d}\Omega_{0}\right)  
  
With this change in the domain of integration, :math:`f` becomes a function of the material coordinates, that is,  
:math:`f(\boldsymbol\phi(\mathbf{X},t), t)\equiv f\circ \boldsymbol\phi`.

Since the domain of integration is now independent of time, we can pull the limit operation
inside the integral and take the limit, which yields

.. math::
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}f\text{d}\Omega=\int\limits_{\Omega_{0}}\cfrac{\partial}{\partial t}( f(\mathbf{X},{t})J(\mathbf{X},{t}))\text{d}\Omega_{0}
  
The partial derivative with respect to time in the integrand is a material time derivative since
the independent space variables are the material coordinates. We next use the product rule for
derivatives on the previous:

.. math::
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}f\text{d}\Omega=\int\limits_{\Omega_{0}}\cfrac{\partial}{\partial t}( f(\mathbf{X},{t})J(\mathbf{X},{t}))\text{d}\Omega_{0}
  =\int\limits_{\Omega_{0}}\left(\cfrac{\partial f(\mathbf{X},{t})}{\partial t}J(\mathbf{X},{t})+ f(\mathbf{X},{t})\cfrac{\partial J(\mathbf{X},{t})}{\partial t}\right)\text{d}\Omega_{0}
  
Bearing in mind that the partial time derivatives are material time derivatives, we can obtain  

.. math:: 
  \cfrac{\text{d}J}{\text{d}t}=J ( \text{div } \mathbf{v})=J ( \nabla \cdot \mathbf{v})=J(\cfrac{\partial v_{k} }{\partial x_{k}})
  
We can now transform the RHS integral to the current domain and change the
independent variables to an Eulerian description, which gives
  
.. math::
  \begin{align}
  &\int\limits_{\Omega_{0}}\left(\cfrac{\partial f(\mathbf{X},{t})}{\partial t}J(\mathbf{X},{t})+ f(\mathbf{X},{t})\cfrac{\partial J(\mathbf{X},{t})}{\partial t}\right)\text{d}\Omega_{0}\\
  = &\int\limits_{\Omega_{0}}\left(\cfrac{\partial f(\mathbf{X},{t})}{\partial t}J(\mathbf{X},{t})+ f(\mathbf{X},{t})J(\mathbf{X},{t}) ( \cfrac{\partial v_{k} }{\partial x_{k}})\right)\text{d}\Omega_{0}\\
  = &\int\limits_{\Omega_{0}}\left(\cfrac{\partial f(\mathbf{X},{t})}{\partial t}+ f(\mathbf{X},{t}) ( \cfrac{\partial v_{k} }{\partial x_{k}})\right)J(\mathbf{X},{t})\text{d}\Omega_{0}\\
  = &\int\limits_{\Omega}\left(\cfrac{\text{d}f(\mathbf{x},{t})}{\text{d} t}+ f(\mathbf{x},{t}) ( \cfrac{\partial v_{k} }{\partial x_{k}})\right)\text{d}\Omega\\
  \end{align}

-
  
.. math::  
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}f\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\text{d}f(\mathbf{x},{t})}{\text{d} t}+ f(\mathbf{x},{t}) ( \cfrac{\partial v_{k} }{\partial x_{k}})\right)\text{d}\Omega\\  
  
This is one form of Reynolds’ transport theorem.  

An alternative form of Reynolds’ transport theorem can be obtained by the definition of the
material time derivative. This gives

.. math::  
  \begin{align}
  \displaystyle\cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}f\text{d}\Omega & = \int\limits_{\Omega}\left(\cfrac{\partial f(\mathbf{x},{t})}{\partial t}+v_{k}\cfrac{\partial f(\mathbf{x},{t})}{\partial x_{k}}+ f(\mathbf{x},{t}) ( \cfrac{\partial v_{k} }{\partial x_{k}})\right)\text{d}\Omega\\
  &  = \int\limits_{\Omega}\left(\cfrac{\partial f(\mathbf{x},{t})}{\partial t}+\cfrac{\partial (v_{k}f(\mathbf{x},{t}))}{\partial x_{k}}\right)\text{d}\Omega
  \end{align}  
  
which can be written in tensor form as  

.. math:: 
  \displaystyle\cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}f\text{d}\Omega 
  = \int\limits_{\Omega}\left(\cfrac{\partial f(\mathbf{x},{t})}{\partial t}+\text{div } \left(\mathbf{v}f(\mathbf{x},{t})\right)\right)\text{d}\Omega
  
Using Gauss’s theorem on the second term of the RHS, which gives

.. math:: 
  \cfrac{\text{d}}{\text{d} t} \int\limits_{\Omega} f \text{d} \Omega=\int\limits_{\Omega} \cfrac{\partial f}{\partial t} \text{d} \Omega+\int\limits_{\Gamma} f v_{i} n_{i} \text{d} \Gamma \quad \text { or } \quad \frac{\text{d}}{\text{d} t} \int\limits_{\Omega} f \text{d} \Omega=\int\limits_{\Omega} \frac{\partial f}{\partial t} \text{d} \Omega+\int\limits_{\Gamma} f \mathbf{v} \cdot \mathbf{n} \text{d} \Gamma
  
where the product :math:`f\mathbf{v}` is assumed to be :math:`C_{0}` in :math:`\Omega`. Reynolds’ transport theorem, which in the
above has been given for a scalar, applies to a tensor of any order. Thus to apply it to a firstorder tensor (vector) :math:`g_{k}`, replace :math:`f` by :math:`g_{k}`, which gives  

.. math:: 
  \cfrac{\text{d}}{\text{d} t} \int\limits_{\Omega} g_{k} \text{d} \Omega=\int\limits_{\Omega}\left(\frac{\partial g_{k}}{\partial t}+\frac{\partial\left(v_{i} g_{k}\right)}{\partial x_{i}}\right) \text{d} \Omega
  
Mass Conservation
-----------------------
The mass :math:`m(\Omega)` of a material domain :math:`\Omega` is given by

.. math::
  m(\Omega)=\int\limits_{\Omega}\rho(\mathbf{x},t)\text{d}\Omega
  
where :math:`\rho(\mathbf{x},t)` is the density. Mass conservation requires that the mass of any material domain
be constant, since no material flows through the boundaries of a material domain and we are
not considering mass to energy conversion. Therefore, according to the principle of mass
conservation, the material time derivative of :math:`m(\Omega)` vanishes, that is,  

.. math::
  \cfrac{\text{d}m}{\text{d} t}=\cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}\rho\text{d}{\Omega}=0
  
Applying Reynolds’ theorem, yields
 
.. math::
  \int\limits_{\Omega}\left(\cfrac{\text{d}\rho}{\text{d} t} +\rho\text{ div }\mathbf{v}\right) \text{d} \Omega=0
  
Since this holds for any subdomain :math:`\Omega`, so:  

.. math::
  \cfrac{\text{d}\rho}{\text{d} t} +\rho\text{ div }\mathbf{v}=0
  
- 
 
.. math::
  \begin{align}
  \displaystyle \cfrac{\text{d}\rho}{\text{d} t} +\rho\text{ div }\mathbf{v} & = 0\\
  \displaystyle \cfrac{\text{d}\rho}{\text{d} t} +\rho\cfrac{\partial v_{i}}{\partial x_{i}} & = 0\\
  \displaystyle \cfrac{\text{d}\rho}{\text{d} t} +\rho v_{i,i} & = 0\\
  \displaystyle \dot{\rho} +\rho v_{i,i} & = 0\\
  \end{align}  
  
The above is the equation of mass conservation, often called the continuity equation. It is a
first-order partial differential equation.

The continuity equation can be written in the form

.. math::
  \cfrac{\partial\rho}{\partial t}+v_{i} \cfrac{\partial\rho}{\partial x_{i}}+\rho\cfrac{\partial v_{i}}{\partial x_{i}}
  = \cfrac{\partial\rho}{\partial t}+ \cfrac{\partial(\rho v_{i})}{\partial x_{i}}
  = \cfrac{\partial\rho}{\partial t}+ (\rho v_{i})_{,i}=0
  
This is called the conservative form of the mass conservation equation. It is often preferred in
computational fluid dynamics because discretizations are thought to enforce mass
conservation more accurately.  

Conservation of Linear Momentum
-----------------------------------------
The equation emanating from the principle of linear momentum conservation is a key equation
in nonlinear finite element procedures. Linear momentum conservation is equivalent to
Newton’s second law of motion, which relates the forces acting on a body to its acceleration.
The principle is often called the momentum conservation principle, or the balance of
momentum principle.

We will here state the principle in integral form and then derive an equivalent partial
differential equation. We consider an arbitrary domain :math:`\Omega` with boundary :math:`\Gamma` subjected to body
forces :math:`\rho \mathbf{b}` and to surface tractions :math:`\mathbf{t}`, where :math:`\mathbf{b}` is a force per unit mass and :math:`\mathbf{t}` is a force per unit
area. The total force is given by

.. math::
  \mathbf{f}(t)=\int\limits_{\Omega}\rho\mathbf{b}(\mathbf{x},t)\text{d}{\Omega}+\int\limits_{\Omega}\mathbf{t}(\mathbf{x},t)\text{d}{\Gamma}
  
The linear momentum is given by

.. math::
  \mathbf{p}(t)=\int\limits_{\Omega}\rho\mathbf{v}(\mathbf{x},t)\text{d}{\Omega}
  
Newton’s second law of motion for a continuum, the momentum conservation principle,
states that the material time derivative of the linear momentum equals the net force.

.. math::
  \cfrac{\text{d}\mathbf{p}}{\text{d} t}=\mathbf{f}
  \Rightarrow \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}\rho\mathbf{v}(\mathbf{x},t)\text{d}{\Omega}
  =\int\limits_{\Omega}\rho\mathbf{b}(\mathbf{x},t)\text{d}{\Omega}
  +\int\limits_{\Omega}\mathbf{t}(\mathbf{x},t)\text{d}{\Gamma}
  
Reynolds’ transport theorem applied to the LHS integral gives

.. math::
  \begin{align}
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}\rho\mathbf{v}\text{d}{\Omega} & = \int\limits_{\Omega}\left(\cfrac{\text{d}}{\text{d} t}(\rho\mathbf{v})+(\text{div }\mathbf{v})(\rho\mathbf{v})\right)\text{d}{\Omega}\\
  &=\int\limits_{\Omega}\left(\rho\cfrac{\text{d}\mathbf{v}}{\text{d} t}+\mathbf{v}\cfrac{\text{d}\rho}{\text{d} t}+(\text{div }\mathbf{v})(\rho\mathbf{v})\right)\text{d}{\Omega}\\
  &=\int\limits_{\Omega}\left(\rho\cfrac{\text{d}\mathbf{v}}{\text{d} t}+\mathbf{v}(\cfrac{\text{d}\rho}{\text{d} t}+\rho\text{ div }\mathbf{v})\right)\text{d}{\Omega}\\
  \end{align}  
  
The term multiplying the velocity in the RHS of the above can be recognized as the continuity equation, which vanishes, giving  

.. math::
  \cfrac{\text{d}\rho}{\text{d} t} +\rho\text{ div }\mathbf{v}=0
  
-  
  
.. math::  
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}\rho\mathbf{v}\text{d}{\Omega} 
  &=\int\limits_{\Omega}\rho\cfrac{\text{d}\mathbf{v}}{\text{d} t}\text{d}{\Omega}\\ 

To convert the second term on the RHS to a domain integral, we invoke Cauchy’s
relation and Gauss’s theorem in sequence, giving  

.. math::  
  \int\limits_{\Gamma}\mathbf{t}\text{d}{\Gamma}=\int\limits_{\Gamma}\mathbf{n}\cdot\boldsymbol{\sigma}\text{d}{\Gamma}
  =\int\limits_{\Omega}\nabla\cdot\boldsymbol{\sigma}\text{d}{\Omega} \quad\times
  
- 

.. math::
  \int\limits_{\Gamma}\mathbf{t}\text{d}{\Gamma}=\int\limits_{\Gamma}\boldsymbol{\sigma}\cdot\mathbf{n}\text{d}{\Gamma}
  =\int\limits_{\Omega}\text{div }\boldsymbol{\sigma}\text{d}{\Omega}
  =\int\limits_{\Omega}\nabla\cdot ([\boldsymbol{\sigma}]^{\text{T}})\text{d}{\Omega}
  \quad\checkmark
  
-  
  
.. math::  
  \boldsymbol\sigma=
  \begin{bmatrix}
  \sigma_{11}& \sigma_{12} & \sigma_{13}\\
  \sigma_{21}& \sigma_{22} & \sigma_{23}\\
  \sigma_{31}& \sigma_{32} & \sigma_{33}\\
  \end{bmatrix}   
  
- 

.. math::
  \text{div}\mathbf{T} = \nabla \cdot (\mathbf{T}^{\text T})
  =\cfrac{\partial T_{ij}}{\partial x_{j}}=
  \begin{bmatrix}
   \cfrac{\partial {T_{11}}}{\partial x_{1}}
  +\cfrac{\partial {T_{12}}}{\partial x_{2}}
  +\cfrac{\partial {T_{13}}}{\partial x_{3}}\\
   \cfrac{\partial {T_{21}}}{\partial x_{1}}
  +\cfrac{\partial {T_{22}}}{\partial x_{2}}
  +\cfrac{\partial {T_{23}}}{\partial x_{3}}\\
   \cfrac{\partial {T_{31}}}{\partial x_{1}}
  +\cfrac{\partial {T_{32}}}{\partial x_{2}}
  +\cfrac{\partial {T_{33}}}{\partial x_{3}}\\
  \end{bmatrix}
  
- 

.. math::
  \text{div }\boldsymbol{\sigma} = \nabla \cdot (\boldsymbol{\sigma}^{\text T})
  =\cfrac{\partial {\sigma}_{ij}}{\partial x_{j}}=
  \begin{bmatrix}
   \cfrac{\partial {{\sigma}_{11}}}{\partial x_{1}}
  +\cfrac{\partial {{\sigma}_{12}}}{\partial x_{2}}
  +\cfrac{\partial {{\sigma}_{13}}}{\partial x_{3}}\\
   \cfrac{\partial {{\sigma}_{21}}}{\partial x_{1}}
  +\cfrac{\partial {{\sigma}_{22}}}{\partial x_{2}}
  +\cfrac{\partial {{\sigma}_{23}}}{\partial x_{3}}\\
   \cfrac{\partial {{\sigma}_{31}}}{\partial x_{1}}
  +\cfrac{\partial {{\sigma}_{32}}}{\partial x_{2}}
  +\cfrac{\partial {{\sigma}_{33}}}{\partial x_{3}}\\
  \end{bmatrix}  

-  
  
.. math::  
  \int\limits_{\Gamma}{t}_{j}\text{d}{\Gamma}=\int\limits_{\Gamma}{n}_{i}\cdot{\sigma}_{ij}\text{d}{\Gamma}
  =\int\limits_{\Omega}\cfrac{\partial {\sigma}_{ij}}{\partial x_{i}}\text{d}{\Omega} \quad\times\\
  
-  
  
.. math::  
  t_{j}=n_{1}\sigma_{1j}+n_{2}\sigma_{2j}+n_{3}\sigma_{3j} \quad\times\\  
  
-  
  
.. math::  
  \begin{array}{c}
  t_{1}=n_{1}\sigma_{11}+n_{2}\sigma_{21}+n_{3}\sigma_{31}\\
  t_{2}=n_{1}\sigma_{12}+n_{2}\sigma_{22}+n_{3}\sigma_{32}\\
  t_{3}=n_{1}\sigma_{13}+n_{2}\sigma_{23}+n_{3}\sigma_{33}\\
  \end{array}  \quad\times\\
  
-  
  
.. math::  
  \int\limits_{\Gamma}{t}_{i}\text{d}{\Gamma}=\int\limits_{\Gamma}{\sigma}_{ij}\cdot{n}_{j}\text{d}{\Gamma}
  =\int\limits_{\Omega}\cfrac{\partial {\sigma}_{ij}}{\partial x_{j}}\text{d}{\Omega} \quad\checkmark 
  
-  
  
.. math::  
  t_{i}=\sigma_{i1}n_{1}+\sigma_{i2}n_{2}+\sigma_{i3}n_{3}  \quad\checkmark \\  
  
-  
  
.. math:: 
  \begin{array}{c}
  t_{1}=\sigma_{11}n_{1}+\sigma_{12}n_{2}+\sigma_{13}n_{3}\\
  t_{2}=\sigma_{21}n_{1}+\sigma_{22}n_{2}+\sigma_{23}n_{3}\\
  t_{3}=\sigma_{31}n_{1}+\sigma_{32}n_{2}+\sigma_{33}n_{3}\\
  \end{array}    \quad\checkmark \\  
  
  
Note that since the normal is to the left on the boundary integral, the divergence is to the left
and contracts with the first index on the stress tensor. When the divergence operator acts on
the first index of the stress tensor it is called the left divergence operator and is placed to the
left. When it acts on the second index, it is placed to the right and called the right divergence.
Since the Cauchy stress is symmetric, the left and right divergence operators have the same
effect. However, in contrast to linear continuum mechanics, in nonlinear continuum mechanics
it is important to become accustomed to placing the divergence operator where it belongs 
because some stress tensors, such as the nominal stress, are not symmetric. When the stress is
not symmetric, the left and right divergence operators lead to different results. In this book we
use the convention that the divergence and gradient operators are placed on the left and the
normal appears on the left in the surface integrals.

.. math::  
  \int\limits_{\Omega}\left(\rho\cfrac{\text{d}\mathbf{v}}{\text{d} t}-\rho\mathbf{b}-\nabla\cdot\boldsymbol\sigma\right)\text{d}{\Omega}=0\quad\times\\  

-  
  
.. math::   
  \int\limits_{\Omega}\left(\rho\cfrac{\text{d}\mathbf{v}}{\text{d} t}-\rho\mathbf{b}-\nabla\cdot([\boldsymbol{\sigma}]^{\text{T}})\right)\text{d}{\Omega}=0\quad\checkmark \\

Therefore, if the integrand is :math:`C_{1}`, the above formula holds for an arbitrary domain, yields

.. math::  
  \rho\cfrac{\text{d}\mathbf{v}}{\text{d} t}=\nabla\cdot\boldsymbol\sigma+\rho\mathbf{b}\equiv \text{div}\boldsymbol\sigma+\rho\mathbf{b} \quad\times\\  

-  

.. math::  
  \rho\cfrac{\text{d}\mathbf{v}}{\text{d} t}=\nabla\cdot([\boldsymbol{\sigma}]^{\text{T}})+\rho\mathbf{b}\equiv \text{div}\boldsymbol\sigma+\rho\mathbf{b} \quad\checkmark \\
  
-  
  
.. math::   
  \rho\cfrac{\text{d}{v}_{i}}{\text{d} t}=\cfrac{\partial \sigma_{ji}}{\partial x_{j}}+\rho b_{i}\quad\times\\  
  
-  
  
.. math::  
  \rho\cfrac{\text{d}{v}_{i}}{\text{d} t}=\cfrac{\partial \sigma_{ij}}{\partial x_{j}}+\rho b_{i}\quad\checkmark \\
  
This is called the momentum equation; it is also called the balance of linear momentum
equation. The LHS term represents the change in momentum, since it is a product of the
acceleration and the density; it is also called the inertial or kinetic term. The first term on
the RHS is the net resultant internal force per unit volume due to the divergence of the
stress field.  

This form of the momentum equation is applicable to both Lagrangian and Eulerian descriptions. In a Lagrangian description, the dependent variables are assumed to be functions of the
Lagrangian coordinates :math:`\mathbf{X}` and time :math:`t`, so the momentum equation is

.. math::  
  \rho(\mathbf{X},t)\cfrac{\partial \mathbf{v}(\mathbf{X},t)}{\partial t}
  =\text{div }\boldsymbol\sigma(\boldsymbol\phi^{-1}(\mathbf{x},t),t)
  +\rho(\mathbf{X},t)\mathbf{b}(\mathbf{X},t)
  
Note that the stress must be expressed as a function of the Eulerian coordinates through the
inverse of the motion :math:`\boldsymbol\phi^{-1}(\mathbf{x},t)` so that the spatial divergence of the stress field can be evaluated,
but it is considered a function of :math:`\mathbf{X}` and time :math:`t`, :math:`\boldsymbol\sigma(\mathbf{X},t)`. The material derivative of the velocity
with respect to time becomes a partial derivative with respect to time when the
independent variables are changed from :math:`\mathbf{x}` to :math:`\mathbf{X}`. 

This would not be considered a true Lagrangian description in classical texts on continuum
mechanics because of the appearance of the derivative with respect to Eulerian coordinates.
However, the essential feature of a Lagrangian description is that the independent variables
are the Lagrangian (material) coordinates. This requirement is met by this, and we will see in
the development of the updated Lagrangian finite element method that this form of the
momentum equation can be discretized with a Lagrangian mesh.

In an Eulerian description, the material derivative of the velocity is written out
and all variables are considered functions of the Eulerian coordinates:

.. math::  
  \rho(\mathbf{x},t)\left(\cfrac{\partial \mathbf{v}(\mathbf{x},t)}{\partial t}+(\mathbf{v}(\mathbf{x},t)\cdot \nabla )\mathbf{v}(\mathbf{x},t)\right )
  =\text{div }\boldsymbol\sigma(\mathbf{x},t)
  +\rho(\mathbf{x},t)\mathbf{b}(\mathbf{x},t)
  
or

.. math::  
  \rho\left(\cfrac{\partial {v}_{i}}{\partial t}+({v}_{j}\cdot\cfrac{\partial}{\partial x_{j}} ){v}_{i}\right )=\cfrac{\partial\sigma_{ij}}{\partial x_{j}} +\rho b_{i}
  
-

.. math::  
  \rho\left(\cfrac{\partial {v}_{i}}{\partial t}+{v}_{j}\cdot\cfrac{\partial {v}_{i}}{\partial x_{j}} \right )=\cfrac{\partial\sigma_{ij}}{\partial x_{j}} +\rho b_{i}\\
  
-
  
.. math:: 
  \begin{align}
  \cfrac{\text{d} \mathbf{v}}{\text{d} t} & = \cfrac{\partial \mathbf{v}}{\partial t}+[(\text{grad }\mathbf{v})][\mathbf{v}]\\
  \cfrac{\text{d} \mathbf{v}}{\text{d} t} & = \cfrac{\partial \mathbf{v}}{\partial t}+[\nabla \mathbf{v}]^{\text{T}}[\mathbf{v}] \\
  \cfrac{\text{d} \mathbf{v}}{\text{d} t} & = \cfrac{\partial \mathbf{v}}{\partial t}+ (\mathbf{v}\cdot \nabla)\mathbf{v} \\
  \cfrac{\text{d} \boldsymbol{\alpha}}{\text{d} t} & = \cfrac{\partial \boldsymbol{\alpha}}{\partial t}+ (\mathbf{v}\cdot \nabla)\boldsymbol{\alpha} \\
  \end{align}  
  
-
  
.. math:: 
  \cfrac{\text{d}(f(\mathbf{x},t)g(\mathbf{x},t))}{\text{d} t}
  =g(\mathbf{x},t)\cfrac{\text{d}f(\mathbf{x},t)}{\text{d} t}+f(\mathbf{x},t)\cfrac{\text{d}g(\mathbf{x},t)}{\text{d} t}  
  
-
  
.. math:: 
  \cfrac{\text{d}f(\mathbf{x},t)}{\text{d} t}=\cfrac{\partial f(\mathbf{x},t)}{\partial  t}+v_{i}\cfrac{\partial f(\mathbf{x},t)}{\partial  x_{i}}

-
  
.. math:: 
  \cfrac{\text{d}(f\cdot g)}{\text{d} t}=\cfrac{\partial(f\cdot g)}{\partial  t}+v_{i}\cfrac{\partial (f\cdot g)}{\partial  x_{i}}  

Reynolds’ transport theorem  

.. math:: 
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}f\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\text{d}(f)}{\text{d} t}+(f)\cfrac{\partial v_{i}}{\partial  x_{i}}\right)\text{d}\Omega  
  
-
  
.. math:: 
  \begin{array}{c}
  \displaystyle  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(f\cdot g_{1})\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\text{d}(f\cdot g_{1})}{\text{d} t}+(f\cdot g_{1})\cfrac{\partial v_{i}}{\partial  x_{i}}\right)\text{d}\Omega\\
  \displaystyle  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(f\cdot g_{2})\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\text{d}(f\cdot g_{2})}{\text{d} t}+(f\cdot g_{2})\cfrac{\partial v_{i}}{\partial  x_{i}}\right)\text{d}\Omega\\
  \displaystyle  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(f\cdot g_{3})\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\text{d}(f\cdot g_{3})}{\text{d} t}+(f\cdot g_{3})\cfrac{\partial v_{i}}{\partial  x_{i}}\right)\text{d}\Omega
  \end{array}  

-
  
.. math:: 
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(f\cdot \mathbf{g})\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\text{d}(f\cdot \mathbf{g})}{\text{d} t}+(f\cdot \mathbf{g})\cfrac{\partial v_{i}}{\partial  x_{i}}\right)\text{d}\Omega  
  
Partial derivative

.. math:: 
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(f)\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\partial(f)}{\partial t}+v_{i}\cfrac{\partial (f)}{\partial  x_{i}}+(f)\cfrac{\partial v_{i}}{\partial  x_{i}}\right)\text{d}\Omega
  =\int\limits_{\Omega}\left(\cfrac{\partial(f)}{\partial t}+\cfrac{\partial (v_{i}\cdot f)}{\partial  x_{i}}\right)\text{d}\Omega
 
-
 
.. math:: 
  \begin{array}{c}
  \displaystyle  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(f\cdot g_{1})\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\partial(f\cdot g_{1})}{\partial t}+\cfrac{\partial (v_{i}\cdot f\cdot g_{1})}{\partial  x_{i}}\right)\text{d}\Omega\\
  \displaystyle  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(f\cdot g_{2})\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\partial(f\cdot g_{2})}{\partial t}+\cfrac{\partial (v_{i}\cdot f\cdot g_{2})}{\partial  x_{i}}\right)\text{d}\Omega\\
  \displaystyle  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(f\cdot g_{3})\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\partial(f\cdot g_{3})}{\partial t}+\cfrac{\partial (v_{i}\cdot f\cdot g_{3})}{\partial  x_{i}}\right)\text{d}\Omega\\
  \end{array}  
  
-
 
.. math:: 
  \displaystyle  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(f\cdot \mathbf{g})\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\partial(f\cdot \mathbf{g})}{\partial t}+\cfrac{\partial (v_{i}\cdot f\cdot \mathbf{g})}{\partial  x_{i}}\right)\text{d}\Omega\\

-
 
.. math:: 
  \cfrac{\text{d}(f\cdot \mathbf{g})}{\text{d} t}+(f\cdot \mathbf{g})\cfrac{\partial v_{i}}{\partial  x_{i}}=
  \cfrac{\partial(f\cdot \mathbf{g})}{\partial t}+\cfrac{\partial (v_{i}\cdot f\cdot \mathbf{g})}{\partial  x_{i}}  
  
-
 
.. math::
  \cfrac{\text{d}(f\cdot \mathbf{g})}{\text{d} t}+(f\cdot \mathbf{g})\cfrac{\partial v_{i}}{\partial  x_{i}}=
  \cfrac{\partial(f\cdot \mathbf{g})}{\partial t}+\cfrac{\partial (v_{i}\cdot f\cdot \mathbf{g})}{\partial  x_{i}}
  =\cfrac{\partial(f\cdot \mathbf{g})}{\partial t}+( f\cdot \mathbf{g})\cfrac{\partial (v_{i})}{\partial  x_{i}}
  +(v_{i})\cfrac{\partial ( f\cdot \mathbf{g})}{\partial  x_{i}}  
  
-
 
.. math:: 
  \cfrac{\text{d}(f\cdot \mathbf{g})}{\text{d} t}
  =\cfrac{\partial(f\cdot \mathbf{g})}{\partial t}+v_{i}\cfrac{\partial ( f\cdot \mathbf{g})}{\partial  x_{i}}  
  
-
 
.. math:: 
  \cfrac{\text{d}(f\cdot \mathbf{g})}{\text{d} t}
  =\cfrac{\partial(f\cdot \mathbf{g})}{\partial t}
  +v_{1}\cfrac{\partial ( f\cdot \mathbf{g})}{\partial  x_{1}}
  +v_{2}\cfrac{\partial ( f\cdot \mathbf{g})}{\partial  x_{2}}
  +v_{3}\cfrac{\partial ( f\cdot \mathbf{g})}{\partial  x_{3}}\\  

-
 
.. math:: 
  \cfrac{\text{d}(f\cdot \mathbf{g})}{\text{d} t}
  =\cfrac{\partial(f\cdot \mathbf{g})}{\partial t}+(\mathbf{v}\cdot\nabla ) ( f\cdot \mathbf{g})
  
-
 
.. math:: 
  \cfrac{\text{d}(f\cdot \mathbf{g})}{\text{d} t}
  =\cfrac{\partial(f\cdot \mathbf{g})}{\partial t}+(\mathbf{v}\cdot\nabla )\begin{bmatrix}
   ( f\cdot {g}_{1})\\( f\cdot {g}_{2})\\( f\cdot {g}_{3})
  \end{bmatrix}
  =\cfrac{\partial(f\cdot \mathbf{g})}{\partial t}+\begin{bmatrix}
   (\mathbf{v}\cdot\nabla )( f\cdot {g}_{1})\\(\mathbf{v}\cdot\nabla )( f\cdot {g}_{2})\\(\mathbf{v}\cdot\nabla )( f\cdot {g}_{3})
  \end{bmatrix}  
  
For scalar

.. math:: 
  (\mathbf{v}\cdot\nabla)(f)=(\mathbf{v})\cdot(\nabla f)=v_{i}\cfrac{\partial f}{\partial x_{i}}
  =v_{1}\cfrac{\partial f}{\partial x_{1}}
  +v_{2}\cfrac{\partial f}{\partial x_{2}}
  +v_{3}\cfrac{\partial f}{\partial x_{3}} 

For vector 

.. math:: 
  \cfrac{\text{d}(\mathbf{g})}{\text{d} t}
  =\cfrac{\partial(\mathbf{g})}{\partial t}+(\mathbf{v}\cdot\nabla ) ( \mathbf{g}) 
  
- `Vector calculus <https://en.wikipedia.org/wiki/Vector_calculus_identities>`_

.. math:: 
  \nabla \mathbf{a}=\left[\cfrac{\partial a_{j}}{\partial x_{i}}\right]_{ij}
  =\begin{bmatrix}
  \cfrac{\partial a_{1}}{\partial x_{1}}& \cfrac{\partial a_{2}}{\partial x_{1}} & \cfrac{\partial a_{3}}{\partial x_{1}}\\
  \cfrac{\partial a_{1}}{\partial x_{2}}& \cfrac{\partial a_{2}}{\partial x_{2}} & \cfrac{\partial a_{3}}{\partial x_{2}}\\
  \cfrac{\partial a_{1}}{\partial x_{3}}& \cfrac{\partial a_{2}}{\partial x_{3}} & \cfrac{\partial a_{3}}{\partial x_{3}}\\
  \end{bmatrix}
  
- 
 
.. math:: 
  \text{grad } \mathbf{a}=\left[\cfrac{\partial a_{i}}{\partial x_{j}}\right]_{ij}
  =\begin{bmatrix}
  \cfrac{\partial a_{1}}{\partial x_{1}}& \cfrac{\partial a_{1}}{\partial x_{2}} & \cfrac{\partial a_{1}}{\partial x_{3}}\\
  \cfrac{\partial a_{2}}{\partial x_{1}}& \cfrac{\partial a_{2}}{\partial x_{2}} & \cfrac{\partial a_{2}}{\partial x_{3}}\\
  \cfrac{\partial a_{3}}{\partial x_{1}}& \cfrac{\partial a_{3}}{\partial x_{2}} & \cfrac{\partial a_{3}}{\partial x_{3}}\\
  \end{bmatrix}=(\nabla \mathbf{a})^{\text{T}} 
  
- 
 
.. math:: 
  [\text{grad } \mathbf{a}]\cdot  \mathbf{b}
  =\begin{bmatrix}
  \cfrac{\partial a_{1}}{\partial x_{1}}& \cfrac{\partial a_{1}}{\partial x_{2}} & \cfrac{\partial a_{1}}{\partial x_{3}}\\
  \cfrac{\partial a_{2}}{\partial x_{1}}& \cfrac{\partial a_{2}}{\partial x_{2}} & \cfrac{\partial a_{2}}{\partial x_{3}}\\
  \cfrac{\partial a_{3}}{\partial x_{1}}& \cfrac{\partial a_{3}}{\partial x_{2}} & \cfrac{\partial a_{3}}{\partial x_{3}}\\
  \end{bmatrix}
  \begin{bmatrix}
  b_{1}\\b_{2}\\b_{3}\\
  \end{bmatrix}
  =\begin{bmatrix}
  b_{1}\cfrac{\partial a_{1}}{\partial x_{1}}
  +b_{2}\cfrac{\partial a_{1}}{\partial x_{2}}
  +b_{3}\cfrac{\partial a_{1}}{\partial x_{3}}\\
   b_{1}\cfrac{\partial a_{2}}{\partial x_{1}}
  +b_{2}\cfrac{\partial a_{2}}{\partial x_{2}}
  +b_{3}\cfrac{\partial a_{2}}{\partial x_{3}}\\
   b_{1}\cfrac{\partial a_{3}}{\partial x_{1}}
  +b_{2}\cfrac{\partial a_{3}}{\partial x_{2}}
  +b_{3}\cfrac{\partial a_{3}}{\partial x_{3}}\\
  \end{bmatrix}
  
- 
 
.. math:: 
  [\text{grad }(\rho \mathbf{v})]\cdot  \mathbf{v}
  =\begin{bmatrix}
  \cfrac{\partial (\rho v_{1})}{\partial x_{1}}& \cfrac{\partial (\rho v_{1})}{\partial x_{2}} & \cfrac{\partial (\rho v_{1})}{\partial x_{3}}\\
  \cfrac{\partial (\rho v_{2})}{\partial x_{1}}& \cfrac{\partial (\rho v_{2})}{\partial x_{2}} & \cfrac{\partial (\rho v_{2})}{\partial x_{3}}\\
  \cfrac{\partial (\rho v_{3})}{\partial x_{1}}& \cfrac{\partial (\rho v_{3})}{\partial x_{2}} & \cfrac{\partial (\rho v_{3})}{\partial x_{3}}\\
  \end{bmatrix}
  \begin{bmatrix}
  v_{1}\\v_{2}\\v_{3}\\
  \end{bmatrix}
  =\begin{bmatrix}
   v_{1}\cfrac{\partial (\rho v_{1})}{\partial x_{1}}
  +v_{2}\cfrac{\partial (\rho v_{1})}{\partial x_{2}}
  +v_{3}\cfrac{\partial (\rho v_{1})}{\partial x_{3}}\\
   v_{1}\cfrac{\partial (\rho v_{2})}{\partial x_{1}}
  +v_{2}\cfrac{\partial (\rho v_{2})}{\partial x_{2}}
  +v_{3}\cfrac{\partial (\rho v_{2})}{\partial x_{3}}\\
   v_{1}\cfrac{\partial (\rho v_{3})}{\partial x_{1}}
  +v_{2}\cfrac{\partial (\rho v_{3})}{\partial x_{2}}
  +v_{3}\cfrac{\partial (\rho v_{3})}{\partial x_{3}}\\
  \end{bmatrix}  
  
- 
 
.. math:: 
  \begin{align}
  (\mathbf{v}\cdot\nabla ) ( \rho\mathbf{v})
  & = v_{1}\cfrac{\partial( \rho\mathbf{v})}{\partial x_{1}}
  +v_{2}\cfrac{\partial( \rho\mathbf{v})}{\partial x_{2}}
  +v_{3}\cfrac{\partial( \rho\mathbf{v})}{\partial x_{3}}\\
  &=\begin{bmatrix}
  \displaystyle v_{1}\cfrac{\partial( \rho{v}_{1})}{\partial x_{1}}
  +v_{2}\cfrac{\partial( \rho{v}_{1})}{\partial x_{2}}
  +v_{3}\cfrac{\partial( \rho{v}_{1})}{\partial x_{3}}\\
  \displaystyle v_{1}\cfrac{\partial( \rho{v}_{2})}{\partial x_{1}}
  +v_{2}\cfrac{\partial( \rho{v}_{2})}{\partial x_{2}}
  +v_{3}\cfrac{\partial( \rho{v}_{2})}{\partial x_{3}}\\
  \displaystyle v_{1}\cfrac{\partial( \rho{v}_{3})}{\partial x_{1}}
  +v_{2}\cfrac{\partial( \rho{v}_{3})}{\partial x_{2}}
  +v_{3}\cfrac{\partial( \rho{v}_{3})}{\partial x_{3}}\\
  \end{bmatrix}
  \end{align}
  
- 
 
.. math::   
  [\text{grad } \mathbf{g}]\cdot  \mathbf{v}=(\mathbf{v}\cdot\nabla ) ( \mathbf{g}) =v_{j}\cfrac{\partial g_{i}}{\partial x_{j}}  

- 
 
.. math::   
  \begin{align}
  \cfrac{\text{d}(f\cdot\mathbf{g})}{\text{d} t} & = \cfrac{\partial(f\cdot\mathbf{g})}{\partial t}+(\mathbf{v}\cdot\nabla ) ( f\cdot\mathbf{g}) \\ & = \cfrac{\partial(f\cdot\mathbf{g})}{\partial t}+[\text{grad } (f\cdot\mathbf{g})]\cdot  \mathbf{v}
  \end{align}
  
- 
 
.. math::
  \begin{align}
  \cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t} & = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+(\mathbf{v}\cdot\nabla ) ( \rho\mathbf{v}) \\ & = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+[\text{grad } (\rho\mathbf{v})]\cdot  \mathbf{v}
  \end{align}  
  
- 
 
.. math::
  \cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t}  = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+(\mathbf{v}\cdot\nabla ) ( \rho\mathbf{v}) \quad\checkmark\\  
  
- 
 
.. math::
  \cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t} = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+[\text{grad } (\rho\mathbf{v})]\cdot  \mathbf{v} \quad\checkmark\\  
  
- 
 
.. math::
  \cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t} = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+\mathbf{v}\cdot[\text{grad } (\rho\mathbf{v})] \quad\times\\  
  
This is called the conservative form of the momentum equation. In the conservative form, the
momentum per unit volume :math:`\rho\mathbf{v}` is a dependent variable. This variant of the momentum equation
is said to observe momentum conservation more accurately  

.. math::
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}(\rho\mathbf{v})\text{d}\Omega=\int\limits_{\Omega}(\rho\cfrac{\text{d}\mathbf{v}}{\text{d}t})\text{d}\Omega
  
but

.. math::
  \cfrac{\text{d}(\rho\mathbf{v})}{\text{d}t}=\rho\cfrac{\text{d}\mathbf{v}}{\text{d}t}+\mathbf{v}\cfrac{\text{d}\rho}{\text{d}t}\ne \rho\cfrac{\text{d}\mathbf{v}}{\text{d}t}  

Divergence
  
.. math::  
  \text{div}(\mathbf{a}\otimes \mathbf{b})=[\text{grad}(\mathbf{a})]\cdot\mathbf{b}+\mathbf{a}\text{ div}(\mathbf{b})  
  
Let :math:`\mathbf{a}=\rho\mathbf{v}`, :math:`\mathbf{v}=\mathbf{v}`, then

.. math::  
  \text{div}((\rho\mathbf{v})\otimes \mathbf{v})=[\text{grad}(\rho\mathbf{v})]\cdot\mathbf{v}+(\rho\mathbf{v})\text{ div}(\mathbf{v})\\
  
-

.. math::  
  \begin{align}
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}\rho\mathbf{v}\text{d}{\Omega} & = \int\limits_{\Omega}\left(\cfrac{\text{d}}{\text{d} t}(\rho\mathbf{v})+(\text{div }\mathbf{v})(\rho\mathbf{v})\right)\text{d}{\Omega}\\
  & = \int\limits_{\Omega}\left(\cfrac{\partial(\rho\mathbf{v})}{\partial t}+[\text{grad } (\rho\mathbf{v})]\cdot  \mathbf{v}+(\text{div }\mathbf{v})(\rho\mathbf{v})\right)\text{d}{\Omega}\\
  & = \int\limits_{\Omega}\left(\cfrac{\partial(\rho\mathbf{v})}{\partial t}+\text{div}((\rho\mathbf{v})\otimes \mathbf{v})\right)\text{d}{\Omega}\\
  & = \int\limits_{\Omega}\left(\cfrac{\partial(\rho\mathbf{v})}{\partial t}+\text{div}(\rho\mathbf{v}\otimes \mathbf{v})\right)\text{d}{\Omega}\\
  \end{align}
  
-

.. math::  
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}\rho\mathbf{v}\text{d}{\Omega}
  = \int\limits_{\Omega}\left(\cfrac{\partial(\rho\mathbf{v})}{\partial t}+\text{div}(\rho\mathbf{v}\otimes \mathbf{v})\right)\text{d}{\Omega}\\
  
-

.. math::
  \cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t}\ne\cfrac{\partial(\rho\mathbf{v})}{\partial t}+\text{div}(\rho\mathbf{v}\otimes \mathbf{v})  
  
-

.. math::
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}\rho\mathbf{v}\text{d}{\Omega}
  = \int\limits_{\Omega}\left(\cfrac{\partial(\rho\mathbf{v})}{\partial t}+\text{div}(\rho\mathbf{v}\otimes \mathbf{v})\right)\text{d}{\Omega}
  =\int\limits_{\Omega}\left(\text{div }\boldsymbol\sigma +\rho \mathbf{b}\right)\text{d}{\Omega}\\  
  
-

.. math::
  \cfrac{\partial(\rho\mathbf{v})}{\partial t}+\text{div}(\rho\mathbf{v}\otimes \mathbf{v})=\text{div }\boldsymbol\sigma +\rho \mathbf{b}  
  
-

.. math:: 
  \begin{align}
  \displaystyle\cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t} & = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+(\mathbf{v}\cdot\nabla ) ( \rho\mathbf{v}) \quad\checkmark\\
  \displaystyle\cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t} & = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+[\text{grad } (\rho\mathbf{v})]\cdot  \mathbf{v} \quad\checkmark\\
  \displaystyle\cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t} & = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+\mathbf{v}\cdot[\text{grad } (\rho\mathbf{v})] \quad\times\\
  \end{align}  
  
-

.. math::  
  \begin{align}
  \displaystyle\cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t}& = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+(\mathbf{v}\cdot\nabla ) ( \rho\mathbf{v})=\text{div }\boldsymbol\sigma +\rho \mathbf{b} \quad\times\\
  \displaystyle\cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t}& = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+[\text{grad } (\rho\mathbf{v})]\cdot  \mathbf{v}=\text{div }\boldsymbol\sigma +\rho \mathbf{b} \quad\times\\
  \displaystyle\cfrac{\text{d}(\rho\mathbf{v})}{\text{d} t}& = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+\mathbf{v}\cdot[\text{grad } (\rho\mathbf{v})]=\text{div }\boldsymbol\sigma +\rho \mathbf{b} \quad\times\\
  \end{align}  
  
-

.. math::  
  \begin{align}
  \cfrac{\text{d}}{\text{d} t}(\rho\mathbf{v})+(\text{div }\mathbf{v})(\rho\mathbf{v}) &= \cfrac{\partial(\rho\mathbf{v})}{\partial t}+[\text{grad } (\rho\mathbf{v})]\cdot  \mathbf{v}+(\text{div }\mathbf{v})(\rho\mathbf{v})\\
   & = \cfrac{\partial(\rho\mathbf{v})}{\partial t}+\text{div }(\rho\mathbf{v}\otimes\mathbf{v})\\
   & = \text{div }\boldsymbol\sigma +\rho \mathbf{b} \quad\checkmark\\  
  \end{align}

  
Reynolds’ Theorem for a Density-Weighted Integrand
------------------------------------------------------
The material time derivative of an integral in which the integrand is a product of the density and a function f is given by

.. math::  
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}(\rho f)\text{d}\Omega=\int\limits_{\Omega}(\rho\cfrac{\text{d}f}{\text{d}t})\text{d}\Omega

This holds for a tensor of any order and is a consequence of Reynolds’ theorem and mass
conservation; it is another variant of Reynolds’ theorem.

Conservation of Energy
------------------------------------------------------
We consider thermomechanical processes where the only sources of energy are mechanical
work and heat. The principle of conservation of energy, that is, the energy balance principle,
states that the rate of change of total energy is equal to the work done by the body forces and
surface tractions plus the heat energy delivered to the body by the heat flux and other sources
of heat. The internal energy per unit volume is denoted by :math:`\rho E_{int}` where :math:`E_{int}` is the internal
energy per unit mass. The heat flux per unit area is denoted by a vector :math:`\mathbf{q}`, in units of power
per area, and the heat source per unit volume is denoted by :math:`\rho s`. The conservation of energy
then requires that the rate of change of the total energy in the body, which includes both
internal energy and kinetic energy, equals the power of the applied forces and the energy
added to the body by heat conduction and any heat sources.

The rate of change of the total energy in the body is given by

.. math:: 
  \begin{align}
  \displaystyle P^{tot} & = P^{int}+P^{kin}\\
  \displaystyle P^{int} & = \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}(\rho E^{int})\text{d}\Omega\\
  \displaystyle P^{kin} & = \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}(\cfrac{1}{2}\rho \mathbf{v}\cdot\mathbf{v})\text{d}\Omega\\
  \end{align}
  
where :math:`P^{int}` denotes the rate of change of internal energy and :math:`P^{kin}` the rate of change of the
kinetic energy. The rate of work by the body forces in the domain and the tractions on the
surface is  

.. math:: 
  \begin{align}
  \displaystyle P^{ext}  & = \int\limits_{\Omega}(\mathbf{v}\cdot\rho \mathbf{b})\text{d}\Omega+
  \int\limits_{\Gamma}(\mathbf{v}\cdot \mathbf{t})\text{d}\Gamma\\
  & = \int\limits_{\Omega}({v}_{i}\cdot\rho {b}_{i})\text{d}\Omega+
  \int\limits_{\Gamma}({v}_{i}\cdot {t}_{i})\text{d}\Gamma\\\\
  \end{align}
  
The power supplied by heat sources :math:`s` and the heat flux :math:`\mathbf{q}` is  

.. math:: 
  \begin{align}
  \displaystyle P^{heat}  & = \int\limits_{\Omega}(\rho s)\text{d}\Omega
  -\int\limits_{\Gamma}(\mathbf{n}\cdot \mathbf{q})\text{d}\Gamma\\
   & = \int\limits_{\Omega}(\rho s)\text{d}\Omega
  -\int\limits_{\Gamma}({n}_{i}\cdot {q}_{i})\text{d}\Gamma\\
  \end{align}
  
where the sign of the heat flux term is negative since positive heat flow is out of the body.

The statement of the conservation of energy is 

.. math:: 
  P^{tot} = P^{ext}+P^{heat}
  
that is, the rate of change of the total energy in the body (consisting of the internal and kinetic
energies) is equal to the rate of work by the external forces and rate of work provided by heat flux
and energy sources. This is known as the first law of thermodynamics. The disposition of the
internal energy depends on the material. In an elastic material, it is stored as elastic internal energy
and is fully recoverable upon unloading. In an elastic–plastic material, some of the internal energy
is converted to heat and some is dissipated in changes of the internal structure of the material.

.. math:: 
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}(\rho E^{int}+\cfrac{1}{2}\rho \mathbf{v}\cdot\mathbf{v})\text{d}\Omega
  =\int\limits_{\Omega}(\mathbf{v}\cdot\rho \mathbf{b})\text{d}\Omega+
  \int\limits_{\Gamma}(\mathbf{v}\cdot \mathbf{t})\text{d}\Gamma
  +\int\limits_{\Omega}(\rho s)\text{d}\Omega
  -\int\limits_{\Gamma}(\mathbf{n}\cdot \mathbf{q})\text{d}\Gamma
  
The Divergence Theorem-vector 

.. math::
  \int\limits_{S}\mathbf{v}\cdot\mathbf{n}\text{d}S=\int\limits_{V}\text{div }\mathbf{v}\text{d}V

-  

.. math::
  \int\limits_{S}{v}_{i}\cdot{n}_{i}\text{d}S=\int\limits_{V}\cfrac{\partial {v}_{i}}{\partial {x}_{i}}\text{d}V\\


Applying Cauchy’s law and Gauss’s theorem to the traction boundary integrals yields:

.. math:: 
  \int\limits_{\Gamma} n_{i} \sigma_{i j} v_{j} d \Gamma
  =\int\limits_{\Omega}\left(\cfrac{\partial (\sigma_{i j} v_{j})}{\partial x_{i}}\right)d \Omega \\
  =\int\limits_{\Omega}\left(\sigma_{i j}\cfrac{\partial ( v_{j})}{\partial x_{i}}+v_{j}\cfrac{\partial ( \sigma_{i j})}{\partial x_{i}}\right)d \Omega \\

-
  
.. math::   
  \int\limits_{\Gamma} n_{i} \sigma_{i j} v_{j} d \Gamma
  =\int\limits_{\Omega}\left(\cfrac{\partial (\sigma_{i j} v_{j})}{\partial x_{i}}\right)d \Omega \\
  =\int\limits_{\Omega}\left(\text{div }(\boldsymbol{\sigma}\cdot\mathbf{v})\right)d \Omega \\
  
then
  
.. math:: 
  \int\limits_{\Gamma} \mathbf{v} \cdot \mathbf{t} d \Gamma  =\int\limits_{\Gamma} \mathbf{n} \cdot \boldsymbol{\sigma} \cdot \mathbf{v} d \Gamma=\int\limits_{\Gamma} n_{i} \sigma_{i j} v_{j} d \Gamma
  =\int\limits_{\Omega}\left(\cfrac{\partial (\sigma_{i j} v_{j})}{\partial x_{i}}\right)d \Omega
  =\int\limits_{\Omega}\left(\text{div }(\boldsymbol{\sigma}\cdot\mathbf{v})\right)d \Omega \\
  
-

.. math:: 
  \int\limits_{\Gamma}(\mathbf{n}\cdot \mathbf{q})\text{d}\Gamma
  =\int\limits_{\Omega}\text{div }\mathbf{ q}\text{d}\Omega
  =\int\limits_{\Omega}\nabla\cdot\mathbf{ q}\text{d}\Omega
  =\int\limits_{\Omega}\cfrac{\partial q_{i}}{\partial x_{i}}\text{d}\Omega

-

.. math:: 
  \int\limits_{\Gamma} \mathbf{v} \cdot \mathbf{t} d \Gamma-\int\limits_{\Gamma}(\mathbf{n}\cdot \mathbf{q})\text{d}\Gamma
  = \int\limits_{\Omega}\left(\text{div }(\boldsymbol{\sigma}\cdot\mathbf{v})-\text{div }\mathbf{ q}\right)d \Omega\\
  
-

.. math::
  \rho E=\rho E^{int}+\cfrac{1}{2}\rho \mathbf{v}\cdot\mathbf{v}  

-
  
.. math:: 
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}f\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\text{d}(f)}{\text{d} t}+(f)\cfrac{\partial v_{i}}{\partial  x_{i}}\right)\text{d}\Omega  
  
-
  
.. math:: 
  \cfrac{\text{d}}{\text{d} t}\int\limits_{\Omega}(\rho E)\text{d}\Omega=\int\limits_{\Omega}\left(\cfrac{\text{d}(\rho E)}{\text{d} t}+(\rho E)\cfrac{\partial v_{i}}{\partial  x_{i}}\right)\text{d}\Omega  

-

.. math:: 
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}(\rho E)\text{d}\Omega
  =\int\limits_{\Omega}(\mathbf{v}\cdot\rho \mathbf{b}+\rho s)\text{d}\Omega+
  \int\limits_{\Gamma}(\mathbf{v}\cdot \mathbf{t})\text{d}\Gamma
  -\int\limits_{\Gamma}(\mathbf{n}\cdot \mathbf{q})\text{d}\Gamma

-

.. math:: 
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}(\rho E)\text{d}\Omega
  =\int\limits_{\Omega}(\mathbf{v}\cdot\rho \mathbf{b}+\rho s+\text{div }(\boldsymbol{\sigma}\cdot\mathbf{v})-\text{div }\mathbf{ q})\text{d}\Omega
  
-

.. math:: 
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}f\text{d}\Omega 
  = \int\limits_{\Omega}\left(\cfrac{\partial f}{\partial t}+\text{div } \left(\mathbf{v}f\right)\right)\text{d}\Omega
  
-

.. math::   
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}(\rho E)\text{d}\Omega 
  = \int\limits_{\Omega}\left(\cfrac{\partial (\rho E)}{\partial t}+\text{div } \left((\rho E)\mathbf{v}\right)\right)\text{d}\Omega
  
-  
  
.. math:: 
  \begin{align}
  \cfrac{\text{d}}{\text{d}t}\int\limits_{\Omega}(\rho E)\text{d}\Omega
  & = \int\limits_{\Omega}\left(\cfrac{\partial (\rho E)}{\partial t}+\text{div } \left((\rho E)\mathbf{v}\right)\right)\text{d}\Omega \\
  & = \int\limits_{\Omega}(\text{div }(\boldsymbol{\sigma}\cdot\mathbf{v})-\text{div }\mathbf{ q}+\mathbf{v}\cdot\rho \mathbf{b}+\rho s)\text{d}\Omega\\  
  \end{align} 
  
-  
  
.. math::  
  \begin{array}{c}
  \text{div }(\mathbf{T}\cdot\mathbf{a})
  = \text{tr}(\mathbf{T}\text{ grad } \mathbf{a})+\mathbf{a}\cdot \text{div}(\mathbf{T}^{\text{T}})\\
  \text{div }(\boldsymbol{\sigma}\cdot\mathbf{v})
  = \text{tr}(\boldsymbol{\sigma}\text{ grad } \mathbf{v})+\mathbf{v}\cdot \text{div}(\boldsymbol{\sigma}^{\text{T}})\\
  \end{array}  
  
-  
  
.. math::  
  \cfrac{\partial (\rho E)}{\partial t}+\text{div } \left((\rho E)\mathbf{v}\right)=\text{div }(\boldsymbol{\sigma}\cdot\mathbf{v})-\text{div }\mathbf{ q}+\mathbf{v}\cdot\rho \mathbf{b}+\rho s  
  
-  
  
.. math::    
  \boldsymbol\sigma=\boldsymbol\tau-p\mathbf{I}  
  
-  
  
.. math::
  \text{div }(\boldsymbol{\sigma}\cdot\mathbf{v})=\text{div }(\boldsymbol{\tau}\cdot\mathbf{v}-p\mathbf{v})\\
  
-  
  
.. math::
  \cfrac{\partial (\rho E)}{\partial t}+\text{div } \left((\rho E+p)\mathbf{v}\right)=\text{div }(\boldsymbol{\tau}\cdot\mathbf{v})-\text{div }\mathbf{ q}+\mathbf{v}\cdot\rho \mathbf{b}+\rho s  
  
-  
  
.. math::
  \cfrac{\partial (\rho E)}{\partial t}+\text{div } \left(\rho H\mathbf{v}\right)=\text{div }(\boldsymbol{\tau}\cdot\mathbf{v})-\text{div }\mathbf{ q}+\mathbf{v}\cdot\rho \mathbf{b}+\rho s  
  
where

.. math::
  \rho H=\rho E+p; \quad H=E+p/\rho

  
  

