Navier–Stokes equations
==================================

In physics, the Navier–Stokes equations (/nævˈjeɪ stoʊks/ nav-YAY STOHKS) are partial differential equations which describe the motion of viscous fluid substances, named after French engineer and physicist Claude-Louis Navier and Anglo-Irish physicist and mathematician George Gabriel Stokes. They were developed over several decades of progressively building the theories, from 1822 (Navier) to 1842-1850 (Stokes).

Control Mass(CM) &Control Volume (CV)
--------------------------------------

Conservation laws can be derived by considering a given quantity of matter or 
control mass (CM) and its extensive properties, such as mass, momentum and 
energy. This approach is used to study the dynamics of solid bodies, where the 
CM (sometimes called the system) is easily identified. In fluid flows, however, 
it is difficult to follow a parcel of matter. It is more convenient to deal with 
the flow within a certain spatial region we call a control volume (CV), rather 
than in a parcel of matter which quickly passes through the region of interest. 
This method of analysis is called the control volume approach.  

Examples are density :math:`\rho` (mass per unit volume) 
and velocity :math:`v \mathbf{v} \boldsymbol{v}` (momentum per unit mass).

If :math:`\phi` is any conserved intensive property (for mass conservation, :math:`\phi=1`;  for 
momentum conservation, :math:`\phi=\boldsymbol{v}`;  for conservation of a scalar, :math:`\phi` represents the 
conserved property per unit mass), then the corresponding extensive property 
:math:`\Phi` can be expressed as:

.. math::    
  \Phi =\int_{\Omega _{CM} } \rho \phi d\Omega
  
where :math:`\Omega_{CM}` stands for volume occupied by the CM.  
  
Reynolds' transport theorem
------------------------------

.. math::

  \cfrac{\mathrm{d}}{\mathrm{d}t}\int_{\Omega(t)} \mathbf{f}~\text{dV} = 
  \int_{\Omega(t)} \frac{\partial \mathbf{f}}{\partial t}~\text{dV} + \int_{\partial \Omega(t)} (\mathbf{v}_{\mathrm{b}}\cdot\mathbf{n})\mathbf{f}~\text{dA}~ ~
  
.. math::  
  \frac{\mathrm{d}}{\mathrm{d} t} \int_{\Omega_{\mathrm{CM}}} \rho \phi \mathrm{d} \Omega=
  \frac{\mathrm{d}}{\mathrm{d} t} \int_{\Omega_{\mathrm{CV}}} \rho \phi \mathrm{d} \Omega+
  \int_{S_{\mathrm{CV}}} \rho \phi\left(\boldsymbol{v}-\boldsymbol{v}_{\mathrm{b}}\right) \cdot n \mathrm{~d} S
  
where :math:`\Omega_{\mathrm{CV}}` is the CV volume, :math:`S_{\mathrm{CV}}` is the surface enclosing CV, :math:`\boldsymbol{n}` is the unit 
vector orthogonal to :math:`S_{\mathrm{CV}}` and directed outwards, :math:`\boldsymbol{v}` is the fluid velocity and :math:`\boldsymbol{v}_{\mathrm{b}}` 
is the velocity with which the CV surface is moving. For a fixed CV, which 
we shall be considering most of the time, :math:`\boldsymbol{v}_{\mathrm{b}}=0` and the first derivative 
on the right hand side becomes a local (partial) derivative. This equation 
states that the rate of change of the amount of the property in the control 
mass, :math:`\Phi`, is the rate of change of the property within the control volume plus 
the net flux of it through the CV boundary due to fluid motion relative to 
CV boundary. The last term is usually called the convective (or sometimes, 
advective) flux of :math:`\phi` through the CV boundary. If the CV moves so that its 
boundary coincides with the boundary of a control mass, then :math:`\boldsymbol{v}=\boldsymbol{v}_{\mathrm{b}}`  and 
this term will be zero as required. 
A detailed derivation of this equation is given in in many textbooks on 
fluid dynamics (e.g. in Bird et al., 1962; Fox and McDonald, 1982) and will not 
be repeated here. The mass, momentum and scalar conservation equations 
will be presented in the next three sections. For convenience, a fixed CV will 
be considered; :math:`\Omega` represents the CV volume and S its surface.  


Conservation Equations
------------------------------

Conservation Laws
`````````````````````````````
One group of the fundamental equations of continuum mechanics arises from the conservation
laws. These equations must always be satisfied by physical systems. Four conservation laws
relevant to thermomechanical systems are considered here:

#. Conservation of mass
#. Conservation of linear momentum, often called conservation of momentum
#. Conservation of energy
#. Conservation of angular momentum.

The conservation laws are also known as balance laws, for example, the conservation of
energy is often called the balance of energy.

Material Time Derivative of an Integral and Reynolds’ Transport Theorem
```````````````````````````````````````````````````````````````````````````````````````
The material time derivative of an integral is the rate of change of an integral on a material
domain. A material domain moves with the material, so that the material points on the
boundary remain on the boundary and no mass flux occurs across the boundaries. A material
domain is analogous to a Lagrangian mesh; a Lagrangian element or group of Lagrangian elements is a nice example of a material domain. The various forms for material time derivatives
of integrals are called Reynolds’ transport theorem.

The material time derivative of an integral is defined by

.. math::
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }^{}f(\mathbf{x},t)\mathrm{d}\Omega 
  =\lim_{\Delta  t \to 0}\cfrac{1}{\Delta  t }\left ( \int\limits_{\Omega_{\tau+\Delta t} }f(\mathbf{x},{\tau+\Delta t})\mathrm{d}\Omega-\int\limits_{\Omega_{\tau} }f(\mathbf{x},\tau )\mathrm{d}\Omega \right )
  
where :math:`\Omega_{\tau}`
t is the spatial domain at time :math:`{\tau}` and :math:`\Omega_{\tau+\Delta t}` is the spatial domain occupied by the same
material points at time :math:`{\tau+\Delta t}`. The notation on the left-hand side is a little confusing because it
appears to refer to a single spatial domain. However, in this notation, which is standard, the material derivative on the integral implies that the domain Ω is a material domain. We now
transform both integrals on the RHS to the reference domain:

.. math::
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }f(\mathbf{x},t)\mathrm{d}\Omega 
  =\lim_{\Delta  t \to 0}\cfrac{1}{\Delta  t }\left ( \int\limits_{\Omega_{0} }f(\mathbf{X},{\tau+\Delta t})J(\mathbf{X},{\tau+\Delta t})\mathrm{d}\Omega_{0}-\int\limits_{\Omega_{0} }f(\mathbf{X},\tau )J(\mathbf{X},\tau )\mathrm{d}\Omega_{0} \right )
  
With this change in the domain of integration, :math:`{f}` becomes a function of the material coordinates, that is, :math:`f(\mathbf{\Phi}(\mathbf{X},t),t)\equiv f\circ \mathbf{\Phi}`.

Since the domain of integration is now independent of time, we can pull the limit operation
inside the integral and take the limit, which yields  

.. math::
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }f(\mathbf{x},t)\mathrm{d}\Omega 
  =\int\limits_{\Omega_{0} }\frac{\partial }{\partial t}[f(\mathbf{X},t)J(\mathbf{X},t)] \mathrm{d}\Omega_{0}
  
The partial derivative with respect to time in the integrand is a material time derivative since
the independent space variables are the material coordinates. We next use the product rule for
derivatives on the previous:

.. math::
  \begin{align}
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }f(\mathbf{x},t)\mathrm{d}\Omega 
  =\int\limits_{\Omega_{0} }\frac{\partial }{\partial t}[f(\mathbf{X},t)J(\mathbf{X},t)] \mathrm{d}\Omega_{0}\\
  =\int\limits_{\Omega_{0} }(\frac{\partial f(\mathbf{X},t)}{\partial t}J(\mathbf{X},t)+f(\mathbf{X},t)\frac{\partial J(\mathbf{X},t)}{\partial t}) \mathrm{d}\Omega_{0}\\
  \end{align}

-

.. math::  
  \begin{align}
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }f(\mathbf{x},t)\mathrm{d}\Omega 
  =\int\limits_{\Omega_{0} }\frac{\partial }{\partial t}[f(\mathbf{X},t)J(\mathbf{X},t)] \mathrm{d}\Omega_{0}\\
  =\int\limits_{\Omega_{0} }(\frac{\partial f(\mathbf{X},t)}{\partial t}J(\mathbf{X},t)+f(\mathbf{X},t)\frac{\partial J(\mathbf{X},t)}{\partial t}) \mathrm{d}\Omega_{0}\\
  =\int\limits_{\Omega_{0} }(\frac{\partial f(\mathbf{X},t)}{\partial t}J(\mathbf{X},t)+f(\mathbf{X},t)J(\mathbf{X},t)\frac{\partial v_{k}}{\partial x_{k}}) \mathrm{d}\Omega_{0}\\
  \end{align}
  
We can now transform the RHS integral to the current domain and change the
independent variables to an Eulerian description, which gives

.. math:: 
  \begin{align}
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }f(\mathbf{x},t)\mathrm{d}\Omega 
  =\int\limits_{\Omega }(\frac{D f(\mathbf{x},t)}{D t}+f(\mathbf{x},t)\frac{\partial v_{k}}{\partial x_{k}}) \mathrm{d}\Omega\\
  \end{align}
  
where we have used :math:`{D f(\mathbf{x},t)}/{D t}\equiv \partial f(\mathbf{X},t)/\partial t`. This is one form of
Reynolds’ transport theorem. 
 
.. math:: 
  \cfrac{\mathrm{D}  f}{\mathrm{D} t}=\cfrac{\partial f}{\partial t}+v_{i} \cfrac{\partial f}{\partial x_{i}}=\cfrac{\partial f}{\partial t}+\mathbf{v} \cdot \nabla f=\cfrac{\partial f}{\partial t}+\mathbf{v} \cdot \operatorname{grad} f \\
  
-  
  
.. math::   
  \begin{array}{l}
  \cfrac{\mathrm{D} f(\mathbf{x},t)}{\mathrm{D} t}&=\cfrac{\partial f(\mathbf{x},t)}{\partial t}+v_{i}(\mathbf{x},t) \cfrac{\partial f(\mathbf{x},t)}{\partial x_{i}}\\
  &=\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\mathbf{v}(\mathbf{x},t) \cdot \nabla f(\mathbf{x},t)\\
  &=\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\mathbf{v}(\mathbf{x},t) \cdot \operatorname{grad} f(\mathbf{x},t) \\
  \end{array}  

-  
  
.. math:: 
  \begin{align}
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }f(\mathbf{x},t)\mathrm{d}\Omega 
  &=\int\limits_{\Omega }(\frac{D f(\mathbf{x},t)}{D t}+f(\mathbf{x},t)\frac{\partial v_{k}}{\partial x_{k}}) \mathrm{d}\Omega\\
  &=\int\limits_{\Omega }(\cfrac{\partial f(\mathbf{x},t)}{\partial t}+v_{k}(\mathbf{x},t) \cfrac{\partial f(\mathbf{x},t)}{\partial x_{k}}+f(\mathbf{x},t)\frac{\partial v_{k}}{\partial x_{k}}) \mathrm{d}\Omega\\
  &=\int\limits_{\Omega }(\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\cfrac{\partial (v_{k}(\mathbf{x},t)f(\mathbf{x},t))}{\partial x_{k}}) \mathrm{d}\Omega\\
  \end{align}
  
which can be written in tensor form as

.. math:: 
  \begin{align}
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }f(\mathbf{x},t)\mathrm{d}\Omega 
  &=\int\limits_{\Omega }(\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\mathrm{div}\ (\mathbf{v}(\mathbf{x},t)\cdot f(\mathbf{x},t)) \mathrm{d}\Omega\\
  &=\int\limits_{\Omega }(\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\nabla\cdot(\mathbf{v}(\mathbf{x},t)\cdot f(\mathbf{x},t)) \mathrm{d}\Omega\\
  \end{align}

-
  
.. math::
  \begin{align}
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }f(\mathbf{x},t)\mathrm{d}\Omega 
  &=\int\limits_{\Omega }(\cfrac{\partial f(\mathbf{x},t)}{\partial t}) \mathrm{d}\Omega
  +\int\limits_{\Omega }(\nabla\cdot(\mathbf{v}(\mathbf{x},t)\cdot f(\mathbf{x},t)) \mathrm{d}\Omega\\
  &=\int\limits_{\Omega }(\cfrac{\partial f(\mathbf{x},t)}{\partial t}) \mathrm{d}\Omega
  +\int\limits_{\Gamma }f(\mathbf{x},t)(\mathbf{v}(\mathbf{x},t)\cdot \mathbf{n}(\mathbf{x},t)) \mathrm{d}\Gamma \\
  &=\int\limits_{\Omega }(\cfrac{\partial f(\mathbf{x},t)}{\partial t}) \mathrm{d}\Omega
  +\int\limits_{\Gamma }f(\mathbf{x},t)({v}_{k}(\mathbf{x},t){n}_{k}(\mathbf{x},t)) \mathrm{d}\Gamma \\
  \end{align}
  
Mass Conservation
-------------------

The mass :math:`m(\Omega)` of a material domain:math:`\Omega` is given by

.. math::
  m(\Omega)=\int\limits_{\Omega }\rho (\mathbf{x},t) \mathrm{d}\Omega  
  
where :math:`\rho (\mathbf{x},t)` is the density. Mass conservation requires that the mass of any material domain
be constant, since no material flows through the boundaries of a material domain and we are
not considering mass to energy conversion. Therefore, according to the principle of mass
conservation, the material time derivative of :math:`m(\Omega)` vanishes, that is,  

.. math::
  \cfrac{\mathrm{D} m(\Omega)}{\mathrm{D} t} =\cfrac{\mathrm{D}}{\mathrm{D} t}\int\limits_{\Omega }\rho (\mathbf{x},t) \mathrm{d}\Omega=0
  
Applying Reynolds’ theorem

.. math::
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }\rho(\mathbf{x},t)\mathrm{d}\Omega 
  =\int\limits_{\Omega }(\frac{D \rho(\mathbf{x},t)}{D t}+\rho(\mathbf{x},t)\frac{\partial v_{k}}{\partial x_{k}}) \mathrm{d}\Omega=0\\

-

.. math::
  \frac{\mathrm{D} }{\mathrm{D} t}\int\limits_{\Omega }\rho(\mathbf{x},t)\mathrm{d}\Omega 
  =\int\limits_{\Omega }(\frac{D \rho(\mathbf{x},t)}{D t}+\rho(\mathbf{x},t)\nabla \mathbf{v}) \mathrm{d}\Omega=0\\  
  
-

.. math::  
  \frac{D \rho(\mathbf{x},t)}{D t}+\rho(\mathbf{x},t)\frac{\partial v_{k}}{\partial x_{k}}=0
  
- 
 
.. math::  
  \frac{D \rho(\mathbf{x},t)}{D t}+\rho(\mathbf{x},t)\nabla \cdot\mathbf{v}=0  
