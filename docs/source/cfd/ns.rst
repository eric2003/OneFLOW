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
