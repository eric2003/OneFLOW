Conservation Of Mass
==================================

- Fluid Mechanics Fundamentals and Applications, Fourth Edition, Yunus Cengel, John Cimbala

pp-445

Derivation Using the Divergence Theorem
``````````````````````````````````````````

Conservation of mass for a CV:

.. math::
  \cfrac{\text{d}m_{sys}}{\text{d}t}=\cfrac{\text{d}}{\text{d}t}\int_{\text{CV}}\rho \text{d}V
  +\int_{\text{CS}}\rho (\mathbf{v}-\mathbf{v}_{\text{CS}})\cdot\mathbf{n}\text{d}A
 
-  

.. math::
  \int_{\text{CV}}\cfrac{\partial \rho}{\partial t}dV+\int_{\text{CS}}\rho \vec{V}\cdot \vec{n}dA=0

- 
  
.. math::
  \cfrac{\text{d}m_{sys}}{\text{d}t}=\cfrac{\text{d}}{\text{d}t}\int_{\text{CV}}\rho \text{d}V
  +\int_{\text{CS}}\rho (\mathbf{v}-\mathbf{v}_{\text{CS}})\cdot\mathbf{n}\text{d}A=
  \int_{\text{CV}}\cfrac{\partial \rho}{\partial t}dV+\int_{\text{CS}}\rho \vec{V}\cdot \vec{n}dA=0
  
  
Divergence theorem:

.. math::
  \int_{V}\nabla\cdot\vec{G}=\oint_{A}\vec{G}\cdot\vec{n}dA

-

.. math::
  \int_{\text{CV}}\cfrac{\partial \rho}{\partial t}dV+\int_{\text{CV}}\nabla\cdot(\rho \vec{V})dV=0 

We now combine the two volume integrals into one,

.. math::
  \int_{\text{CV}}\left[\cfrac{\partial \rho}{\partial t}+\nabla\cdot(\rho \vec{V})\right]dV=0

Continuity equation:

.. math::
  \cfrac{\partial \rho}{\partial t}+\nabla\cdot(\rho \vec{V})=0
  
Alternative Form of the Continuity Equation  

.. math::
  \cfrac{\partial \rho}{\partial t}+\nabla\cdot(\rho \vec{V})=\underbrace{\cfrac{\partial \rho}{\partial t}
  +\vec{V}\cdot\nabla\rho}_{\text{Material derivative of } \rho}+\rho\nabla\cdot\vec{V}=0

-
  
.. math::
  \cfrac{\text{d}\rho}{\text{d} t}+\rho\nabla\cdot\vec{V}=0
