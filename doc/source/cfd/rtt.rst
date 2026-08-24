Reynolds' transport theorem
==================================

- Fluid Mechanics Fundamentals and Applications, Fourth Edition, Yunus Cengel, John Cimbala

pp-167

RTT, fixed CV:

.. math::
  \cfrac{\text{d}B_{sys}}{\text{d}t}=\cfrac{\text{d}}{\text{d}t}\int_{\text{CV}}\rho b\text{d}V
  +\int_{\text{CS}}\rho b\mathbf{v}\cdot\mathbf{n}\text{d}A
  
Alternate RTT, fixed CV:  

.. math::
  \cfrac{\text{d}B_{sys}}{\text{d}t}=\int_{\text{CV}}\cfrac{\partial}{\partial t} (\rho b)\text{d}V
  +\int_{\text{CS}}\rho b\mathbf{v}\cdot\mathbf{n}\text{d}A
  
RTT, nonfixed CV:

.. math::
  \cfrac{\text{d}B_{sys}}{\text{d}t}=\cfrac{\text{d}}{\text{d}t}\int_{\text{CV}}\rho b\text{d}V
  +\int_{\text{CS}}\rho b(\mathbf{v}-\mathbf{v}_{\text{CS}})\cdot\mathbf{n}\text{d}A
  
Alternate RTT, nonfixed CV:

.. math::
  \cfrac{\text{d}B_{sys}}{\text{d}t}=\int_{\text{CV}}\cfrac{\partial}{\partial t} (\rho b)\text{d}V
  +\int_{\text{CS}}\rho b\mathbf{v}\cdot\mathbf{n}\text{d}A
    
  