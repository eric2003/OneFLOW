Rusanov Riemann Solver
==================================

#. `CFD Julia: A Learning Module Structuring an Introductory Course on Computational Fluid Dynamics <https://www.mdpi.com/2311-5521/4/3/159/>`_



Rusanov Riemann Solver
----------------------------
The Rusanov scheme uses maximum local wave propagation speed to compute the flux as follows

.. math::
  f_{i+1/2}=\frac{1}{2}(f^{L}_{i+1/2}+f^{R}_{i+1/2})-\frac{1}{2}c_{i+1/2}(u^{R}_{i+1/2}-u^{L}_{i+1/2})
  
where :math:`f^{L}` is the flux component using the right reconstructed state :math:`f^{L}_{i+1/2}=f(u^{L}_{i+1/2})` and :math:`f^{R}`
is the flux component using the right reconstructed state :math:`f^{R}_{i+1/2}=f(u^{R}_{i+1/2})`. Here, :math:`c_{i+1/2}` 
is the local wave propagation speed which is obtained by taking the maximum absolute value of the eigenvalues corresponding to
the Jacobian, :math:`\cfrac{\partial f}{\partial u}=u`, between cells :math:`i` and :math:`i+1` can be obtained as

.. math::
  c_{i+1/2}=\text{max}(|u_{i}|,|u_{i+1}|)
  
or in a wider stencil as shown in Equation   

.. math::
  c_{i+1/2}=\left | \cfrac{\partial f}{\partial u} \right |=max(|u_{i-2}|,|u_{i-1}|,|u_{i}|,|u_{i+1}|,|u_{i+2}|) 
  
The Riemann solver based on Rusanov scheme is simple compared to Roe’s Riemann solver and HLLC based Riemann solver. For Euler equations, instead of solving one equation (as in inviscid Burgers equation), now we have to follow the procedure for three equations
(i.e., density, velocity, and energy). We need to approximate wave propagation speed at the interface :math:`c_{i+1/2}` to compute flux at the interface. For Rusanov scheme, we simply use maximum eigenvalue of the Jacobian matrix as the wave propagation speed. We have 
:math:`c_{i+1/2}=\text{max}(|\bar{u}|,|\bar{u}+\bar{a}|,|\bar{u}-\bar{a}|)`, where :math:`\bar{u}` and :math:`\bar{a}` are are computed using Roe averaging.

The one-dimensional Euler equations in its conservative form can be written as

.. math::
  \cfrac{\partial \mathbf{U}}{\partial t}+\cfrac{\partial \mathbf{F}}{\partial x}=0
  
where

.. math::
  \mathbf{U}=\left[\begin{array}{c}
   \rho\\ \rho u \\\rho E
  \end{array}\right]\quad
  \mathbf{F}=\left[\begin{array}{c}
   \rho u\\ \rho uu + p \\\rho H u
  \end{array} \right]
  
Rusanov Flux
----------------------------  

.. math::
  \mathbf{F}_{i+1/2}=\cfrac{1}{2}(\mathbf{F}_{L}+\mathbf{F}_{R})-\cfrac{1}{2}S^{+}(\mathbf{U}_{R}-\mathbf{U}_{L})
  
Actually, the above speed is bounded by

.. math::
  S^{+}=\text{max}\left\{ |u_{L}|+a_{L},|u_{R}|+a_{R}\right\}  

where
  
.. math::
  \mathbf{U}_{L}=\left[\begin{array}{l}
  \rho_{L}\\\rho_{L} u_{L}\\\rho_{L} E_{L}
  \end{array}\right]\quad 
  \mathbf{U}_{R}=\left[\begin{array}{l}
  \rho_{R}\\\rho_{R} u_{R}\\\rho_{R} E_{R}
  \end{array}\right]
  
-

.. math::
  \mathbf{F}_{L}=\left[\begin{array}{l}
  \rho_{L}u_{L}\\\rho_{L} u^{2}_{L}+p_{L}\\\rho_{L} H_{L}u_{L}
  \end{array}\right]\quad
  \mathbf{F}_{R}=\left[\begin{array}{l}
  \rho_{R}u_{R}\\\rho_{R} u^{2}_{R}+p_{R}\\\rho_{R} H_{R}u_{R}
  \end{array}\right]  