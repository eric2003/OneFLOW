HLLC Riemann Solver
==================================

#. `CFD Julia: A Learning Module Structuring an Introductory Course on Computational Fluid Dynamics <https://www.mdpi.com/2311-5521/4/3/159/>`_
#. `Eigenstructure and Approximate Riemann Solvers for Hyperbolic Conservation Laws <https://www3.nd.edu/~dbalsara/Numerical-PDE-Course/Appendix_LesHouches/LesHouches_Lecture_5_Approx_RS.pdf>`_
#. `Damien Furfaro: A simple HLLC-type Riemann solver for compressible non-equilibrium two-phase flows <https://www.youtube.com/watch?v=x7g82GaYSXM/>`_
#. `Numerical Fluid Dynamics <https://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2011/>`_



Reference
----------------------------
- Riemann Solvers and Numerical Methods for Fluid Dynamics

HLLC Riemann Solver
----------------------------
HLLC Riemann Solver used lower and upper bounds on the characteristics speeds in the solution of Riemann problem involving left and right states. These bounds are approximated as

.. math::
  \begin{array}{c}
  S_{L}=\text{min}(u_{L},u_{R})-\text{max}(a_{L},a_{R})\\
  S_{R}=\text{max}(u_{L},u_{R})+\text{max}(a_{L},a_{R})\\
  \end{array}
  
where :math:`S_{L}` and :math:`S_{R}` are the lower and upper bound for the left and right state characteristics speed. For HLLC
scheme we also include middle wave of speed :math:`S_{*}` given by

.. math::
  S_{*}=\cfrac{p_{R}-p_{L}+\rho_{L}u_{L}(S_{L}-u_{L})-\rho_{R}u_{R}(S_{R}-u_{R})}
  {\rho_{L}(S_{L}-u_{L})-\rho_{R}(S_{R}-u_{R})}
  
The mean pressure is given by

.. math::
  P_{LR}=\cfrac{1}{2}(p_{L}+p_{R}+\rho_{L}(S_{L}-u_{L})(S_{*}-u_{L})+\rho_{R}(S_{R}-u_{R})(S_{*}-u_{R}))  
  
The fluxes are computed as

.. math::
  F_{i+1/2}=\left\{
  \begin{array}{ll}
  F^{L},& \text{if } S_{L}\ge 0\\
  F_{*L},& \text{if }S_{L}\le 0\le S_{*}\\
  F_{*R},& \text{if }S_{*}\le 0\le S_{R}\\
  F^{R},& \text{if }S_{R}\le 0\\
  \end{array}
  \right.
  
-
  
.. math::
  \begin{array}{l}
  \mathbf{U}_{*L}=\cfrac{S_{L}\mathbf{U}_{L}-\mathbf{F}_{L}+P_{LR}\mathbf{D}_{*}}{S_{L}-S_{*}}\\
  \mathbf{U}_{*R}=\cfrac{S_{R}\mathbf{U}_{R}-\mathbf{F}_{R}+P_{LR}\mathbf{D}_{*}}{S_{R}-S_{*}}\\
  \mathbf{F}_{*L}=\cfrac{S_{*}(S_{L}\mathbf{U}_{L}-\mathbf{F}_{L})+S_{L}P_{LR}\mathbf{D}_{*}}{S_{L}-S_{*}}\\
  \mathbf{F}_{*R}=\cfrac{S_{*}(S_{R}\mathbf{U}_{R}-\mathbf{F}_{R})+S_{R}P_{LR}\mathbf{D}_{*}}{S_{R}-S_{*}}\\
  \end{array}  

-
  
.. math::
  \mathbf{U}_{*K}=\cfrac{S_{K}\mathbf{U}_{K}-\mathbf{F}_{K}+P_{LR}\mathbf{D}_{*}}{S_{K}-S_{*}}  
  
-
  
.. math::
  \mathbf{F}_{*K}=\cfrac{S_{*}(S_{K}\mathbf{U}_{K}-\mathbf{F}_{K})+S_{K}P_{LR}\mathbf{D}_{*}}{S_{K}-S_{*}}\\  

-

.. math::
  \mathbf{U}=\begin{bmatrix}
  \rho\\\rho u\\\rho E
  \end{bmatrix}\quad 
  \mathbf{F}=\begin{bmatrix}
  \rho u\\\rho u^2+p\\\rho H u
  \end{bmatrix}\quad 
  \mathbf{D}_{*}=\begin{bmatrix}
  0\\1\\S_{*}
  \end{bmatrix}    

-

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
  