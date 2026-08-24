One-Dimensional Euler Solver
==================================

#. `CFD Julia: A Learning Module Structuring an Introductory Course on Computational Fluid Dynamics <https://www.mdpi.com/2311-5521/4/3/159/>`_


The one-dimensional Euler equations in its conservative form can be written as

.. math::
  \cfrac{\partial \mathbf{q}}{\partial t}+\cfrac{\partial \mathbf{F}}{\partial x}=0
  
where

.. math::
  \mathbf{q}=\begin{bmatrix}
   \rho\\ \rho u \\\rho E
  \end{bmatrix}\quad
  \mathbf{F}=\begin{bmatrix}
   \rho u\\ \rho uu + p \\\rho H u
  \end{bmatrix} 
  
in which

.. math::
  \begin{array}{l}
  H=E+\cfrac{p}{\rho}\\
  E=\cfrac{1}{\gamma-1}\cfrac{p}{\rho}+\cfrac{1}{2}\mathbf{v}^2
  \end{array}
  
Here, :math:`\rho` and :math:`p` are the density, and pressure respectively; :math:`u` is the horizontal component of velocity;
:math:`E` and :math:`H` stand for the specific internal energy and specific enthalpy('specific' means per-unit-mass), and :math:`\gamma` is the ratio of specific heats. The above equationscan also
be written as

.. math::
  \cfrac{\partial \mathbf{q}}{\partial t}+\mathbf{A}\cfrac{\partial \mathbf{q}}{\partial x}=0

where :math:`\mathbf{A}=\cfrac{\partial \mathbf{F}}{\partial \mathbf{q}}` is the convective flux Jacobian matrix. The matrix 
:math:`\mathbf{A}` is given as

.. math::
  \mathbf{A}=\cfrac{\partial \mathbf{F}}{\partial \mathbf{q}}=\begin{bmatrix}
  0& 1 & 0\\
  (\gamma-3)\cfrac{u^2}{2}& (3-\gamma)u & (\gamma-1)\\
  (\cfrac{\gamma-1}{2}u^{2}-H)u&H+(1-\gamma)u^{2}  &\gamma u
  \end{bmatrix}
  
or

.. math::
  \mathbf{A}=\cfrac{\partial \mathbf{F}}{\partial \mathbf{q}}=\begin{bmatrix}
  0& 1 & 0\\
  \phi ^{2}-u^2& (3-\gamma)u & (\gamma-1)\\
  (\phi ^{2}-H)u&H+(1-\gamma)u^{2}  &\gamma u
  \end{bmatrix}  

where :math:`\phi^{2}=\cfrac{\gamma-1}{2}q^2`.

in which

.. math::
  \begin{array}{l}
  q^2=u^{2}\quad\quad\quad\quad\quad\quad(1d)\\
  q^2=u^{2}+v^{2}\quad\quad\quad\quad(2d)\\
  q^2=u^{2}+v^{2}+w^{2}\ \ \ \quad(3d)\\
  \end{array}  
  
Eigenstructure

.. math::
 \mathbf{A}=\mathbf{R}\Lambda \mathbf{L}
 
where 

.. math::
  \Lambda=\begin{bmatrix}
  u &  0 & 0\\
  0& u+c & 0\\
  0& 0 &u-c
  \end{bmatrix}
  
-
  
.. math::
  \mathbf{R}=\begin{bmatrix}
  1 &  \beta & \beta\\
  u& \beta (u+c) & \beta (u-c)\\
  \cfrac{\phi^{2}}{\gamma-1}& \beta (H+uc) &\beta (H-uc)
  \end{bmatrix}  
  
-
  
.. math::
  \mathbf{L}=\begin{bmatrix}
  1- \cfrac{\phi^{2}}{c^2}&  (\gamma-1)\cfrac{u}{c^{2}} & -\cfrac{(\gamma-1)}{c^{2}}\\
  \phi^{2}-uc& -(\gamma-1)u+c & \gamma-1\\
  \phi^{2}+uc& -(\gamma-1)u-c &\gamma-1
  \end{bmatrix}
  
where :math:`c` is the speed of sound and is given by :math:`c^{2}=\gamma p/\rho`, and :math:`\beta=1/(2c^{2})`.
The semi-discrete form of the Euler equations can be written as

.. math::
  \cfrac{\partial \mathbf{q}_{i}}{\partial t}+\cfrac{\mathbf{F}_{i+1/2}-\mathbf{F}_{i-1/2}}{\Delta x}=0
  
where :math:`\mathbf{q}_{i}` is cell-centered values stored at nodal points and :math:`\mathbf{F}_{i-1/2}` and :math:`\mathbf{F}_{i+1/2}` are
the fluxes at left and right cell interface.
We use third-order Runge-Kutta numerical method for time integration. We use WENO-5 reconstruction to compute the left and sight states of the fluxes at the interface. 

Roe’s Riemann Solver
---------------------------

.. math::
  \mathbf{F}_{i+1/2}=\cfrac{1}{2}(\mathbf{F}^{L}_{i+1/2}+\mathbf{F}^{R}_{i+1/2})
  -\cfrac{1}{2}\bar{\mathbf{R}}|\bar{\Lambda}|\bar{\mathbf{L}}(\mathbf{q}^{R}_{i+1/2}-\mathbf{q}^{L}_{i+1/2})
  
First we reconstruct the left and right states of :math:`\mathbf{q}` at the interface (i.e., :math:`\mathbf{q}^{L}_{i+1/2}` and :math:`\mathbf{q}^{R}_{i+1/2}` ) using WENO-5 scheme.
Then, we can compute the left and right state of the fluxes (i.e., :math:`\mathbf{F}^{L}`  and :math:`\mathbf{F}^{R}`).
The Roe averaging formulas to compute approximate values for constructing :math:`\bar{\mathbf{R}}` and :math:`\bar{\mathbf{L}}`
are given below

.. math::
  \begin{array}{c}
  \bar{u}=\cfrac{u^{L}\sqrt{\rho^{L}}+u^{R}\sqrt{\rho^{R}}}{\sqrt{\rho^{L}}+\sqrt{\rho^{R}}}\\
  \bar{H}=\cfrac{H^{L}\sqrt{\rho^{L}}+H^{R}\sqrt{\rho^{R}}}{\sqrt{\rho^{L}}+\sqrt{\rho^{R}}}
  \end{array}
  
where the left and right states are computed using WENO-5 reconstruction. The speed of the sound is computed using the below equation  

.. math::
  \bar{c}=\sqrt{(\gamma-1)\left[\bar{H}-\cfrac{1}{2}\bar{u}^{2}\right]}
  
The eigenvalues of the Jacobian matrix are :math:`\lambda_{1}=\bar{u}`, :math:`\lambda_{2}=\bar{u}+\bar{c}`, and :math:`\lambda_{3}=\bar{u}-\bar{c}`.

  
.. math::
  \begin{array}{l}
  \mathbf{q}[0]=\rho\\
  \mathbf{q}[1]=\rho u\\
  \mathbf{q}[2]=\rho E\\
  \end{array}
  
-
  
.. math::
  \begin{array}{l}
  \rho E=\cfrac{p}{\gamma-1}+\cfrac{1}{2}\rho (u^2+v^2+w^2)\\
  (\gamma-1)\rho E=p+\cfrac{(\gamma-1)}{2}\rho (u^2+v^2+w^2)\\
  p=(\gamma-1)\rho E-\cfrac{(\gamma-1)}{2}\rho (u^2+v^2+w^2)\\
  p=(\gamma-1)\rho E-\cfrac{(\gamma-1)}{2} \cfrac{((\rho u)^2+(\rho v)^2+(\rho w)^2)}{\rho}\\
  p=(\gamma-1)\mathbf{q}[2]-\cfrac{(\gamma-1)}{2} \cfrac{(\mathbf{q}[1]^2)}{\mathbf{q}[0]}\\
  p=(\gamma-1)(\mathbf{q}[2]-\cfrac{1}{2} \cfrac{\mathbf{q}[1]^2}{\mathbf{q}[0]})\\
  \end{array}  

-
  
.. math::
  u=\cfrac{\rho u}{\rho}=\mathbf{q}[1]/\mathbf{q}[0]\\  
  
-
  
.. math::
  \begin{array}{l}
  H=E+\cfrac{p}{\rho}\\
  {\rho}H={\rho}E+p\\
  u=\cfrac{\rho u}{\rho}=\mathbf{q}[1]/\mathbf{q}[0]\\
  {\rho}Hu={\rho}Eu+pu=\mathbf{q}[2]\mathbf{q}[1]/\mathbf{q}[0]+p\mathbf{q}[1]/\mathbf{q}[0]\\
  \end{array}  
  
-
  
.. math::
  \mathbf{q}=\begin{bmatrix}
   \rho\\ \rho u \\\rho E
  \end{bmatrix}=\begin{bmatrix}
   \mathbf{q}[0]\\ \mathbf{q}[1] \\\mathbf{q}[2]
  \end{bmatrix}  
  
-
  
.. math::
  \mathbf{F}=\begin{bmatrix}
  \rho u\\ \rho uu + p \\\rho H u
  \end{bmatrix}= \begin{bmatrix}
  \mathbf{q}[1]\\ \cfrac{(\mathbf{q}[1])^{2}}{\mathbf{q}[0]}+p \\
  \mathbf{q}[2]\mathbf{q}[1]/\mathbf{q}[0]+p\mathbf{q}[1]/\mathbf{q}[0]
  \end{bmatrix}
  
-
  
.. math::  
  p=(\gamma-1)(\mathbf{q}[2]-\cfrac{1}{2} \cfrac{\mathbf{q}[1]^2}{\mathbf{q}[0]})\\  