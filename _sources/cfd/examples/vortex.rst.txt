Vortex-Merger Problem
==================================

#. `CFD Julia: A Learning Module Structuring an Introductory Course on Computational Fluid Dynamics <https://www.mdpi.com/2311-5521/4/3/159/>`_
#. `San, O.; Staples, A.E. High-order methods for decaying two-dimensional homogeneous isotropic turbulence. Comput. Fluids 2012, 63, 105–127.`
#. `Meunier, P.; Le Dizes, S.; Leweke, T. Physics of vortex merging. Comptes Rendus Phys. 2005, 6, 431–450.`
#. `San, O.; Staples, A.E. A coarse-grid projection method for accelerating incompressible flow computations. J. Comput. Phys. 2013, 233, 480–508.`
#. `Code Verification for the SENSEI CFD Code <https://asmedigitalcollection.asme.org/verification/article-abstract/8/2/021003/1163559/Code-Verification-for-the-SENSEI-CFD-Code?redirectedFrom=fulltext>`_
#. `Weicheng Xue <https://scholar.google.com/citations?user=tRtBtbIAAAAJ&hl=en>`_


Vortex-Merger
-----------------------------------------------------
We demonstrate the two-dimensional Navier-Stokes solver for the domain with periodic boundary condition using vortex-merger problem. 
The merging process occurs when two vortices of the same sign with parallel axes are within a certain critical distance from each other. 
The two vortices merge to form a single vortex. It is a two-dimensional process and is one of the fundamental processes of fluid motion and it plays a key role in a variety of simulations, such as decaying two-dimensional turbulence, three-dimensional turbulence, and mixing layers.
This phenomenon also occurs in other fields such as astrophysics, meteorology, and geophysics.

We use the Cartesian computational domain :math:`(x,y)\in[0,2\pi]\times[0,2\pi]` and divide it into :math:`128\times 128` grid resolution.
The vorticity equation for the vortex-merger problem is same as the lid-driven cavity problem. We use :math:`Re=2000` and integrate the solution from time
:math:`t=0` to :math:`t=20` with :math:`\Delta t=0.01`. We use third-order Runge-Kutta method for time integration similar to the lid-driven cavity problem.

Taylor–Green vortex
-----------------------------------------------------
The Taylor–Green vortex problem is the two-dimensional,
unsteady flow of a decaying vortex, which is an exact closed form
solution to the incompressible Navier–Stokes equations in Cartesian coordinates. The exact solution of this vortex flow in a
:math:`[0,2\pi]\times[0,2\pi]` domain with periodic boundary conditions is given by

.. math::
  \omega^{e}(x,y,t)=2\kappa \cos(\kappa x)\cos(\kappa y)e^{-2\kappa^{2}t/Re}
  
where :math:`\kappa` is an integer which represents the number of vortices in
each direction. In order to quantify the effective order of accuracy
for each scheme, using the difference between exact and computed
solutions, we compute the discrete :math:`L_{2}` norm as  

.. math::
  \|\omega\|_{L_{2}}=\sqrt{\cfrac{1}{N_xN_y}\sum_{i=1}^{N_{x}}\sum_{j=1}^{N_{y}}|\omega_{i,j}^{e}-\omega_{i,j}|^{2}}
  
In order to test the spatial convergence rate, we first perform
numerical simulations for :math:`Re = 1` and :math:`\kappa = 4` with all the spatial
schemes introduced earlier using the third-order TVD Runge–Kutta
scheme with :math:`\Delta t=10^{-4}`.  We choose this small time step (wherein
the maximum CFL number is around 0.002 for the :math:`128^{2}` resolution case)
to eliminate possible temporal discretization errors due to the Runge–Kutta time integration schemes.
A study with :math:`\Delta t=2\times 10^{-4},10^{-4}`, and :math:`\Delta t=5\times 10^{-5}` shows that the predicted flow
field is independent of :math:`\Delta t`. To further verify whether the simulations are independent of the errors associated with temporal
discretization, we also run the same computations using the fourth-order Runge Kutta scheme, which yields exactly the same results.