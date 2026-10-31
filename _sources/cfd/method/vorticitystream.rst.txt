Vorticity Stream Function
==================================

#. `Vorticity Stream Function Formulation <http://www.fem.unicamp.br/~phoenics/SITE_PHOENICS/Apostilas/CFD-1_U%20Michigan_Hong/Lecture05.pdf>`_

Incompressible N-S Equation in 2D
--------------------------------------

.. math::
  \begin{array}{l}
  \cfrac{\partial u}{\partial t} +u\cfrac{\partial u}{\partial x}+v\cfrac{\partial u}{\partial y}=
  -\cfrac{1}{\rho}\cfrac{\partial p}{\partial x}+\nu\bigg(\cfrac{\partial ^{2}u}{\partial x^{2}}+\cfrac{\partial ^{2}u}{\partial y^{2}}\bigg)\\
  \cfrac{\partial v}{\partial t} +u\cfrac{\partial v}{\partial x}+v\cfrac{\partial v}{\partial y}=
  -\cfrac{1}{\rho}\cfrac{\partial p}{\partial y}+\nu\bigg(\cfrac{\partial ^{2}v}{\partial x^{2}}+\cfrac{\partial ^{2}v}{\partial y^{2}}\bigg)\\
  \end{array}
  
Nondimensionalization
---------------------------  

.. math::
  \begin{array}{l}
  \widetilde{u}=\cfrac{u}{U},\widetilde{x}=\cfrac{x}{L},t=t\cfrac{U}{L},
  \widetilde{p}=\cfrac{p}{\rho U^{2}}\\
  \end{array}  
  
-
  
.. math::
  \begin{array}{l}
  \cfrac{\partial \widetilde{u}}{\partial \widetilde{t}} +\widetilde{u}\cfrac{\partial \widetilde{u}}{\partial \widetilde{x}}+\widetilde{v}\cfrac{\partial \widetilde{u}}{\partial \widetilde{y}}=
  -\cfrac{1}{\rho}\cfrac{\partial p}{\partial \widetilde{x}}+\cfrac{1}{Re}\bigg(\cfrac{\partial ^{2}\widetilde{u}}{\partial \widetilde{x}^{2}}+\cfrac{\partial ^{2}\widetilde{u}}{\partial \widetilde{y}^{2}}\bigg)\\
  \cfrac{\partial \widetilde{v}}{\partial \widetilde{t}} +\widetilde{u}\cfrac{\partial \widetilde{v}}{\partial \widetilde{x}}+\widetilde{v}\cfrac{\partial \widetilde{v}}{\partial \widetilde{y}}=
  -\cfrac{1}{\rho}\cfrac{\partial p}{\partial \widetilde{y}}+\cfrac{1}{Re}\bigg(\cfrac{\partial ^{2}\widetilde{v}}{\partial \widetilde{x}^{2}}+\cfrac{\partial ^{2}\widetilde{v}}{\partial \widetilde{y}^{2}}\bigg)\\
  \end{array}  
  
The Vorticity Equation
---------------------------

.. math::
  \begin{align}
  -\cfrac{\partial}{\partial y} &\bigg[\cfrac{\partial {u}}{\partial {t}} +{u}\cfrac{\partial {u}}{\partial {x}}+{v}\cfrac{\partial {u}}{\partial {y}}  = -\cfrac{1}{\rho}\cfrac{\partial p}{\partial {x}}+\cfrac{1}{Re}\bigg(\cfrac{\partial ^{2}{u}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{u}}{\partial {y}^{2}}\bigg)\bigg]\\
  \cfrac{\partial}{\partial x} &\bigg[\cfrac{\partial {v}}{\partial {t}} +{u}\cfrac{\partial {v}}{\partial {x}}+{v}\cfrac{\partial {v}}{\partial {y}}  = -\cfrac{1}{\rho}\cfrac{\partial p}{\partial {y}}+\cfrac{1}{Re}\bigg(\cfrac{\partial ^{2}{v}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{v}}{\partial {y}^{2}}\bigg)\bigg]\\
  \end{align}

-
  
.. math::
  \begin{array}{l}
  \omega=\cfrac{\partial v}{\partial x} -\cfrac{\partial u}{\partial y}\\
  -\cfrac{\partial}{\partial y} \bigg(\cfrac{\partial ^{2}{u}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{u}}{\partial {y}^{2}}\bigg)
  +\cfrac{\partial}{\partial x}\bigg(\cfrac{\partial ^{2}{v}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{v}}{\partial {y}^{2}}\bigg)
  =\bigg(\cfrac{\partial ^{2}{\omega}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{\omega}}{\partial {y}^{2}}\bigg)
  \end{array}  

-
  
.. math::
  \cfrac{\partial {\omega}}{\partial {t}} +{u}\cfrac{\partial {\omega}}{\partial {x}}+{v}\cfrac{\partial {\omega}}{\partial {y}}  = \cfrac{1}{Re}\bigg(\cfrac{\partial ^{2}{\omega}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{\omega}}{\partial {y}^{2}}\bigg)
  
The Stream Function Equation
-------------------------------
Define the stream function

.. math::
  \begin{array}{l}
  u=\cfrac{\partial \psi}{\partial y},\quad v=-\cfrac{\partial \psi}{\partial x}
  \end{array}  

which automatically satisfies the incompressibility conditions

.. math::
  \cfrac{\partial u}{\partial x} +\cfrac{\partial v}{\partial y}=0
  
by substituting

.. math::
  \cfrac{\partial}{\partial x}\cfrac{\partial \psi}{\partial y} -\cfrac{\partial}{\partial y}\cfrac{\partial \psi}{\partial x}=0\\
  
Substituting  

.. math::
  u=\cfrac{\partial \psi}{\partial y},\quad v=-\cfrac{\partial \psi}{\partial x}
  
into the definition of the vorticity  

.. math::
  \omega=\cfrac{\partial v}{\partial x} -\cfrac{\partial u}{\partial y}\\
  
yields  

.. math::
  \cfrac{\partial ^{2}{\psi}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{\psi}}{\partial {y}^{2}}=-\omega
  
The Navier-Stokes equations in vorticity-stream function form are:  

Advection/diffusion equation

.. math::
  \cfrac{\partial {\omega}}{\partial {t}} +\cfrac{\partial \psi}{\partial y}\cfrac{\partial {\omega}}{\partial {x}}-\cfrac{\partial \psi}{\partial x}\cfrac{\partial {\omega}}{\partial {y}}  = \cfrac{1}{Re}\bigg(\cfrac{\partial ^{2}{\omega}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{\omega}}{\partial {y}^{2}}\bigg)\\

Elliptic equation

.. math::
  \cfrac{\partial ^{2}{\psi}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{\psi}}{\partial {y}^{2}}=-\omega
  
Boundary Conditions for the Stream Function
---------------------------------------------
At the right and the left boundary:

.. math::
  \begin{array}{l}
  u=0\Rightarrow \cfrac{\partial \psi}{\partial y} =0\\
  \Rightarrow \psi=Constant
  \end{array}
  
At the top and the bottom boundary:

.. math::
  \begin{array}{l}
  v=0\Rightarrow \cfrac{\partial \psi}{\partial x} =0\\
  \Rightarrow \psi=Constant
  \end{array}
  
Since the boundaries meet, the constant must be the same on all boundaries  

.. math::
  \psi=Constant
  
The normal velocity is zero since the stream function
is a constant on the wall, but the zero tangential
velocity must be enforced

At the right and left boundary:

.. math::
  v=0\Rightarrow \cfrac{\partial \psi}{\partial x} =0\\
  
At the bottom boundary:

.. math::
  u=0\Rightarrow \cfrac{\partial \psi}{\partial y} =0\\
  
At the top boundary:

.. math::
  u=U_{wall}\Rightarrow \cfrac{\partial \psi}{\partial y} =U_{wall}\\
  
The wall vorticity must be found from the streamfunction.
The stream function is constant on the walls.  

At the right and the left boundary:

.. math::
  \cfrac{\partial ^{2}{\psi}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{\psi}}{\partial {y}^{2}}=-\omega
  \Rightarrow \omega_{wall}=-\cfrac{\partial ^{2}{\psi}}{\partial {x}^{2}}
  
Similarly, at the top and the bottom boundary:

.. math::
  \cfrac{\partial ^{2}{\psi}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{\psi}}{\partial {y}^{2}}=-\omega
  \Rightarrow \omega_{wall}=-\cfrac{\partial ^{2}{\psi}}{\partial {y}^{2}}
  
.. figure:: ../images/stream1.png
   :width: 600
   :align: center
   
   boundary
   
Discretizing the Domain
-----------------------------
Uniform mesh (h=constant)

.. figure:: ../images/stream2.png
   :width: 600
   :align: center
   
   mesh   
   
When using FINITE DIFFERENCE approximations,
the values of :math:`f` are stored at discrete points and the
derivatives of the function are approximated using a
Taylor series:   

Start by expressing the value of :math:`f(x+h)` and :math:`f(x-h)`
in terms of :math:`f(x)`

.. figure:: ../images/stream3.png
   :width: 600
   :align: center

Finite Difference Approximations

.. math::
  \begin{array}{l}
  \cfrac{\partial f(x)}{\partial x}=\cfrac{f(x+h)-f(x-h)}{2h}+\cfrac{\partial^{3} f(x)}{\partial x^{3}}\cfrac{h^{2}}{6}+\cdots \\
  \cfrac{\partial ^{2}f(x)}{\partial x^{2}}=\cfrac{f(x+h)-2f(x)+f(x-h)}{h^{2}}+\cfrac{\partial^{4} f(x)}{\partial x^{4}}\cfrac{h^{2}}{12}+\cdots \\
  \cfrac{\partial f(t)}{\partial t}=\cfrac{f(t+\Delta t)-f(\Delta t)}{\Delta t}+\cfrac{\partial^{2} f(t)}{\partial t^{2}}\cfrac{\Delta t}{2}+\cdots \\
  \end{array} 

For a two-dimensional flow discretize the variables on
a two-dimensional grid

.. figure:: ../images/stream4.png
   :width: 600
   :align: center
   
Laplacian  

.. math::
  \begin{align}
  \cfrac{\partial ^{2}{f}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{f}}{\partial {y}^{2}} & = \cfrac{f_{i+1,j}^{n}-2f_{i,j}^{n}+f_{i-1,j}^{n}}{h^{2}} 
  +\cfrac{f_{i,j+1}^{n}-2f_{i,j}^{n}+f_{i,j-1}^{n}}{h^{2}} \\ & = \cfrac{f_{i+1,j}^{n}+f_{i-1,j}^{n}+f_{i,j+1}^{n}+f_{i,j-1}^{n}-4f_{i,j}^{n}}{h^{2}}
  \end{align} 

-

.. math::
  \cfrac{\partial {\omega}}{\partial {t}} +\cfrac{\partial \psi}{\partial y}\cfrac{\partial {\omega}}{\partial {x}}-\cfrac{\partial \psi}{\partial x}\cfrac{\partial {\omega}}{\partial {y}}  = \cfrac{1}{Re}\bigg(\cfrac{\partial ^{2}{\omega}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{\omega}}{\partial {y}^{2}}\bigg)\\

-

.. math::
  \cfrac{\partial {\omega}}{\partial {t}} =-\cfrac{\partial \psi}{\partial y}\cfrac{\partial {\omega}}{\partial {x}}+\cfrac{\partial \psi}{\partial x}\cfrac{\partial {\omega}}{\partial {y}}  + \cfrac{1}{Re}\bigg(\cfrac{\partial ^{2}{\omega}}{\partial {x}^{2}}+\cfrac{\partial ^{2}{\omega}}{\partial {y}^{2}}\bigg)\\  
  
Using these approximations, the vorticity equation becomes:

.. math::
  \begin{align}
  \cfrac{{\omega}_{i,j}^{n+1}-{\omega}_{i,j}^{n}}{\Delta{t}} & =
  -\bigg(\cfrac{\psi_{i,j+1}^{n}-\psi_{i,j-1}^{n}}{2h}\bigg)\bigg(\cfrac{\omega_{i+1,j}^{n}-\omega_{i-1,j}^{n}}{2h}\bigg)
  +\bigg(\cfrac{\psi_{i+1,j}^{n}-\psi_{i-1,j}^{n}}{2h}\bigg)\bigg(\cfrac{\omega_{i,j+1}^{n}-\omega_{i,j-1}^{n}}{2h}\bigg)\\
  &+\cfrac{1}{Re}\bigg(\cfrac{\omega_{i+1,j}^{n}+\omega_{i-1,j}^{n}+\omega_{i,j+1}^{n}+\omega_{i,j-1}^{n}-4\omega_{i,j}^{n}}{h^{2}} \bigg)
  \end{align}
  
The vorticity at the new time is given by:

.. math::
  \begin{array}{l}
  {\omega}_{i,j}^{n+1}={\omega}_{i,j}^{n}+ {\Delta{t}} \bigg[
  -\bigg(\cfrac{\psi_{i,j+1}^{n}-\psi_{i,j-1}^{n}}{2h}\bigg)\bigg(\cfrac{\omega_{i+1,j}^{n}-\omega_{i-1,j}^{n}}{2h}\bigg)\\
  +\bigg(\cfrac{\psi_{i+1,j}^{n}-\psi_{i-1,j}^{n}}{2h}\bigg)\bigg(\cfrac{\omega_{i,j+1}^{n}-\omega_{i,j-1}^{n}}{2h}\bigg)\\
  +\cfrac{1}{Re}\bigg(\cfrac{\omega_{i+1,j}^{n}+\omega_{i-1,j}^{n}+\omega_{i,j+1}^{n}+\omega_{i,j-1}^{n}-4\omega_{i,j}^{n}}{h^{2}} \bigg)\bigg]
  \end{array}
  
The stream function equation is:  

.. math::
  \cfrac{\partial ^{2}\psi}{\partial x^{2}} +
  \cfrac{\partial ^{2}\psi}{\partial y^{2}} =-\omega

-

.. math::
  \cfrac{\psi_{i+1,j}^{n}+\psi_{i-1,j}^{n}+\psi_{i,j+1}^{n}+\psi_{i,j-1}^{n}-4\psi_{i,j}^{n} }{h^{2}}=-\omega_{i,j}^{n}   
  
Discretized Domain

.. figure:: ../images/stream5.png
   :width: 600
   :align: center
   
   Discretized Domain

Discrete Boundary Condition

.. figure:: ../images/stream6.png
   :width: 600
   :align: center
   
   Discrete Boundary Condition
   
.. math::
  \psi_{i,j=2}=\psi_{i,j=1}+\cfrac{\partial\psi_{i,j=1}}{\partial y} h+\cfrac{\partial ^{2}\psi_{i,j=1}}{\partial y^{2}}\cfrac{h^{2}}{2} +O(h^{3})   
  
Using:

.. math::
  \omega_{wall}=-\cfrac{\partial ^{2}\psi_{i,j=1}}{\partial y^{2}} ;\quad U_{wall}=\cfrac{\partial\psi_{i,j=1}}{\partial y}
  
This becomes:

.. math::
  \psi_{i,j=2}=\psi_{i,j=1}+U_{wall}h-\omega_{wall}\cfrac{h^{2}}{2} +O(h^{3})
  
Solving for the wall vorticity:

.. math::
  \omega_{wall}=(\psi_{i,j=1}-\psi_{i,j=2})\cfrac{2}{h^{2}}+U_{wall}\cfrac{2}{h} +O(h)
  
At the bottom wall (j=1) 

.. math::
  \omega_{wall}=(\psi_{i,j=1}-\psi_{i,j=2})\cfrac{2}{h^{2}}+U_{wall}\cfrac{2}{h} +O(h) 
  
Similarly, at the bottom wall (j=ny):

.. math::
  \omega_{wall}=(\psi_{i,j=ny}-\psi_{i,j=ny-1})\cfrac{2}{h^{2}}-U_{wall}\cfrac{2}{h} +O(h) 

-

.. math::
  \psi_{i,j=ny-1}=\psi_{i,j=ny}+\cfrac{\partial\psi_{i,j=ny}}{\partial y} (-h)+\cfrac{\partial ^{2}\psi_{i,j=ny}}{\partial y^{2}}\cfrac{h^{2}}{2} +O(h^{3})   
  
Using:

.. math::
  \omega_{wall}=-\cfrac{\partial ^{2}\psi_{i,j=ny}}{\partial y^{2}} ;\quad U_{wall}=\cfrac{\partial\psi_{i,j=ny}}{\partial y}  
  
This becomes:

.. math::
  \psi_{i,j=ny-1}=\psi_{i,j=ny}+U_{wall}(-h)-\omega_{wall}\cfrac{h^{2}}{2} +O(h^{3})  
  
Solving for the wall vorticity:  
  
.. math::
  \omega_{wall}=(\psi_{i,j=ny}-\psi_{i,j=ny-1})\cfrac{2}{h^{2}} -U_{wall}\cfrac{2}{h}  +O(h)   
  
At the left wall (i=1):

.. math::
  \psi_{i=2,j}=\psi_{i=1,j}+\cfrac{\partial\psi_{i=1,j}}{\partial x} h+\cfrac{\partial ^{2}\psi_{i=1,j}}{\partial x^{2}}\cfrac{h^{2}}{2} +O(h^{3})   
  
Using:

.. math::
  \omega_{wall}=-\cfrac{\partial ^{2}\psi_{i=1,j}}{\partial x^{2}} ;\quad 0=\cfrac{\partial\psi_{i=1,j}}{\partial x}
  
This becomes:

.. math::
  \psi_{i=2,j}=\psi_{i=1,j}+0h-\omega_{wall}\cfrac{h^{2}}{2} +O(h^{3})
  
Solving for the wall vorticity: 
  
.. math::
  \omega_{wall}=(\psi_{i=1,j}-\psi_{i=2,j})\cfrac{2}{h^{2}} +0 \cfrac{2}{h} +O(h)   

At the right wall (i=nx):
  
.. math::
  \psi_{i=nx-1,j}=\psi_{i=nx,j}+\cfrac{\partial\psi_{i=nx,j}}{\partial x} (-h)+\cfrac{\partial ^{2}\psi_{i=nx,j}}{\partial x^{2}}\cfrac{h^{2}}{2} +O(h^{3})   
  
Using:

.. math::
  \omega_{wall}=-\cfrac{\partial ^{2}\psi_{i=nx,j}}{\partial x^{2}} ;\quad 0=\cfrac{\partial\psi_{i=nx,j}}{\partial x}
  
Solving for the wall vorticity: 

.. math::
  \omega_{wall}=(\psi_{i=nx,j}-\psi_{i=nx-1,j})\cfrac{2}{h^{2}} -0\cfrac{2}{h} +O(h)     
  
Solving the elliptic equation:

.. math::
  \cfrac{\psi_{i+1,j}^{n}+\psi_{i-1,j}^{n}+\psi_{i,j+1}^{n}+\psi_{i,j-1}^{n}-4\psi_{i,j}^{n} }{h^{2}}=-\omega_{i,j}^{n} 
  
Rewrite as

.. math::
  \psi_{i,j}^{n}=\cfrac{1}{4} \bigg(\psi_{i+1,j}^{n}+\psi_{i-1,j}^{n}+\psi_{i,j+1}^{n}+\psi_{i,j-1}^{n}+\omega_{i,j}^{n}{h^{2}} \bigg)  

.. figure:: ../images/stream7.png
   :width: 600
   :align: center
   
If the grid points are done in order, half of the points have
already been updated

.. figure:: ../images/stream8.png
   :width: 600
   :align: center
   
Successive Over Relaxation (SOR)   
  
.. math::
  \psi_{i,j}^{n+1}=\beta\cfrac{1}{4} \bigg(\psi_{i+1,j}^{n}+\psi_{i-1,j}^{n}+\psi_{i,j+1}^{n}+\psi_{i,j-1}^{n}+\omega_{i,j}^{n}{h^{2}} \bigg)+(1-\beta)\psi_{i,j}^{n}
  
Limitations on the time step

.. math::
  \cfrac{v\Delta t}{h^{2}}\le\cfrac{1}{4}  \quad\quad \cfrac{(|u|+|v|)\Delta t}{v} \le2
