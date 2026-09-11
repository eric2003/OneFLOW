Sod's shock-tube problem
==================================

J.D. Anderson, Modern Compressible Flow (1984)

#. `The Sod gasdynamics problem as a tool for benchmarking face flux construction in the finite volume method <https://www.sciencedirect.com/science/article/pii/S2468227620303112/>`_
#. `Exact solution of the 1D riemann problem in Newtonian and relativistic hydrodynamics <https://www.scielo.org.mx/pdf/rmfe/v59n1/v59n1a5.pdf>`_
#. `Sod shock tube calculator <https://github.com/ibackus/sod-shocktube/>`_
#. `Sod Shock Tube Problem - CFD Simulation <https://help.sim-flow.com/validation/sod-shock/>`_

Exact solution for the Sod's shock-tube problem
--------------------------------------------------

.. figure:: ../images/sod1.png
   :width: 800
   :align: center
   
   Schematic of stationary and moving shock waves.
   
Stationary normal shock wave
-----------------------------------  
The continuity, momentum, and energy equations of stationary shock waves are, respectively,

.. math::   
  \begin{array}{c}
  \rho_{1} u_{1}=\rho_{2} u_{2}\\
  p_{1}+\rho_{1} {u}_{1}^{2}=p_{2}+\rho_{2} {u}_{2}^{2}\\
  h_{1}+\cfrac{1}{2}{u}_{1}^{2}=h_{2}+\cfrac{1}{2}{u}_{2}^{2}
  \end{array}
  
where 

.. math:: 
  \begin{array}{l}  
  u_{1} =\text{ velocity of the gas ahead of the shock wave, relative to the wave }\\
  u_{2} =\text{ velocity of gas behind the shock wave, relative to the wave }
  \end{array}
  
MOVING NORMAL SHOCK WAVES  
-----------------------------------  

.. figure:: ../images/sod2.png
   :width: 800
   :align: center
   
   Initial conditions in a pressure-drivenshock tube.
   
.. figure:: ../images/sod3.png
   :width: 800
   :align: center
   
   Flow in a shock tube after the diaphragm is broken.
   
The moving normal-shock continuity, momentum, and energy equations, are  
   
.. math::   
  \begin{array}{c}
  \rho_{1} W=\rho_{2} (W-u_{p})\\
  p_{1}+\rho_{1} W^{2}=p_{2}+\rho_{2} (W-u_{p})^{2}\\
  h_{1}+\cfrac{1}{2}W^{2}=h_{2}+\cfrac{1}{2}(W-u_{p})^{2}
  \end{array}   
  
Let us rearrange thcse equations into a more convenient form. 

.. math::
  W-u_{p}=W\cfrac{\rho_{1}}{\rho_{2}}\\
  
-
  
.. math::
  \begin{array}{c}
  p_{1}+\rho_{1} W^{2}=p_{2}+\rho_{2} (W-u_{p})^{2}=p_{2}+\rho_{2} (W\cfrac{\rho_{1}}{\rho_{2}})^{2}\\
  p_{1}+\rho_{1} W^{2}=p_{2}+ \rho_{1}W^{2} (\cfrac{\rho_{1}}{\rho_{2}})\\
  \end{array}  
  
and rearranging,

.. math::
  p_{2}-p_{1}=\rho_{1} W^{2}(1-\cfrac{\rho_{1}}{\rho_{2}})

-
  
.. math::
  \begin{array}{l}
  W^{2}=\cfrac{p_{2}-p_{1}}{\rho_{1}(1-\cfrac{\rho_{1}}{\rho_{2}})}\\
  W^{2}=\cfrac{(p_{2}-p_{1})\rho_{2}}{\rho_{1}(\rho_{2}-\rho_{1})}
  =\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{\rho_{2}}{\rho_{1}}\right)\\
  W^{2}=\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{\rho_{2}}{\rho_{1}}\right)\\
  \end{array}  
  
-
  
.. math::
  \begin{array}{c}
  W= (W-u_{p})\left(\cfrac{\rho_{2}}{\rho_{1}}\right)\\
  W^{2}=(W-u_{p})^{2}\left(\cfrac{\rho_{2}}{\rho_{1}}\right)^{2}\\
  W^{2}=\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{\rho_{2}}{\rho_{1}}\right)\\
  (W-u_{p})^{2}\left(\cfrac{\rho_{2}}{\rho_{1}}\right)=\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}
  \end{array} 
  
-
  
.. math::
  (W-u_{p})^{2}=\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{\rho_{1}}{\rho_{2}}\right)

-
  
.. math::
  h=e+p/\rho  
  
-
  
.. math::
  \begin{align}
  h_{1}+\cfrac{1}{2}W^{2} & = h_{2}+\cfrac{1}{2}(W-u_{p})^{2}\\
  &\Rightarrow e_{1}+\cfrac{p_{1}}{\rho_{1}}
  +\cfrac{1}{2}\left[\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{\rho_{2}}{\rho_{1}}\right)\right]\\ & = e_{2}+\cfrac{p_{2}}{\rho_{2}}
  +\cfrac{1}{2}\left[\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{\rho_{1}}{\rho_{2}}\right)\right]
  \end{align} 
  
The above equation algebraically simplifies to

.. math::  
  \begin{align}
  e_{2}-e_{1}=
  \cfrac{p_{1}}{\rho_{1}}-\cfrac{p_{2}}{\rho_{2}}
  +\cfrac{1}{2}\left[\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{\rho_{2}}{\rho_{1}}\right)\right]
  -\cfrac{1}{2}\left[\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{\rho_{1}}{\rho_{2}}\right)\right]
  \end{align}
  
-
  
.. math::
  \begin{align}
  e_{2}-e_{1}=
  \cfrac{p_{1}}{\rho_{1}}-\cfrac{p_{2}}{\rho_{2}}
  +\cfrac{1}{2}\left[\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{\rho_{2}}{\rho_{1}}-\cfrac{\rho_{1}}{\rho_{2}}\right)\right]\\
  e_{2}-e_{1}=
  \cfrac{p_{1}}{\rho_{1}}-\cfrac{p_{2}}{\rho_{2}}
  +\cfrac{1}{2}\left[\cfrac{p_{2}-p_{1}}{\rho_{2}-\rho_{1}}\left(\cfrac{(\rho_{2})^{2}-(\rho_{1})^{2}}{\rho_{1}\rho_{2}}\right)\right]\\
  e_{2}-e_{1}=
  \cfrac{p_{1}}{\rho_{1}}-\cfrac{p_{2}}{\rho_{2}}
  +\cfrac{1}{2}\left[\cfrac{(p_{2}-p_{1})(\rho_{2}+\rho_{1})}{\rho_{1}\rho_{2}}\right]\\
  \end{align}  
  
-
  
.. math::
  \begin{align}
  e_{2}-e_{1}=
  \cfrac{p_{1}}{\rho_{1}}-\cfrac{p_{2}}{\rho_{2}}
  +\cfrac{1}{2}\left[{(p_{2}-p_{1})\left(\cfrac{1}{\rho_{1}}+\cfrac{1}{\rho_{2}}\right)}\right]\\
  e_{2}-e_{1}=
  p_{1}\left(\cfrac{1}{\rho_{1}}-\cfrac{1}{2}\left(\cfrac{1}{\rho_{1}}+\cfrac{1}{\rho_{2}}\right)\right)
  +p_{2}\left(-\cfrac{1}{\rho_{2}}+\cfrac{1}{2}\left(\cfrac{1}{\rho_{1}}+\cfrac{1}{\rho_{2}}\right)\right)\\
  e_{2}-e_{1}=
  \cfrac{1}{2}p_{1}\left(\cfrac{1}{\rho_{1}}-\cfrac{1}{\rho_{2}}\right)
  -\cfrac{1}{2}p_{2}\left(\cfrac{1}{\rho_{2}}-\cfrac{1}{\rho_{1}}\right)\\
  e_{2}-e_{1}=
  \cfrac{1}{2}(p_{1}+p_{2})\left(\cfrac{1}{\rho_{1}}-\cfrac{1}{\rho_{2}}\right)
  \end{align}    
  
Hugoniot equation
--------------------------  

.. math::
  e_{2}-e_{1}=
  \cfrac{1}{2}(p_{1}+p_{2})\left(\cfrac{1}{\rho_{1}}-\cfrac{1}{\rho_{2}}\right)
  
or  
 
.. math::
  e_{2}-e_{1}=
  \cfrac{1}{2}(p_{1}+p_{2})\left(\nu_{1}-\nu_{2}\right)
  
where :math:`\nu` is specific volume, It is the reciprocal of density :math:`\rho`.
  
Hugoniot equation, and is identically the same form for a stationary shock. In hindsight, this is to be expected; the Hugoniot
equation relates changes of thermodynamic variables across a normal shock wave,
and these are physically independent of whether or not the shock is moving.  

Let us specialize to the case of a calorically perfect gas. In this case, :math:`e=c_{v}T` and :math:`\nu=RT/p`, hence

.. math::
  \begin{array}{l}
  c_{v}=\cfrac{1}{\gamma-1}R\\
  e=c_{v}T=\cfrac{1}{\gamma-1}RT\\
  \end{array}
  
-
  
.. math::  
  \begin{array}{l}
  \cfrac{1}{\gamma-1}R(T_{2}-T_{1})=
    \cfrac{1}{2}(p_{1}+p_{2})\left(\cfrac{RT_{1}}{p_{1}}-\cfrac{RT_{2}}{p_{2}}\right)\\
  \cfrac{1}{\gamma-1}R(\cfrac{T_{2}}{T_{1}}-1)=
    \cfrac{1}{2}(p_{1}+p_{2})\left(\cfrac{R}{p_{1}}-\cfrac{R}{p_{2}}\cfrac{T_{2}}{T_{1}}\right)\\
  \cfrac{1}{\gamma-1}(\cfrac{T_{2}}{T_{1}}-1)=
    \cfrac{1}{2}(p_{1}+p_{2})\left(\cfrac{1}{p_{1}}-\cfrac{1}{p_{2}}\cfrac{T_{2}}{T_{1}}\right)\\
  (\cfrac{T_{2}}{T_{1}}-1)=
    \cfrac{1}{2}({\gamma-1})(p_{1}+p_{2})\left(\cfrac{1}{p_{1}}-\cfrac{1}{p_{2}}\cfrac{T_{2}}{T_{1}}\right)\\
  \cfrac{T_{2}}{T_{1}}=
    1+\cfrac{1}{2}({\gamma-1})(p_{1}+p_{2})\left(\cfrac{1}{p_{1}}-\cfrac{1}{p_{2}}\cfrac{T_{2}}{T_{1}}\right)\\
  \cfrac{T_{2}}{T_{1}}=
    1+\cfrac{1}{2}({\gamma-1})(p_{1}+p_{2})\left(\cfrac{1}{p_{1}}\right)
  -\cfrac{1}{2}({\gamma-1})(p_{1}+p_{2})\left(\cfrac{1}{p_{2}}\cfrac{T_{2}}{T_{1}}\right)
  \end{array}  
  
-
  
.. math:: 
  \begin{array}{l}
  \cfrac{T_{2}}{T_{1}}\left(1+\cfrac{1}{2}({\gamma-1})\left(1+\cfrac{p_{1}}{p_{2}}\right)\right)=
    1+\cfrac{1}{2}({\gamma-1})\left(1+\cfrac{p_{2}}{p_{1}}\right)\\
  \cfrac{T_{2}}{T_{1}}=
    \cfrac{\left(1+\cfrac{1}{2}({\gamma-1})\left(1+\cfrac{p_{2}}{p_{1}}\right)\right)}
  {\left(1+\cfrac{1}{2}({\gamma-1})\left(1+\cfrac{p_{1}}{p_{2}}\right)\right)}\\
  \cfrac{T_{2}}{T_{1}}=
    \cfrac{\left(1+\cfrac{1}{2}({\gamma-1})+\cfrac{1}{2}({\gamma-1})\left(\cfrac{p_{2}}{p_{1}}\right)\right)}
  {\left(1+\cfrac{1}{2}({\gamma-1})+\cfrac{1}{2}({\gamma-1})\left(\cfrac{p_{1}}{p_{2}}\right)\right)}\\
  \cfrac{T_{2}}{T_{1}}=
    \cfrac{\left(\cfrac{1}{2}({\gamma+1})+\cfrac{1}{2}({\gamma-1})\left(\cfrac{p_{2}}{p_{1}}\right)\right)}
  {\left(\cfrac{1}{2}({\gamma+1})+\cfrac{1}{2}({\gamma-1})\left(\cfrac{p_{1}}{p_{2}}\right)\right)}\\
  \cfrac{T_{2}}{T_{1}}=
    \cfrac{\left(({\gamma+1})+({\gamma-1})\left(\cfrac{p_{2}}{p_{1}}\right)\right)}
  {\left(({\gamma+1})+({\gamma-1})\left(\cfrac{p_{1}}{p_{2}}\right)\right)}\\
  \end{array}  
  
-
  
.. math:: 
  \begin{array}{c}
  \cfrac{T_{2}}{T_{1}}=
    \cfrac{\left(\cfrac{\gamma+1}{\gamma-1}+\left(\cfrac{p_{2}}{p_{1}}\right)\right)}
  {\left(\cfrac{\gamma+1}{\gamma-1}+\left(\cfrac{p_{1}}{p_{2}}\right)\right)}\\
  \cfrac{T_{2}}{T_{1}}=
    \cfrac{\left(\cfrac{p_{2}}{p_{1}}\right)\left(\cfrac{\gamma+1}{\gamma-1}+\left(\cfrac{p_{2}}{p_{1}}\right)\right)}
  {\left(\cfrac{p_{2}}{p_{1}}\right)\left(\cfrac{\gamma+1}{\gamma-1}+\left(\cfrac{p_{1}}{p_{2}}\right)\right)}\\
  \cfrac{T_{2}}{T_{1}}=
    \cfrac{\left(\cfrac{p_{2}}{p_{1}}\right)\left(\cfrac{\gamma+1}{\gamma-1}+\left(\cfrac{p_{2}}{p_{1}}\right)\right)}
  {\left(\cfrac{\gamma+1}{\gamma-1}\left(\cfrac{p_{2}}{p_{1}}\right)+1\right)}\\
  \cfrac{T_{2}}{T_{1}}=\cfrac{p_{2}}{p_{1}}
    \left(\cfrac{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}
  {1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}\right)\\
  \end{array}  
  
Similarly,

.. math:: 
  \begin{array}{l}
  p=\rho RT\Rightarrow T=\cfrac{p}{\rho R}\\
  \cfrac{T_{2}}{T_{1}}=\cfrac{p_{2}}{p_{1}}
    \left(\cfrac{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}
  {1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}\right)\\
  \cfrac{T_{2}}{T_{1}}=\cfrac{\rho_{2} RT_{2}}{\rho_{1} RT_{1}}
    \left(\cfrac{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}
  {1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}\right)\\
  1=\cfrac{\rho_{2} }{\rho_{1}}
    \left(\cfrac{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}
  {1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}\right)\\
  \end{array} 

-
  
.. math:: 
  \begin{array}{l}
  \cfrac{{1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}}{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}=\cfrac{\rho_{2} }{\rho_{1}}\\
  \cfrac{\rho_{2} }{\rho_{1}}=\cfrac{{1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}}{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}\\
  \end{array}  
  
Define the moving shock Mach number as  

.. math::
  M_{s}=\cfrac{W}{a_{1}}
  
INCIDENT AND REFLECTED EXPANSION WAVES  
-----------------------------------------------
.. figure:: ../images/sod6.png
   :width: 800
   :align: center
   
   The :math:`C_{+}` and :math:`C_{-}` characteristics for a centered expansion wave(on an xt diagram).
   
.. math::
  \cfrac{a}{a_{4}}=1-\cfrac{\gamma-1}{2}\left(\cfrac{u}{a_{4}} \right)

To obtain the variation of properties in a centered expansion wave as a funclion
of x and t. The equation of any C characteristic is

.. math::
  \cfrac{dx}{dt}=u-a

or, because the characteristic ib a straight line through the origin

.. math::
  x=(u-a)t

-

.. math::
  \begin{array}{l}
  \cfrac{a}{a_{4}}=1-\cfrac{\gamma-1}{2}\left(\cfrac{u}{a_{4}} \right)\\  
  a=a_{4}-\cfrac{\gamma-1}{2}u\\
  u-a=u-a_{4}+\cfrac{\gamma-1}{2}u\\
  u-a=-a_{4}+\cfrac{\gamma+1}{2}u\\
  \end{array}
  
-

.. math::
  \begin{array}{l}
  x=(u-a)t\\
  x=(-a_{4}+\cfrac{\gamma+1}{2}u)t\\
  \cfrac{x}{t}=(-a_{4}+\cfrac{\gamma+1}{2}u)\\
  \cfrac{x}{t}+a_{4}=\cfrac{\gamma+1}{2}u\\
  \cfrac{2}{\gamma+1}\left( {\cfrac{x}{t}+a_{4}}\right)=u\\
  u=\cfrac{2}{\gamma+1}\left( {\cfrac{x}{t}+a_{4}}\right)\\
  \end{array}


SHOCK TUBE RELATIONS
------------------------------------
.. figure:: ../images/sod3.png
   :width: 800
   :align: center
   
   Flow in a shock tube after the diaphragm is broken.
   
.. figure:: ../images/sod5.png
   :width: 600
   :align: center   
   
   Schematic shock tube problem with pressure distribution for pre- and post-diaphragm removal.
   
.. figure:: ../images/sod4.png
   :width: 600
   :align: center   
   
   Diagram of the shock, expansion waves and contact surface. 
   
   
.. math::   
  \begin{array}{l}
  p_{3}=p_{2}\\
  u_{3}=u_{2}=u_{p}
  \end{array}   

-

.. math::
  u_{p}=u_{2}=\cfrac{a_{1}}{\gamma_{1}}\left(\cfrac{p_{2}}{p_{1}}-1\right)\left(\cfrac{\cfrac{2\gamma_{1}}{\gamma_{1}+1}}{\cfrac{p_{2}}{p_{1}}+\cfrac{\gamma_{1}-1}{\gamma_{1}+1}}\right)^{1/2}
  
between the head and tail of the expansion wave

.. math::
  \cfrac{p_{3}}{p_{4}}=\left[1-\cfrac{\gamma_{4}}{2}\left(\cfrac{u_{3}}{u_{4}}\right)\right]^{\cfrac{2\gamma_{4}}{\gamma_{4}-1}}
  
:math:`\cfrac{p2}{p1}`:

.. math::  
  \cfrac{p_{4}}{p_{1}}=\cfrac{p_{2}}{p_{1}}
  \left\{1-\cfrac{(\gamma_{4}-1)\left(\cfrac{a_{1}}{a_{4}}\right)\left(\cfrac{p_{2}}{p_{1}}-1\right)}
  {\sqrt{2\gamma_{1}\left[2\gamma_{1}+(\gamma_{1}+1)\left(\cfrac{p_{2}}{p_{1}}-1\right)\right]}}\right\}
  ^{\cfrac{-2\gamma_{4}}{\gamma_{4}-1}}  

The analysis of the flow of a calorically perfect gas in a shock tube is now
straightforward. For a given diaphragm pressure ratio :math:`p4/p1`:

1. Calculate :math:`p2/p1`. This defines the strength of the incident shock wave.

.. math::  
  \cfrac{p_{4}}{p_{1}}=\cfrac{p_{2}}{p_{1}}
  \left\{1-\cfrac{(\gamma_{4}-1)\left(\cfrac{a_{1}}{a_{4}}\right)\left(\cfrac{p_{2}}{p_{1}}-1\right)}
  {\sqrt{2\gamma_{1}\left[2\gamma_{1}+(\gamma_{1}+1)\left(\cfrac{p_{2}}{p_{1}}-1\right)\right]}}\right\}
  ^{\cfrac{-2\gamma_{4}}{\gamma_{4}-1}} 
  
2. Calculate all other incident shock properties:

.. math::  
  \cfrac{T_{2}}{T_{1}}=\cfrac{p_{2}}{p_{1}}
    \left(\cfrac{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}
  {1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}\right)\\
  
-
  
.. math::   
  \cfrac{\rho_{2} }{\rho_{1}}=\cfrac{{1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}}{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}\\

-
  
.. math::   
  W=a_{1}\sqrt{\cfrac{\gamma+1}{2\gamma}\left(\cfrac{p_{2}}{p_{1}}-1\right)+1}  
  
-

.. math::
  u_{p}=u_{2}=\cfrac{a_{1}}{\gamma_{1}}\left(\cfrac{p_{2}}{p_{1}}-1\right)\left(\cfrac{\cfrac{2\gamma_{1}}{\gamma_{1}+1}}{\cfrac{p_{2}}{p_{1}}+\cfrac{\gamma_{1}-1}{\gamma_{1}+1}}\right)^{1/2}
  
-

.. math::   
  \begin{array}{l}
  p_{3}=p_{2}\\
  u_{3}=u_{2}=u_{p}
  \end{array}   
    

3. Calculate :math:`p_{3}/p_{4}=(p_{3}/p_{1})/(p_{4}/p_{1})=(p_{2}/p_{1})/(p_{4}/p_{1})`. This defines the strength of the incident expansion wave.
4. All other thermodynamic properties immediately behind the expansion wave can be found from the isentropic relations

.. math::
  \cfrac{p_{3}}{p_{4}}=\left(\cfrac{\rho_{3}}{\rho_{4}}\right)^{\gamma}=\left(\cfrac{T_{3}}{T_{4}}\right)^{\cfrac{\gamma}{\gamma -1}}
  
5. Calculate the local properties inside the expansion wave:

.. math::
  \cfrac{a}{a_{4}}=1-\cfrac{\gamma-1}{2}\left(\cfrac{u}{a_{4}} \right)
  
-
  
.. math::  
  u=\cfrac{2}{\gamma+1}\left(a_{4}+\cfrac{x}{t} \right )   
  
Equation holds for the region between the head and tail of the centered expansion wave  
i.e., :math:`-a_{4}\le \cfrac{x}{t}\le u_{3}-a_{3}`

.. math::  
  \begin{array}{l}
  x_{\text{head}} = x_{\text{diaphragm}}+(u_{4}-a_{4})\times t\\
  x_{\text{tail}} = x_{\text{diaphragm}}+(u_{3}-a_{3})\times t\\
  x_{\text{head}}-x_{\text{tail}}=\left\{(u_{4}-a_{4})-(u_{3}-a_{3})\right\}\times t\\
  \end{array}
  
-
  
.. math::
  \begin{array}{l}
  u=\cfrac{2}{\gamma+1}\left(a_{4}+\cfrac{x}{t} \right )\\   
  u=\cfrac{2a_{4}}{\gamma+1}+\cfrac{2}{\gamma+1}\cfrac{x}{t}\\  
  u_{\text{head}}=\cfrac{2a_{4}}{\gamma+1}+\cfrac{2}{\gamma+1}\cfrac{x_{\text{head}}}{t}\\ 
  u-u_{\text{head}}=\cfrac{2}{\gamma+1}\cfrac{x-x_{\text{head}}}{t}\\ 
  u-u_{\text{head}}\propto x-x_{\text{head}}\\ 
  \end{array} 
  
-
  
.. math::
  \begin{array}{l}
  u_{\text{head}}=\cfrac{2a_{4}}{\gamma+1}+\cfrac{2}{\gamma+1}\cfrac{x_{\text{head}}}{t}\\ 
  u_{\text{tail}}=\cfrac{2a_{4}}{\gamma+1}+\cfrac{2}{\gamma+1}\cfrac{x_{\text{tail}}}{t}\\ 
  u_{\text{tail}}-u_{\text{head}}=\cfrac{2}{\gamma+1}\cfrac{x_{\text{tail}}-x_{\text{head}}}{t}\\ 
  u-u_{\text{head}}=\cfrac{2}{\gamma+1}\cfrac{x-x_{\text{head}}}{t}\\ 
  \cfrac{u-u_{\text{head}}}{u_{\text{tail}}-u_{\text{head}}} =\cfrac{x-x_{\text{head}}}{x_{\text{tail}}-x_{\text{head}}} \\
  \end{array}  

Code implementation details
--------------------------------
Set initial states (non-dimensional).

.. math:: 
  \begin{array}{l}
    p_{4} = 1.0;\\
    r_{4} = 1.0;\\
    u_{4} = 0.0;\\
  \\
    p_{1} = 0.1;\\
    r_{1} = 0.125;\\
    u_{1} = 0.0;\\
  \end{array}
  
Set dimensions of shocktube. 

.. math:: 
  \begin{array}{c}
    x_{l} = 0.0;\\
    x_{r} = 1.0;\\
    x_{d} = 0.5;\\
  \end{array}
  
Calc acoustic velocities. 

.. math:: 
  \begin{array}{c}
    a_{1} = \sqrt{ \gamma \cfrac{p_{1}}{\rho_{1}}}\\
    a_{4} = \sqrt{ \gamma \cfrac{p_{4}}{\rho_{4}}}\\
  \end{array}
  
Use a Newton-secant iteration to compute p2p1.  

.. math::  
  \cfrac{p_{4}}{p_{1}}=\cfrac{p_{2}}{p_{1}}
  \left\{1-\cfrac{(\gamma_{4}-1)\left(\cfrac{a_{1}}{a_{4}}\right)\left(\cfrac{p_{2}}{p_{1}}-1\right)}
  {\sqrt{2\gamma_{1}\left[2\gamma_{1}+(\gamma_{1}+1)\left(\cfrac{p_{2}}{p_{1}}-1\right)\right]}}\right\}
  ^{\cfrac{-2\gamma_{4}}{\gamma_{4}-1}} 
  
Calculate all other incident shock properties:

.. math::  
  \cfrac{T_{2}}{T_{1}}=\cfrac{p_{2}}{p_{1}}
    \left(\cfrac{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}
  {1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}\right)\\  
  
-
  
.. math::   
  \cfrac{\rho_{2} }{\rho_{1}}=\cfrac{{1+\cfrac{\gamma+1}{\gamma-1}\cfrac{p_{2}}{p_{1}}}}{\cfrac{\gamma+1}{\gamma-1}+\cfrac{p_{2}}{p_{1}}}\\  
  
Calculate shock-wave speed.  

.. math::   
  W=a_{1}\sqrt{\cfrac{\gamma+1}{2\gamma}\left(\cfrac{p_{2}}{p_{1}}-1\right)+1} 
  
Calculate Shock location.  

.. math::   
  x_{\text{shock}} = x_{\text{diaphragm}}+W\times t;\\
  
Calculate State 2. 

.. math::   
  \begin{array}{l}
  p_{2}=\cfrac{p_{2}}{p_{1}}p_{1}\\
  {\rho}_{2}=\cfrac{{\rho}_{2}}{{\rho}_{1}}{\rho}_{1}\\
  a_{2}=\sqrt{\gamma \cfrac{p_{2}}{\rho_{2}} }
  \end{array}
  
Calculate State 3.   

.. math::   
  p_{3} = p_{2}\\
  
Isentropic between 3 and 4.  

.. math::  
  \begin{array}{l}
  \cfrac{\rho_{3}}{\rho_{4}}=\left(\cfrac{p_{3}}{p_{4}}\right)^{\cfrac{1}{\gamma}}\\
  a_{3}=\sqrt{\gamma \cfrac{p_{3}}{\rho_{3}} }
  \end{array} 
  
Calculate the speed of contact discontinuity.  

.. math::  
  \begin{array}{l}
  u_{p}=\cfrac{a_{1}}{\gamma_{1}}\left(\cfrac{p_{2}}{p_{1}}-1\right)\left(\cfrac{\cfrac{2\gamma_{1}}{\gamma_{1}+1}}{\cfrac{p_{2}}{p_{1}}+\cfrac{\gamma_{1}-1}{\gamma_{1}+1}}\right)^{\cfrac{1}{2}}\\
  u_{2}=u_{p}\\
  u_{3}=u_{p}\\
  \end{array}
  
Calculate Mach numbers. 

.. math::  
  \begin{array}{l}
  M_{1}=\cfrac{u_{1}}{a_{1}}\\
  M_{2}=\cfrac{u_{2}}{a_{2}}\\
  M_{3}=\cfrac{u_{3}}{a_{3}}\\
  M_{4}=\cfrac{u_{4}}{a_{4}}\\
  \end{array} 
  
Calculate the location of contact discontinuity.  

.. math::  
  x_{\text{contact discontinuity}}=x_{\text{diaphragm}}+u_{p}\times t\\
  
Calculate the location of expansion region.

.. math::  
  \begin{align}
  x_{\text{expansion head}} & = x_{\text{diaphragm}}+(u_{4}-a_{4})\times t\\
  x_{\text{expansion   tail}} & = x_{\text{diaphragm}}+(u_{3}-a_{3})\times t\\
  \end{align}
  
Calculate all other expansion wave properties: 
 
.. math::
  x_{\text{expansion}}=x_{\text{expansion head}}+dx_{\text{expansion}}*\cfrac{i}{n_{\text{ expansion points}}}
  
Let   

.. math::
  x=x_{\text{expansion}}
  
then  

.. math::
  \begin{array}{l}
  u(x)=u_{\text{head}}+(u_{\text{tail}}-u_{\text{head}})\cfrac{x-x_{\text{head}}}{x_{\text{tail}}-x_{\text{head}}} \\
  u(x)=u_{4}+(u_{3}-u_{4})\cfrac{x-x_{\text{head}}}{x_{\text{tail}}-x_{\text{head}}} \\
  u_{4}=0\\
  u(x)=(u_{3})\cfrac{x-x_{\text{head}}}{x_{\text{tail}}-x_{\text{head}}} \\
  \end{array}

-
  
.. math::  
  \begin{array}{l}
  \cfrac{p(x)}{p_{4}}=\left[1-\cfrac{\gamma-1}{2}\left(\cfrac{u(x)}{a_{4}} \right)\right]^{\cfrac{2\gamma }{\gamma-1} }\\
  \cfrac{\rho(x)}{\rho_{4}}=\left[1-\cfrac{\gamma-1}{2}\left(\cfrac{u(x)}{a_{4}} \right)\right]^{\cfrac{2 }{\gamma-1} }\\
  M(x)=\cfrac{u(x)}{\sqrt{\gamma \cfrac{p(x)}{\rho (x)} } } 
  \end{array}
 
Newton-secant method
----------------------------
Use a Newton-secant iteration to compute p2p1.

.. math::  
  \cfrac{p_{4}}{p_{1}}=\cfrac{p_{2}}{p_{1}}
  \left\{1-\cfrac{(\gamma_{4}-1)\left(\cfrac{a_{1}}{a_{4}}\right)\left(\cfrac{p_{2}}{p_{1}}-1\right)}
  {\sqrt{2\gamma_{1}\left[2\gamma_{1}+(\gamma_{1}+1)\left(\cfrac{p_{2}}{p_{1}}-1\right)\right]}}\right\}
  ^{\cfrac{-2\gamma_{4}}{\gamma_{4}-1}} 

Let :math:`x=p_{2}/p_{1}`
Initialize x for starting guess 

.. math:: 
  \begin{array}{l}
  x=0.9\cfrac{p_{4}}{p_{1}}\\
  \end{array}
  
-

.. math:: 
  f=\cfrac{p_{4}}{p_{1}}-x\times\left\{1-\cfrac{(\gamma-1)(\cfrac{a_{1}}{a_{4}} )(x-1)}{\sqrt{2\gamma [2\gamma+(\gamma +1)(x-1)]}} \right \}
 ^{-\cfrac{2\gamma }{\gamma -1} }
 
Perturb :math:`x`

.. math:: 
  \hat{x}=0.95x
  
Begin iteration

.. math:: 
  \begin{array}{l}
  \text{iter} = 0\\
  \text{itmax} = 20\\
  \end{array}
  
.. math::   
  \begin{align}
    \text{while True:}\\
        iter & = iter + 1\\
  \hat{f}&=\cfrac{p_{4}}{p_{1}}-\hat{x}\times\left\{1-\cfrac{(\gamma-1)(\cfrac{a_{1}}{a_{4}} )(\hat{x}-1)}{\sqrt{2\gamma [2\gamma+(\gamma +1)(\hat{x}-1)]}} \right \}
  ^{-\cfrac{2\gamma }{\gamma -1} }\\
  \text{if}&\text{ abs}( \hat{f} ) \le \text{tol or iter} \ge \text{itmax:}\\
  \quad &\quad\quad\text{   break}\\
       \widetilde{x}&= \hat{x}- \hat{f} * ( \hat{x}- x) / ( \hat{f} - f);\\
        x&= \hat{x};\\
        f&= \hat{f};\\
        \hat{x}&= \widetilde{x};\\
  \end{align}  
  
