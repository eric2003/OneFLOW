Arbitrary Lagrangian-Eulerian
==================================
 
Arbitrary Lagrangian-Eulerian (ALE)
---------------------------------------------------------------------------
#. `Arbitrary Lagrangian-Eulerian (ALE) and Computational Fluid Dynamics (CFD) <https://2021.help.altair.com/2021/hwsolvers/rad/topics/solvers/rad/arbitrary_lagrangian_eulerian_computational_fluid_dynamics_c.htm>`_
#. `ALE formula <https://max.book118.com/html/2016/1212/71106288.shtm>`_
#. `Arbitrary Lagrangian-Eulerian Methods <https://ww2.lacan.upc.edu/scientificPublications/files/pdfs/2017-ECM-DHPR-blanc.pdf>`_
#. `Arbitrary Lagrangian-Eulerian Methods(Wiley Online Library) <https://onlinelibrary.wiley.com/doi/full/10.1002/0470091355.ecm009>`_
#. `Nonlinear finite elements for continua and structures, Second edition <https://www.wiley.com/en-nz/Nonlinear+Finite+Elements+for+Continua+and+Structures%2C+2nd+Edition-p-9781118632703>`_
#. `image-to-latex <https://github.com/kingyiusuen/image-to-latex/>`_


Eulerian and Lagrangian Descriptions
---------------------------------------------------------------------------
In the mathematics and continuum mechanics literature (cf. Marsden and Hughes, 1983),
different symbols are often used for the same field when it is expressed in terms of different
independent variables, that is, when the description is Eulerian or Lagrangian. In this
convention, the function which in an Eulerian description is :math:`f(\mathbf{x},t)` is denoted by :math:`F(\mathbf{X},t)` in a
Lagrangian description. The two functions are related by

.. math::
  F(\mathbf{X},t)=f(\mathbf{\Phi}(\mathbf{X},t),t) \quad or \quad F=f\circ  \mathbf{\Phi}

Lagrangian and Eulerian viewpoints
---------------------------------------------------------------------------

Two domains are commonly used in continuum mechanics: the material domain :math:`R_{X} \subset \mathbb{R}^{n_{\text {sd }}}`,
with :math:`n_{\text {sd }}` spatial dimensions, made up of material particles :math:`\mathbf{X}`, and the spatial domain :math:`R_{\mathbf{x}}`,
consisting of spatial points :math:`\mathbf{x}`.

.. image:: images/ale1.png
   :width: 600

.. math::
  \begin{aligned}
    \boldsymbol{\varphi}: R_{\mathbf{X}} \times\left[t_{0}, t_{\text {final }}]\right. & \longrightarrow R_{\mathbf{x}} \times\left[t_{0}, t_{\text {final }}]\right. \\
    (\mathbf{X}, t) & \longmapsto \boldsymbol{\varphi}(\mathbf{X}, t)=(\mathbf{x}, t)
  \end{aligned}
  
which allows us to link :math:`\mathbf{X}` and :math:`\mathbf{x}` in time by the law of motion, namely

.. math::
  \mathbf{x}=\mathbf{x}(\mathbf{X}, t), \quad t=t
  
for two dimensions:

.. math:: 
  \mathbf{x}=\begin{bmatrix}
   x_{1}\\ x_{2}
  \end{bmatrix} \quad or \quad\mathbf{x}=\begin{bmatrix}
   x\\y
  \end{bmatrix} 
  
-

.. math:: 
  \mathbf{X}=\begin{bmatrix}
  X_{1}\\ X_{2}
  \end{bmatrix} \quad or \quad\mathbf{X}=\begin{bmatrix}
  X\\Y
  \end{bmatrix}
  
-

.. math:: 
  \begin{align}
  x_{1}=x_{1}(X_{1},X_{2},t)=\varphi_{1}(X_{1},X_{2},t)\\
  x_{2}=x_{2}(X_{1},X_{2},t)=\varphi_{2}(X_{1},X_{2},t)
  \end{align} 

The matrix

.. math::
  J(u,v)=\cfrac{\partial (x,y)}{\partial (u,v)} =\cfrac{\partial (x(u,v),y(u,v))}{\partial (u,v)}=\begin{bmatrix}
  \cfrac{\partial x(u,v)}{\partial u} & \cfrac{\partial x(u,v)}{\partial v}\\
  \cfrac{\partial y(u,v)}{\partial u}&\cfrac{\partial y(u,v)}{\partial v}
  \end{bmatrix}

-
  
.. math::
  J(u,v)=\cfrac{\partial (x,y)}{\partial (u,v)} =\begin{bmatrix}
  \cfrac{\partial x}{\partial u} & \cfrac{\partial x}{\partial v}\\
  \cfrac{\partial y}{\partial u}&\cfrac{\partial y}{\partial v}
  \end{bmatrix} 
  
-
  
.. math::
  J(X_{1},X_{2})=\cfrac{\partial (x_{1},x_{2})}{\partial (X_{1},X_{2})} =\begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial X_{1}} & \cfrac{\partial x_{1}}{\partial X_{2}}\\
  \cfrac{\partial x_{2}}{\partial X_{1}}&\cfrac{\partial x_{2}}{\partial X_{2}}
  \end{bmatrix} 
  
-
  
.. math::  
  J(X,Y)=\cfrac{\partial (x,y)}{\partial (X,Y)} =\begin{bmatrix}
  \cfrac{\partial x}{\partial X} & \cfrac{\partial x}{\partial Y}\\
  \cfrac{\partial y}{\partial X} & \cfrac{\partial y}{\partial Y}
  \end{bmatrix}
  
-
  
.. math:: 
  J(\mathbf{X})=\cfrac{\partial (\mathbf{x})}{\partial (\mathbf{X})} =
  \cfrac{\partial \mathbf{x}}{\partial \mathbf{X}} =\begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial X_{1}} & \cfrac{\partial x_{1}}{\partial X_{2}}\\
  \cfrac{\partial x_{2}}{\partial X_{1}} & \cfrac{\partial x_{2}}{\partial X_{2}}
  \end{bmatrix} 
  
-
  
.. math:: 
  (\mathbf{x},t)={\stackrel\frown{\mathbf{x}} }=\begin{bmatrix}
  x_{1}\\x_{2}\\t
 \end{bmatrix}
 
-
  
.. math::
  \begin{align}
  x_{1}=x_{1}(X_{1},X_{2},t)=\varphi_{1}(X_{1},X_{2},t)\\
  x_{2}=x_{2}(X_{1},X_{2},t)=\varphi_{2}(X_{1},X_{2},t)\\
  t=\varphi_{3}(X_{1},X_{2},t)\\
  \end{align} 
  
-
  
.. math:: 
  J(X_{1},X_{2},t)=\cfrac{\partial (x_{1},x_{2},t)}{\partial (X_{1},X_{2},t)} =\begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial X_{1}} & \cfrac{\partial x_{1}}{\partial X_{2}}& \cfrac{\partial x_{1}}{\partial t}\\
  \cfrac{\partial x_{2}}{\partial X_{1}} & \cfrac{\partial x_{2}}{\partial X_{2}}& \cfrac{\partial x_{2}}{\partial t}\\
  \cfrac{\partial t}{\partial X_{1}} & \cfrac{\partial t}{\partial X_{2}}& \cfrac{\partial t}{\partial t}\\
  \end{bmatrix} 
  
-
  
.. math:: 
  J(X_{1},X_{2},t)=\cfrac{\partial (x_{1},x_{2},t)}{\partial (X_{1},X_{2},t)} =\begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial X_{1}} & \cfrac{\partial x_{1}}{\partial X_{2}}& v_{1}\\
  \cfrac{\partial x_{2}}{\partial X_{1}} & \cfrac{\partial x_{2}}{\partial X_{2}}& v_{2}\\
  0 & 0 & 1\\
  \end{bmatrix} 
  
which explicitly states the particular nature of :math:`\boldsymbol{\varphi}` : first, the spatial coordinates :math:`\mathbf{x}` depend both
on the material particle, :math:`\mathbf{X}` , and time :math:`t` , and, second, physical time is measured by the same
variable :math:`t` in both material and spatial domains. For every fixed instant :math:`t`, the mapping :math:`\boldsymbol{\varphi}` defines
a configuration in the spatial domain. It is convenient to employ a matrix representation for
the gradient of :math:`\boldsymbol{\varphi}`,

.. math::
  \frac{\partial ({\varphi}_{1},{\varphi}_{2},t)}{\partial(X_{1},X_{2}, t)}=
  \frac{\partial \boldsymbol{\varphi}}{\partial(\mathbf{X}, t)}=\left(\begin{array}{cc}
  \frac{\partial \mathbf{x}}{\partial \mathbf{X}} & \mathbf{v} \\
  \mathbf{0}^{\mathrm{T}} & 1
  \end{array}\right)  
  
where :math:`\mathbf{0}^{\mathrm{T}}` is a null row-vector and the material velocity :math:`\mathbf{v}` is

.. math::
  \mathbf{v}(\mathbf{X}, t)=\left.\frac{\partial \mathbf{x}(\mathbf{X}, t)}{\partial t}\right|_{\mathbf{X}}  
  
-

.. math::  
  \mathbf{v}=\begin{bmatrix}
   v_{1}\\v_{2}
  \end{bmatrix}=
  \begin{bmatrix}
   \cfrac{\partial x_{1}(X_{1},X_{2},t)}{\partial t}\Bigg|_{\mathbf{X}} \\
   \cfrac{\partial x_{2}(X_{1},X_{2},t)}{\partial t}\Bigg|_{\mathbf{X}} 
  \end{bmatrix}  
  
ALE kinematical description
---------------------------------------------------------------------------

.. image:: images/ale2.png
   :width: 500
   
In the ALE description of motion, neither the material configuration :math:`R_{\mathbf{X}}` nor the spatial
configuration :math:`R_{\mathbf{x}}` is taken as the reference. Thus, a third domain is needed: the referential
configuration :math:`R_{\boldsymbol{\chi}}` where reference coordinates :math:`\boldsymbol{\chi}` are introduced to identify the grid points.
Figure above shows these domains and the one-to-one transformations relating the configurations.
The referential domain :math:`R_{\boldsymbol{\chi}}` is mapped into the material and spatial domains by :math:`\boldsymbol{\Psi}` and :math:`\boldsymbol{\varphi}`
respectively. The particle motion :math:`\boldsymbol{\varphi}`  may then be expressed as :math:`\boldsymbol{\varphi}=\boldsymbol{\Phi} \circ \boldsymbol{\Psi}^{-1}` , clearly showing
that, of course, the three mappings :math:`\boldsymbol{\Psi}`, :math:`\boldsymbol{\Phi}` , and :math:`\boldsymbol{\varphi}` are not independent. 
  
.. math::
  \begin{aligned}
  \boldsymbol{\Phi}: R_{\boldsymbol{\chi}} \times\left[t_0, t_{\text {final }}]\right. & \longrightarrow R_{\mathbf{x}} \times\left[t_0, t_{\text {final }}]\right. \\
  (\boldsymbol{\chi}, t) & \longmapsto \boldsymbol{\Phi}(\boldsymbol{\chi}, t)=(\mathbf{x}, t)
  \end{aligned}
  
-  
  
.. math::
  \begin{align}
  x_{1}=x_{1}(\chi_{1},\chi_{2},t)=\phi_{1}(\chi_{1},\chi_{2},t)\\
  x_{2}=x_{2}(\chi_{1},\chi_{2},t)=\phi_{2}(\chi_{1},\chi_{2},t)\\
  \end{align} 

-
  
.. math::  
  J(\chi_{1},\chi_{2})=\cfrac{\partial (x_{1},x_{2})}{\partial (\chi_{1},\chi_{2})} =\begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial \chi_{1}} & \cfrac{\partial x_{1}}{\partial \chi_{2}}\\
  \cfrac{\partial x_{2}}{\partial \chi_{1}} & \cfrac{\partial x_{2}}{\partial \chi_{2}}
  \end{bmatrix}   
  
-
  
.. math::  
  J(\chi_{1},\chi_{2},t)=\cfrac{\partial (x_{1},x_{2},t)}{\partial (\chi_{1},\chi_{2},t)} =\begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial \chi_{1}} & \cfrac{\partial x_{1}}{\partial \chi_{2}}& \cfrac{\partial x_{1}}{\partial t}\\
  \cfrac{\partial x_{2}}{\partial \chi_{1}} & \cfrac{\partial x_{2}}{\partial \chi_{2}}& \cfrac{\partial x_{2}}{\partial t}\\
  \cfrac{\partial t}{\partial \chi_{1}} & \cfrac{\partial t}{\partial \chi_{2}}& \cfrac{\partial t}{\partial t}\\
  \end{bmatrix}

-
  
.. math:: 
   J(\chi_{1},\chi_{2},t)=\cfrac{\partial (x_{1},x_{2},t)}{\partial (\chi_{1},\chi_{2},t)} =\begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial \chi_{1}} & \cfrac{\partial x_{1}}{\partial \chi_{2}}& \hat{v}_{1}\\
  \cfrac{\partial x_{2}}{\partial \chi_{1}} & \cfrac{\partial x_{2}}{\partial \chi_{2}}& \hat{v}_{2}\\
  0 & 0& 1\\
  \end{bmatrix}

-
  
.. math:: 
  \mathbf{\hat{v}}=\begin{bmatrix}
   \hat{v}_{1}\\\hat{v}_{2}
  \end{bmatrix}=\begin{bmatrix}
   \hat{v}_{1}(\chi_{1},\chi_{2},t)\\\hat{v}_{2}(\chi_{1},\chi_{2},t)
  \end{bmatrix}=
  \begin{bmatrix}
   \cfrac{\partial x_{1}(\chi_{1},\chi_{2},t)}{\partial t}\Bigg|_{\boldsymbol{\chi}} \\
   \cfrac{\partial x_{2}(\chi_{1},\chi_{2},t)}{\partial t}\Bigg|_{\boldsymbol{\chi}} 
  \end{bmatrix}  
  
and its gradient is

.. math::
  \frac{\partial \boldsymbol{\Phi}}{\partial(\boldsymbol{\chi}, t)}=\left(\begin{array}{ll}
  \frac{\partial \mathbf{x}}{\partial \boldsymbol{\chi}} & \hat{\mathbf{v}} \\
  \mathbf{0}^{\mathrm{T}} & 1
  \end{array}\right)
  
where now, the mesh velocity

.. math::
  \hat{\mathbf{v}}(\boldsymbol{\chi}, t)=\left.\frac{\partial \mathbf{x}(\boldsymbol{\chi}, t)}{\partial t}\right|_{\boldsymbol{\chi}}
  
is involved. Note that both the material and the mesh move with respect to the laboratory.
Thus, the corresponding material and mesh velocities have been defined by deriving the
equations of material motion and mesh motion respectively with respect to time.

Finally, regarding :math:`\boldsymbol{\Psi}`, it is convenient to represent directly its inverse :math:`\boldsymbol{\Psi}^{-1}`,

.. math::
  \begin{aligned}
  \boldsymbol{\Psi}^{-1}: R_{\mathbf{X}} \times\left[t_0, t_{\text {final }}]\right. & \longrightarrow R_{\boldsymbol{\chi}} \times\left[t_0, t_{\text {final }}]\right. \\
  (\mathbf{X}, t) & \longmapsto \boldsymbol{\Psi}^{-1}(\mathbf{X}, t)=(\boldsymbol{\chi}, t)
  \end{aligned}
  
-
  
.. math::
  \begin{align}
  \chi_{1}=\chi_{1}(X_{1},X_{2},t)=\psi^{-1}_{1}(X_{1},X_{2},t)\\
  \chi_{2}=\chi_{2}(X_{1},X_{2},t)=\psi^{-1}_{2}(X_{1},X_{2},t)\\
  \end{align} 
  
-
  
.. math::
  J(X_{1},X_{2},t)=\cfrac{\partial (\chi_{1},\chi_{2},t)}{\partial (X_{1},X_{2},t)} =\begin{bmatrix}
  \cfrac{\partial \chi_{1}}{\partial X_{1}} & \cfrac{\partial \chi_{1}}{\partial X_{2}}& \cfrac{\partial \chi_{1}}{\partial t}\\
  \cfrac{\partial \chi_{2}}{\partial X_{1}} & \cfrac{\partial \chi_{2}}{\partial X_{2}}& \cfrac{\partial \chi_{2}}{\partial t}\\
  \cfrac{\partial t}{\partial X_{1}} & \cfrac{\partial t}{\partial X_{2}}& \cfrac{\partial t}{\partial t}\\
  \end{bmatrix} 
  
-
  
.. math::  
  J(X_{1},X_{2},t)=\cfrac{\partial (\chi_{1},\chi_{2},t)}{\partial (X_{1},X_{2},t)} =\begin{bmatrix}
  \cfrac{\partial \chi_{1}}{\partial X_{1}} & \cfrac{\partial \chi_{1}}{\partial X_{2}}& w_{1}\\
  \cfrac{\partial \chi_{2}}{\partial X_{1}} & \cfrac{\partial \chi_{2}}{\partial X_{2}}& w_{2}\\
  0 & 0& 1\\
  \end{bmatrix}
  
-
  
.. math::  
  \mathbf{w}=\begin{bmatrix}
  {w}_{1}\\{w}_{2}
  \end{bmatrix}=\begin{bmatrix}
  {w}_{1}(X_{1},X_{2},t)\\{w}_{2}(X_{1},X_{2},t)
  \end{bmatrix}=
  \begin{bmatrix}
  \left.\cfrac{\partial \chi_{1}(X_{1},X_{2},t)}{\partial t}\right|_{\mathbf{X}} \\
  \left.\cfrac{\partial \chi_{2}(X_{1},X_{2},t)}{\partial t}\right|_{\mathbf{X}} 
  \end{bmatrix}  
  
-
  
.. math::    
  {\mathbf{w}}(\mathbf{X}, t)=\left.\frac{\partial \boldsymbol{\chi}(\mathbf{X}, t)}{\partial t}\right|_{\mathbf{X}}  
  
and its gradient is

.. math::
  \frac{\partial \boldsymbol{\Psi}^{-1}}{\partial(\mathbf{X}, t)}=\left(\begin{array}{cc}
  \frac{\partial \boldsymbol{\chi}}{\partial \mathbf{X}} & \mathbf{w} \\
  \mathbf{0}^{\mathrm{T}} & 1
  \end{array}\right)
  
where the velocity :math:`\mathbf{w}` is defined as  

.. math::    
  {\mathbf{w}}(\mathbf{X}, t)=\left.\frac{\partial \boldsymbol{\chi}(\mathbf{X}, t)}{\partial t}\right|_{\mathbf{X}} 
  
and can be interpreted as the particle velocity in the referential domain, since it measures
the time variation of the referential coordinate :math:`\boldsymbol{\chi}` holding the material particle :math:`\mathbf{X}` fixed. The
relation between velocities :math:`\mathbf{v}`, :math:`\hat{\mathbf{v}}`, and :math:`\mathbf{w}` can be obtained by differentiating :math:`\boldsymbol{\varphi}=\boldsymbol{\Phi} \circ \boldsymbol{\Psi}^{-1}`,

.. math::
  \begin{aligned}
  x_{1}={\varphi}_{1}(X_{1},X_{2},t)=\phi_{1}(\chi_{1}(X_{1},X_{2},t),\chi_{2}(X_{1},X_{2},t),t)
  =\phi_{1}(\psi^{-1}_{1}(X_{1},X_{2},t),\psi^{-1}_{2}(X_{1},X_{2},t),t)\\
  x_{2}={\varphi}_{2}(X_{1},X_{2},t)=\phi_{2}(\chi_{1}(X_{1},X_{2},t),\chi_{2}(X_{1},X_{2},t),t)
  =\phi_{2}(\psi^{-1}_{1}(X_{1},X_{2},t),\psi^{-1}_{2}(X_{1},X_{2},t),t)\\
  \end{aligned}
  
-
  
.. math::
  \begin{aligned}
  \cfrac{\partial {\varphi}_{1}(X_{1},X_{2},t)}{\partial X_{1}} =
  \cfrac{\partial \phi_{1}(\chi_{1}(X_{1},X_{2},t),\chi_{2}(X_{1},X_{2},t),t)}{\partial \chi_{1}}\cdot \frac{\partial \psi^{-1}_{1}(X_{1},X_{2},t)}{\partial X_{1}}\\ 
  +\cfrac{\partial \phi_{1}(\chi_{1}(X_{1},X_{2},t),\chi_{2}(X_{1},X_{2},t),t)}{\partial \chi_{2}}\cdot \frac{\partial \psi^{-1}_{2}(X_{1},X_{2},t)}{\partial X_{1}} 
  \end{aligned}
  
 
-

.. math::  
  \begin{aligned}
  \cfrac{\partial {\varphi}_{i}(X_{1},X_{2},t)}{\partial X_{j}} =
  \cfrac{\partial \phi_{i}(\chi_{1}(X_{1},X_{2},t),\chi_{2}(X_{1},X_{2},t),t)}{\partial \chi_{k}}\cdot \frac{\partial \psi^{-1}_{k}(X_{1},X_{2},t)}{\partial X_{j}}\\ 
  \end{aligned}  
  
-

.. math:: 
  \begin{bmatrix}
  \cfrac{\partial {\varphi}_{i}(X_{1},X_{2},t)}{\partial X_{j}}
  \end{bmatrix}=\cfrac{\partial {\varphi}(X_{1},X_{2},t)}{\partial (X_{1},X_{2},t)}=\frac{\partial \boldsymbol{\varphi}(\mathbf{X}, t)}{\partial(\mathbf{X}, t)}
  
 
-

.. math::  
  \cfrac{\partial {\varphi}(X_{1},X_{2},t)}{\partial (X_{1},X_{2},t)}=\begin{bmatrix}
  \cfrac{\partial {\varphi}_{1}(X_{1},X_{2},t)}{\partial X_{1}} & \cfrac{\partial {\varphi}_{1}(X_{1},X_{2},t)}{\partial X_{2}} &\cfrac{\partial {\varphi}_{1}(X_{1},X_{2},t)}{\partial t} \\
  \cfrac{\partial {\varphi}_{2}(X_{1},X_{2},t)}{\partial X_{1}} & \cfrac{\partial {\varphi}_{2}(X_{1},X_{2},t)}{\partial X_{2}} &\cfrac{\partial {\varphi}_{2}(X_{1},X_{2},t)}{\partial t} \\
  \cfrac{\partial t}{\partial X_{1}} & \cfrac{\partial t}{\partial X_{2}} &\cfrac{\partial t}{\partial t} \\
  \end{bmatrix}
  
-

.. math::   
  \begin{bmatrix}
  \cfrac{\partial \phi_{i}(\chi_{1}(X_{1},X_{2},t),\chi_{2}(X_{1},X_{2},t),t)}{\partial \chi_{k}}
  \end{bmatrix}=\cfrac{\partial {\phi}(\chi_{1},\chi_{2},t)}{\partial (\chi_{1},\chi_{2},t)}=
  \frac{\partial \boldsymbol{\Phi}(\boldsymbol{\chi}, t)}{\partial(\boldsymbol{\chi}, t)}  
  
-

.. math::  
  \cfrac{\partial {\phi}(\chi_{1},\chi_{2},t)}{\partial (\chi_{1},\chi_{2},t)}=\begin{bmatrix}
  \cfrac{\partial {\phi}_{1}(\chi_{1},\chi_{2},t)}{\partial \chi_{1}} & \cfrac{\partial {\phi}_{1}(\chi_{1},\chi_{2},t)}{\partial \chi_{2}} &\cfrac{\partial {\phi}_{1}(\chi_{1},\chi_{2},t)}{\partial t} \\
  \cfrac{\partial {\phi}_{2}(\chi_{1},\chi_{2},t)}{\partial \chi_{1}} & \cfrac{\partial {\phi}_{2}(\chi_{1},\chi_{2},t)}{\partial \chi_{2}} &\cfrac{\partial {\phi}_{2}(\chi_{1},\chi_{2},t)}{\partial t} \\
  \cfrac{\partial t}{\partial \chi_{1}} & \cfrac{\partial t}{\partial \chi_{2}} &\cfrac{\partial t}{\partial t} \\
  \end{bmatrix}  

-

.. math:: 
  \begin{bmatrix}
  \cfrac{\partial \psi^{-1}_{k}(X_{1},X_{2},t)}{\partial X_{j}}
  \end{bmatrix}=\cfrac{\partial \psi^{-1}(X_{1},X_{2},t)}{\partial (X_{1},X_{2},t)}=
  \cfrac{\partial \boldsymbol{\Psi}^{-1}(\boldsymbol{X}, t)}{\partial(\boldsymbol{X}, t)}

-

.. math:: 
  \cfrac{\partial {\psi}^{-1}(X_{1},X_{2},t)}{\partial (X_{1},X_{2},t)}=\begin{bmatrix}
  \cfrac{\partial {\psi}^{-1}_{1}(X_{1},X_{2},t)}{\partial X_{1}} & \cfrac{\partial {\psi}^{-1}_{1}(X_{1},X_{2},t)}{\partial X_{2}} &\cfrac{\partial {\psi}^{-1}_{1}(X_{1},X_{2},t)}{\partial t} \\
  \cfrac{\partial {\psi}^{-1}_{2}(X_{1},X_{2},t)}{\partial X_{1}} & \cfrac{\partial {\psi}^{-1}_{2}(X_{1},X_{2},t)}{\partial X_{2}} &\cfrac{\partial {\psi}^{-1}_{2}(X_{1},X_{2},t)}{\partial t} \\
  \cfrac{\partial t}{\partial X_{1}} & \cfrac{\partial t}{\partial X_{2}} &\cfrac{\partial t}{\partial t} \\
  \end{bmatrix}
  
-

.. math:: 
  \cfrac{\partial {\varphi}(X_{1},X_{2},t)}{\partial (X_{1},X_{2},t)}=
  \cfrac{\partial {\phi}(\chi_{1},\chi_{2},t)}{\partial (\chi_{1},\chi_{2},t)} \cdot 
  \cfrac{\partial \psi^{-1}(X_{1},X_{2},t)}{\partial (X_{1},X_{2},t)}  
  

-

.. math::
  \begin{aligned}
  \cfrac{\partial \boldsymbol{\varphi}}{\partial(\mathbf{X}, t)}(\mathbf{X}, t) & =\frac{\partial \boldsymbol{\Phi}}{\partial(\boldsymbol{\chi}, t)}\left(\boldsymbol{\Psi}^{-1}(\mathbf{X}, t)\right) \cfrac{\partial \boldsymbol{\Psi}^{-1}}{\partial(\mathbf{X}, t)}(\mathbf{X}, t) \\
  & =\cfrac{\partial \boldsymbol{\Phi}}{\partial(\boldsymbol{\chi}, t)}(\boldsymbol{\chi}, t) \quad \cfrac{\partial \boldsymbol{\Psi}^{-1}}{\partial(\mathbf{X}, t)}(\mathbf{X}, t)
  \end{aligned}
  
  
or, in matrix format:

.. math::  
  \begin{bmatrix}
  \cfrac{\partial \mathbf{x}}{\partial \mathbf{X}} & \mathbf{v} \\
  \mathbf{0}^T & 1
  \end{bmatrix}=
  \begin{bmatrix}
  \cfrac{\partial \mathbf{x}}{\partial \boldsymbol{\chi}} & \hat{\mathbf{v}} \\
  \mathbf{0}^T & 1
  \end{bmatrix}\cdot 
  \begin{bmatrix}
  \cfrac{\partial \boldsymbol{\chi}}{\partial \mathbf{X}} & \mathbf{w} \\
  \mathbf{0}^T & 1
  \end{bmatrix}  
  
which yields, after block multiplication,

.. math::
  \mathbf{v}=\hat{\mathbf{v}}+\cfrac{\partial \mathbf{x}}{\partial \boldsymbol\chi} \cdot \mathbf{w}

This equation may be rewritten as

.. math::
  \mathbf{c}:=\mathbf{v}-\hat{\mathbf{v}}=\cfrac{\partial \mathbf{x}}{\partial \boldsymbol\chi} \cdot \mathbf{w}
  
-

.. math::
  \mathbf{v}(\mathbf{X}, t)=\left.\frac{\partial \mathbf{x}(\mathbf{X}, t)}{\partial t}\right|_{\mathbf{X}} 

-
  
.. math::
  \hat{\mathbf{v}}(\boldsymbol{\chi}, t)=\left.\frac{\partial \mathbf{x}(\boldsymbol{\chi}, t)}{\partial t}\right|_{\boldsymbol{\chi}}

-

.. math::    
  {\mathbf{w}}(\mathbf{X}, t)=\left.\frac{\partial \boldsymbol{\chi}(\mathbf{X}, t)}{\partial t}\right|_{\mathbf{X}} 

The convective velocity :math:`\mathbf{c}`, should not be confused with :math:`\mathbf{w}`.
As stated before, :math:`\mathbf{w}` is the particle velocity as seen from the referential domain :math:`R\boldsymbol{\chi}`, whereas :math:`\mathbf{c}` is
the particle velocity relative to the mesh as seen from the spatial domain :math:`R\mathbf{x}` (both :math:`\mathbf{v}` and :math:`\hat{\mathbf{v}}` are
variations of coordinate :math:`\mathbf{x}`). 
  
Material, spatial, and referential time derivatives
---------------------------------------------------------------------------

In order to relate the time derivative in the material, spatial, and referential domains, let
a scalar physical quantity be described by :math:`f(\mathbf{x}, t)`, :math:`f^{*}(\boldsymbol{\chi}, t)`, and :math:`f^{* *}(\mathbf{X}, t)` in the spatial,
referential, and material domains respectively. Stars are employed to emphasize that the
functional forms are, in general, different.

Since the particle motion :math:`\boldsymbol{\varphi}` is a mapping, the spatial description :math:`f(\mathbf{x}, t)`, and the material
description :math:`f^{* *}(\mathbf{X}, t)` of the physical quantity can be related as

.. math::
  f^{* *}(\mathbf{X}, t)=f(\boldsymbol{\varphi}(\mathbf{X}, t),t) \quad or \quad f^{* *}=f \circ \boldsymbol{\varphi}
  
-
  
.. math::  
  f^{* *}(X_{1},X_{2}, t)=f({\varphi}_{1}(X_{1},X_{2}, t),{\varphi}_{2}(X_{1},X_{2}, t), t)  
  
-
  
.. math::  
  \begin{align}
  \cfrac{\partial f^{* *}(X_{1},X_{2}, t)}{\partial X_{1}}=
  \cfrac{\partial f({\varphi}_{1}(X_{1},X_{2}, t),{\varphi}_{2}(X_{1},X_{2}, t), t)}{\partial x_{1}}\cdot  
  \cfrac{\partial{\varphi}_{1}(X_{1},X_{2}, t)}{\partial X_{1}}\\
  +\cfrac{\partial f({\varphi}_{1}(X_{1},X_{2}, t),{\varphi}_{2}(X_{1},X_{2}, t), t)}{\partial x_{2}}\cdot  
  \cfrac{\partial{\varphi}_{2}(X_{1},X_{2}, t)}{\partial X_{1}}\\
  \end{align}
 
-
  
.. math::  
  \begin{align}
  \cfrac{\partial f^{* *}(X_{1},X_{2}, t)}{\partial X_{2}}=
  \cfrac{\partial f({\varphi}_{1}(X_{1},X_{2}, t),{\varphi}_{2}(X_{1},X_{2}, t), t)}{\partial x_{1}}\cdot  
  \cfrac{\partial{\varphi}_{1}(X_{1},X_{2}, t)}{\partial X_{2}}\\
  +\cfrac{\partial f({\varphi}_{1}(X_{1},X_{2}, t),{\varphi}_{2}(X_{1},X_{2}, t), t)}{\partial x_{2}}\cdot  
  \cfrac{\partial{\varphi}_{2}(X_{1},X_{2}, t)}{\partial X_{2}}\\
  \end{align} 

-
  
.. math:: 
  \begin{align}
  \cfrac{\partial f^{* *}(X_{1},X_{2}, t)}{\partial X_{j}}=
  \cfrac{\partial f({\varphi}_{1}(X_{1},X_{2}, t),{\varphi}_{2}(X_{1},X_{2}, t), t)}{\partial x_{k}}\cdot  
  \cfrac{\partial{\varphi}_{k}(X_{1},X_{2}, t)}{\partial X_{j}}\\
  \end{align}  
  
The gradient of this expression can be easily computed as

.. math::
  \frac{\partial f^{* *}}{\partial(\mathbf{X}, t)}(\mathbf{X}, t)=\frac{\partial f}{\partial(\mathbf{x}, t)}(\mathbf{x}, t) \quad \frac{\partial \boldsymbol{\varphi}}{\partial(\mathbf{X}, t)}(\mathbf{X}, t)
  
which is amenable to the matrix form 

.. math::
  \cfrac{\partial f^{* *}(X_{1},X_{2}, t)}{\partial (X_{1},X_{2}, t)}=
  \begin{bmatrix}
    \cfrac{\partial f^{* *}(X_{1},X_{2}, t)}{\partial X_{1}}&
    \cfrac{\partial f^{* *}(X_{1},X_{2}, t)}{\partial X_{2}} &
    \cfrac{\partial f^{* *}(X_{1},X_{2}, t)}{\partial t}
  \end{bmatrix}
  
-
  
.. math::
  \cfrac{\partial f(x_{1},x_{2}, t)}{\partial (x_{1},x_{2}, t)}=
  \begin{bmatrix}
    \cfrac{\partial f(x_{1},x_{2}, t)}{\partial x_{1}}&
    \cfrac{\partial f(x_{1},x_{2}, t)}{\partial x_{2}} &
    \cfrac{\partial f(x_{1},x_{2}, t)}{\partial t}
  \end{bmatrix}
  
-
  
.. math::
  \cfrac{\partial f(x_{1},x_{2}, t)}{\partial (x_{1},x_{2}, t)}=
  \begin{bmatrix}
    \cfrac{\partial f(x_{1},x_{2}, t)}{\partial x_{1}}&
    \cfrac{\partial f(x_{1},x_{2}, t)}{\partial x_{2}} &
    \cfrac{\partial f(x_{1},x_{2}, t)}{\partial t}
  \end{bmatrix}
  
-
  
.. math::
  \cfrac{\partial (x_{1},x_{2}, t)}{\partial (X_{1},X_{2}, t)}=
  \begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial X_{1}}&
  \cfrac{\partial x_{1}}{\partial X_{2}} &
  \cfrac{\partial x_{1}}{\partial t}\\
  \cfrac{\partial x_{2}}{\partial X_{1}}&
  \cfrac{\partial x_{2}}{\partial X_{2}} &
  \cfrac{\partial x_{2}}{\partial t}\\
  \cfrac{\partial t}{\partial X_{1}}&
  \cfrac{\partial t}{\partial X_{2}} &
  \cfrac{\partial t}{\partial t}\\
  \end{bmatrix}
  
which is amenable to the matrix form

.. math::
  \begin{bmatrix}
  \cfrac{\partial f^{**}}{\partial \mathbf{X}}&
  \cfrac{\partial f^{**}}{\partial t}
  \end{bmatrix}=
  \begin{bmatrix}
  \cfrac{\partial f}{\partial \mathbf{x}}&
  \cfrac{\partial f}{\partial t}
  \end{bmatrix}
  \begin{bmatrix}
  \cfrac{\partial \mathbf{x}}{\partial \mathbf{X}}&
  \mathbf{v}\\
  \mathbf{0}^{T}&1
  \end{bmatrix}
  
which renders, after block multiplication, a first expression, which is obvious, that is,
:math:`({\partial f^{**}}/{\partial \mathbf{X}})=({\partial f}/{\partial \mathbf{x}})({\partial \mathbf{x}}/{\partial \mathbf{X}})`); however, the second one is more interesting:  

.. math::
  \cfrac{\partial f^{**}(X_{1},X_{2},t)}{\partial t}=\cfrac{\partial f(x_{1},x_{2},t)}{\partial t}+\cfrac{\partial f(x_{1},x_{2},t)}{\partial \mathbf{x}}\cdot \mathbf{v}(X_{1},X_{2},t)
  
-
  
.. math::
  \cfrac{\partial f^{**}}{\partial t}=\cfrac{\partial f}{\partial t}+\cfrac{\partial f}{\partial \mathbf{x}}\cdot \mathbf{v}
  
Note that this is the well-known equation that relates the material and the spatial time
derivatives. Dropping the stars to ease the notation, this relation is finally cast as

.. math::
  \cfrac{\partial f}{\partial t}\Bigg|_{\mathbf{X}}=\cfrac{\partial f}{\partial t}\Bigg|_{\mathbf{x}}+ \mathbf{v}\cdot \nabla f \quad or \quad 
  \cfrac{\mathrm{d} f}{\mathrm{d} t}=\cfrac{\partial f}{\partial t}+ \mathbf{v}\cdot \nabla f
  
which can be interpreted in the usual way: the variation of a physical quantity for a given
particle :math:`\mathbf{X}` is the local variation plus a convective term taking into account the relative motion
between the material and spatial (laboratory) systems. Moreover, in order not to overload the
rest of the text with notation, except for the specific sections, the material time derivative is
denoted as  

.. math::
  \frac{\mathrm{d} \cdot }{\mathrm{d} t} :=\cfrac{\partial \cdot }{\partial t}\bigg|_{\mathbf{X}}
  
and the spatial time derivative as

.. math::
  \cfrac{\partial  \ \cdot }{\partial t} :=\cfrac{\partial  \ \cdot }{\partial t}\bigg|_{\mathbf{x}}
  
The relation between material and spatial time derivatives is now extended to include the
referential time derivative. With the help of mapping :math:`\mathbf{\Psi}`, the transformation from the referential
description :math:`f^{*}(\boldsymbol{\chi},t)` of the scalar physical quantity to the material description :math:`f^{**}(\mathbf{X},t)` can
be written as 
 
.. math::
  f^{**}=f^{*}\circ \mathbf{\Psi}^{-1}
  
and its gradient can be easily computed as

.. math::
  \cfrac{\partial f^{**}}{\partial (\mathbf X,t)}(\mathbf X,t)=\cfrac{\partial f^{*}}{\partial (\boldsymbol \chi,t)}(\boldsymbol \chi,t)\quad \cfrac{\partial\mathbf{\Psi}^{-1}}{\partial (\mathbf X,t)}(\mathbf X,t)
  
or, in matrix form

.. math::
  \begin{bmatrix}
  \cfrac{\partial f^{**}}{\partial \mathbf{X}}&
  \cfrac{\partial f^{**}}{\partial t}
  \end{bmatrix}=
  \begin{bmatrix}
  \cfrac{\partial f^{*}}{\partial \boldsymbol{\chi}}&
  \cfrac{\partial f^{*}}{\partial t}
  \end{bmatrix}
  \begin{bmatrix}
  \cfrac{\partial \boldsymbol{\chi}}{\partial \mathbf{X}}&
  \mathbf{w}\\
  \mathbf{0}^{T}&1
  \end{bmatrix}
  
which renders, after block multiplication,

.. math::
  \cfrac{\partial f^{**}}{\partial t}=\cfrac{\partial f^{*}}{\partial t}+\cfrac{\partial f^{*}}{\partial \boldsymbol \chi}\cdot \mathbf{w}  
  
Note that this equation relates the material and the referential time derivatives. However,
it also requires the evaluation of the gradient of the considered quantity in the referential
domain. This can be done, but in computational mechanics it is usually easier to work in the
spatial (or material) domain. Moreover, in fluids, constitutive relations are naturally expressed
in the spatial configuration and the Cauchy stress tensor, which will be introduced next, is
the natural measure for stresses. Thus, using the definition of w , the
previous equation may be rearranged into  

.. math::
  \cfrac{\partial f^{**}}{\partial t}=\cfrac{\partial f^{*}}{\partial t}+\cfrac{\partial f}{\partial \mathbf{x}}\cdot \mathbf{c}
  
The fundamental ALE relation between material time derivatives, referential time derivatives,
and spatial gradient is finally cast as (stars dropped)

.. math::
  \cfrac{\partial f}{\partial t}\bigg|_{\mathbf{X}}=\cfrac{\partial f}{\partial t}\bigg|_{\boldsymbol{\chi }}+\cfrac{\partial f}{\partial \mathbf{x}}\cdot \mathbf{c}
  =\cfrac{\partial f}{\partial t}\bigg|_{\boldsymbol{\chi }}+ \mathbf{c}\cdot \nabla f
  
and shows that the time derivative of the physical quantity :math:`f` for a given particle :math:`\mathbf{X}`, that is,
its material derivative, is its local derivative (with the reference coordinate :math:`\boldsymbol{\chi}` held fixed) plus
a convective term taking into account the relative velocity :math:`\mathbf{c}` between the material and the
reference system. This equation is equivalent to equation (19) but in the ALE formulation,
that is, when :math:`(\boldsymbol{\chi},t)` is the reference. 

Time derivative of integrals over moving volumes
-------------------------------------------------- 

Consider thus a material volume Vt bounded by a smooth closed surface :math:`S_{t}` whose points
at time t move with the material velocity :math:`\mathbf{v}=\mathbf{v}(\mathbf{x},t)` where :math:`\mathbf{x}\in S_{t}`. A material volume is a
volume that permanently contains the same particles of the continuum under consideration.
The material time derivative of the integral of a scalar function :math:`f(\mathbf{x}, t)` (note that :math:`f` is defined
in the spatial domain) over the time-varying material volume :math:`V_{t}` is given by the following
well-known expression, often referred to as Reynolds transport theorem (see, for instance,
Belytschko et al, 2000 for a detailed proof):

.. math::
  \frac{\mathrm{d}}{\mathrm{d} t} \int_{V_{t}} f(\mathbf{x}, t) \mathrm{d} V=\int_{V_{c} \equiv V_{t}} \frac{\partial f(\mathbf{x}, t)}{\partial t} \mathrm{~d} V+\int_{S_{c} \equiv S_{t}} f(\mathbf{x}, t) \mathbf{v} \cdot \mathbf{n} \mathrm{d} S
  
which holds for smooth functions :math:`f(\mathbf{x}, t)`. The volume integral in the right-hand side is defined
over a control volume :math:`V_{c}` (fixed in space), which coincides with the moving material volume
:math:`V_{t}` at the considered instant, :math:`t`, in time. Similarly, the fixed control surface :math:`S_{c}` coincides at
time :math:`t` with the closed surface :math:`S_{t}` bounding the material volume :math:`V_{t}`. In the surface integral, :math:`\mathbf{n}`
denotes the unit outward normal to the surface :math:`S_{t}` at time :math:`t`, and :math:`\mathbf{v}` is the material velocity of
points of the boundary :math:`S_{t}`. The first term in the right-hand side is the local
time derivative of the volume integral. The boundary integral represents the flux of the scalar
quantity :math:`f` across the fixed boundary of the control volume :math:`V_{c} \equiv V_{t}`.

Noting that

.. math::
  \int_{S_{c} } f(\mathbf{x}, t) \mathbf{v} \cdot \mathbf{n} \mathrm{d} S=
  \int_{V_{c} } \nabla \cdot (f \mathbf{v})  \mathrm{d} V
  
one obtains the alternative form of Reynolds transport theorem:

.. math::
  \frac{\mathrm{d}}{\mathrm{d} t} \int_{V_{t}} f(\mathbf{x}, t) \mathrm{d} V=\int_{V_{c} \equiv V_{t}} \left ( \frac{\partial f(\mathbf{x}, t)}{\partial t} +\nabla \cdot (f \mathbf{v}) \right ) \mathrm{~d} V
  
Similar forms hold for the material derivative of the volume integral of a vector quantity.
Analogous formulae can be developed in the ALE context, that is, with a referential time
derivative. In this case, however, the characterizing velocity is no longer the material velocity
:math:`\mathbf{v}`, but the grid velocity :math:`\mathbf{\hat{v}}`.  

  
Material Motion, Mesh Displacement, Mesh Velocity,and Mesh Acceleration
---------------------------------------------------------------------------
To establish the integral form of the basic conservation laws for mass, momentum, and energy,
we also need to consider the rate of change of integrals of scalar and vector functions over a
moving volume occupied by fluid

In an ALE method, both the motion of the mesh and the material must be described. The
motion of the material is described as before by:

.. math::
  \mathbf{x}=\boldsymbol{\Phi}(\mathbf{X}, t)
  
where :math:`\mathbf{X}` are the material coordinates. The function :math:`\boldsymbol{\Phi}(\mathbf{X}, t)` maps the body from the initial configuration :math:`\Omega_{0}` to the current or spatial configuration :math:`\Omega`. Although it is called the motion
throughout this book, in this chapter we will often call it the *material motion* to distinguish it
from the *mesh motion*. It is identical to the map used to describe the motion of Lagrangian elements.

In the ALE formulation, we consider another reference domain :math:`\hat{\Omega}` as shown in Figure 7.1.
This domain is called the *referential domain* or the *ALE domain*. The initial values of the
position of particles are denoted by :math:`\boldsymbol{\chi}`, so:

.. math::
  \boldsymbol{\chi}=\boldsymbol{\Phi}(\mathbf{X}, 0)
  
The coordinates :math:`\boldsymbol{\chi}` are called the *referential* or *ALE coordinates*. In most cases :math:`\boldsymbol{\Phi}(\mathbf{X}, 0)=\mathbf{X}`, so
:math:`\boldsymbol{\chi}(\mathbf{X}, 0)=\mathbf{X}` . The referential domain :math:`\hat{\Omega}` is used to describe the motion of the mesh independent
of the motion of the material. In the implementation, the domain :math:`\hat{\Omega}` is used to construct the
initial mesh. It remains coincident with the mesh throughout the computation, so it can also be
considered the computational domain.

The motion of the mesh is described by  

.. math::
  \mathbf{x}=\hat{\boldsymbol{\Phi}}(\boldsymbol{\chi}, t)
  
.. figure:: images/ale3.png
   :width: 800
   :align: center
   
   Maps between Lagrangian, Eulerian and ALE domains
   
This map :math:`\hat{\boldsymbol{\phi}}` plays a crucial role in the ALE finite element formulation. Points :math:`\boldsymbol{\chi}` in the ALE
domain, :math:`\hat{\Omega}` , are mapped to points :math:`\boldsymbol{x}` in the spatial domain, :math:`\Omega` , via this map.
As is apparent from Figure 7.1, (7.2.1) and (7.2.3), we can relate the ALE coordinates to the
material coordinates by a composition of functions:


.. math::
  {\boldsymbol{\chi}}={\hat{\Phi}}^{-1}(\mathbf{x},t)={\hat{\Phi}}^{-1}(\Phi(\mathbf{X},t),t)=\Psi(\mathbf{X},t)\quad{\mathrm{or}}\quad\Psi={\hat{\Phi}}^{-1}\circ\Phi 
  
As can be seen from the previous, the relation between the material coordinates and the ALE
coordinates is a function of time.

The material motion can be expressed as a composition of the mesh motion and the :math:`\Psi`
map:  

.. math::
  {\bf x}={\Phi}({\bf X},t)=\hat{\Phi}(\Psi({\bf X},t),t)\quad\mathrm{or}\quad\Phi=\hat{\Phi}\circ\Psi 
  
As will be seen, in the ALE algorithm the mesh motion is prescribed or computed. The material
motion can then be reconstructed through the above composition of functions if the map :math:`\Psi` is
invertible.

We will now define the displacement, velocity and acceleration of the mesh motion, which
will be called the mesh displacement, mesh velocity and mesh acceleration. The mesh
displacement, :math:`\hat{\mathbf{u}}` , is defined by 

.. math::
  {\hat{\mathbf{u}}}(\boldsymbol{\chi},t)=\mathbf{x}-\boldsymbol{\chi}={\hat{\boldsymbol\Phi}} ( \boldsymbol{\chi}, t ) -\boldsymbol{\chi} 
  
Note the similarity of the above definition to the definition of material displacement, which is
:math:`\mathbf{u} = \mathbf{x} – \mathbf{X}`: the material coordinate in the material description has been replaced by the ALE
referential coordinate to obtain the mesh displacement. The mesh velocity is also defined
analogously to the material velocity:

.. math::
  {\hat{\mathbf{v}}}(\boldsymbol{\chi},t)={\frac{\partial{\hat{\boldsymbol\Phi}}(\boldsymbol{\chi},t)}{\partial t}}\equiv{\frac{\partial{\hat{\boldsymbol\Phi}}}{\partial t}}\Bigg |_{\boldsymbol{\chi}}\equiv{\hat{\Phi}}_{,t}\;\left[\boldsymbol{\chi}\right]
  
In the above, the ALE coordinate :math:`\boldsymbol{\chi}` is fixed; in the expression for the material velocity, the
material coordinate :math:`\mathbf{X}` is fixed. Three notations to be used are shown in (7.2.7). When the
independent variables are explicitly given, we simply use the partial derivative with respect to
time to indicate the mesh velocity. If the independent variables are not explicitly given, we
will designate the coordinate which is fixed either by a subscript following a bar or in brackets
following the subscript ':math:`,t`' as shown previously.

The mesh acceleration is given by 

.. math::
  \hat{\mathbf{a}}=\frac{\partial \hat{\mathbf{v}}(\boldsymbol{\chi}, t)}{\partial t}=\frac{\partial^{2} \hat{\mathbf{u}}(\boldsymbol{\chi}, t)}{\partial t^{2}}=\hat{\mathbf{u}}_{,t t[\boldsymbol{\chi}]} 
  
Neither the mesh acceleration nor the mesh velocity have any physical meaning in an ALE
mesh which is not Lagrangian. When the mesh is Lagrangian, they correspond to the material
velocity and acceleration.

Material Time Derivative and Convective Velocity
---------------------------------------------------------------------------
In ALE descriptions, fields are usually expressed as functions of the ALE coordinates :math:`\boldsymbol{\chi}` and
time :math:`t`. The material time derivative (or total derivative) must then be obtained by the chain
rule, similar to the process used in Section 3.2.5 to obtain the material time derivative in an
Eulerian description. Consider a specific function, :math:`f(\boldsymbol{\chi},t)`. Using the chain rule gives

.. math::
  \cfrac{Df}{Dt} \equiv \dot{f} (\boldsymbol{\chi},t)=\cfrac{\partial {f} (\boldsymbol{\chi},t)}{\partial t}
  + \cfrac{\partial {f} (\boldsymbol{\chi},t)}{\partial {\chi_{i}}}\cfrac{\partial {\chi_{i}} }{\partial t}
  =f_{,t}[\boldsymbol{\chi}] + \cfrac{\partial {f} }{\partial {\chi_{i}}}\cfrac{\partial {\chi_{i}} }{\partial t}

-

.. math::  
  \begin{align}
  {x}_{1}={x}_{1}(\chi_{1},\chi_{2},t)=\hat{\phi}_{1}(\chi_{1},\chi_{2},t)\\
  {x}_{2}={x}_{2}(\chi_{1},\chi_{2},t)=\hat{\phi}_{2}(\chi_{1},\chi_{2},t)
  \end{align}
  
-

.. math::    
  \begin{align}
  {v}_{1}={f}_{1}(x_{1},x_{2},t)={f}_{1}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)=\hat{f}_{1}(\chi_{1},\chi_{2},t)\\
  {v}_{2}={f}_{2}(x_{1},x_{2},t)={f}_{2}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)=\hat{f}_{2}(\chi_{1},\chi_{2},t)\\
  \end{align}  

-

.. math:: 
  \begin{align}
  \cfrac{\partial \hat{f}_{1}(\chi_{1},\chi_{2},t)}{\partial t} &=\cfrac{\partial {f}_{1}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial t}\\
  &+\cfrac{\partial {f}_{1}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{1}}\cdot \cfrac{\partial\hat{\phi}_{1}(\chi_{1},\chi_{2},t)}{\partial t}\\
  &+\cfrac{\partial {f}_{1}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{2}}\cdot \cfrac{\partial\hat{\phi}_{2}(\chi_{1},\chi_{2},t)}{\partial t}
  \end{align}
  
-

.. math:: 
  \begin{align}
  \cfrac{\partial \hat{f}_{1}(\chi_{1},\chi_{2},t)}{\partial t} &=\cfrac{\partial {f}_{1}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial t}\\
  &+\cfrac{\partial {f}_{1}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{1}}\cdot {\hat{v}_{1}(\chi_{1},\chi_{2},t)}\\
  &+\cfrac{\partial {f}_{1}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{2}}\cdot {\hat{v}_{2}(\chi_{1},\chi_{2},t)}
  \end{align}
  
-

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} f(\mathbf{x},t)}{\mathrm{D} t}&=\cfrac{\partial f(\mathbf{x},t)}{\partial t}+v_{i}(\mathbf{x},t) \cfrac{\partial f(\mathbf{x},t)}{\partial x_{i}}\\
  &=\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\mathbf{v}(\mathbf{x},t) \cdot \nabla f(\mathbf{x},t)\\
  &=\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\mathbf{v}(\mathbf{x},t) \cdot \operatorname{grad} f(\mathbf{x},t) \\
  \end{array}  
  
-

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} f(\boldsymbol{\chi},t)}{\mathrm{D} t}\equiv \dot{f} (\boldsymbol{\chi},t)=\cfrac{\partial f(\boldsymbol{\chi},t)}{\partial t}+v_{i}(\mathbf{x},t) \cfrac{\partial f(\boldsymbol{\chi},t)}{\partial {\chi}_{i}}\\
  \end{array}   

-

.. math:: 

  \begin{align}
  {h}_{1}={g}_{1}(\chi_{1},\chi_{2},t)={g}_{1}(\psi_{1}(X_{1},X_{2},t),\psi_{2}(X_{1},X_{2},t),t)=\hat{g}_{1}(X_{1},X_{2},t)\\
  {h}_{2}={g}_{2}(\chi_{1},\chi_{2},t)={g}_{2}(\psi_{1}(X_{1},X_{2},t),\psi_{2}(X_{1},X_{2},t),t)=\hat{g}_{2}(X_{1},X_{2},t)
  \end{align}  

-

.. math:: 
  
  \begin{align}
  \chi_{1}=\psi_{1}(X_{1},X_{2},t)\\
  \chi_{2}=\psi_{2}(X_{1},X_{2},t)\\
  \end{align} 
  
-

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} {g}_{1}(\chi_{1},\chi_{2},t)}{\mathrm{D} t}&\equiv \dot{g}_{1}(\chi_{1},\chi_{2},t)\\
  &=\cfrac{\partial {g}_{1}(\chi_{1},\chi_{2},t)}{\partial t}\\
  &+\cfrac{\partial {g}_{1}(\chi_{1},\chi_{2},t)}{\partial {\chi}_{1}}\cdot \cfrac{\partial \psi_{1}(X_{1},X_{2},t)}{\partial t} \\
  &+\cfrac{\partial {g}_{1}(\chi_{1},\chi_{2},t)}{\partial {\chi}_{2}}\cdot \cfrac{\partial \psi_{2}(X_{1},X_{2},t)}{\partial t} \\
  \end{array}    

We now define the referential particle velocity :math:`w_{i}` by

.. math::
  w_{i}=\cfrac{\partial \Psi_{i}(\mathbf{X},t)}{\partial x} =\cfrac{\partial \chi_{i}}{\partial t}\Bigg|_{[\mathbf{X}]} 
  
-

.. math:: 
  \begin{array}{l}
  &w_{1}(X_{1},X_{2},t)=\cfrac{\partial \psi_{1}(X_{1},X_{2},t)}{\partial t}= \cfrac{\partial \chi_{1}(X_{1},X_{2},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  &w_{2}(X_{1},X_{2},t)=\cfrac{\partial \psi_{2}(X_{1},X_{2},t)}{\partial t}= \cfrac{\partial \chi_{2}(X_{1},X_{2},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  \end{array}
  
-

.. math:: 
   w_{i}(\mathbf{X},t)=\cfrac{\partial \psi_{i}(\mathbf{X},t)}{\partial t}= \cfrac{\partial \chi_{i}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  
-

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} {g}_{1}(\chi_{1},\chi_{2},t)}{\mathrm{D} t}\equiv \dot{g}_{1}(\chi_{1},\chi_{2},t)\\
  =\cfrac{\partial {g}_{1}(\chi_{1},\chi_{2},t)}{\partial t}
  +\cfrac{\partial {g}_{1}(\chi_{1},\chi_{2},t)}{\partial {\chi}_{1}}\cdot w_{1}(X_{1},X_{2},t)
  +\cfrac{\partial {g}_{1}(\chi_{1},\chi_{2},t)}{\partial {\chi}_{2}}\cdot w_{2}(X_{1},X_{2},t) \\
  \cfrac{\mathrm{D} {g}_{2}(\chi_{1},\chi_{2},t)}{\mathrm{D} t}\equiv \dot{g}_{2}(\chi_{1},\chi_{2},t)\\
  =\cfrac{\partial {g}_{2}(\chi_{1},\chi_{2},t)}{\partial t}
  +\cfrac{\partial {g}_{2}(\chi_{1},\chi_{2},t)}{\partial {\chi}_{1}}\cdot w_{1}(X_{1},X_{2},t)
  +\cfrac{\partial {g}_{2}(\chi_{1},\chi_{2},t)}{\partial {\chi}_{2}}\cdot w_{2}(X_{1},X_{2},t) \\  
  \end{array}   
  
-

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} {g}(\chi_{1},\chi_{2},t)}{\mathrm{D} t}\equiv \dot{g}(\chi_{1},\chi_{2},t)\\
  =\cfrac{\partial {g}(\chi_{1},\chi_{2},t)}{\partial t}
  +\cfrac{\partial {g}(\chi_{1},\chi_{2},t)}{\partial {\chi}_{1}}\cdot w_{1}(X_{1},X_{2},t)
  +\cfrac{\partial {g}(\chi_{1},\chi_{2},t)}{\partial {\chi}_{2}}\cdot w_{2}(X_{1},X_{2},t) \\
  \end{array}   
  
-

.. math::   
  \begin{array}{l}
  \cfrac{\mathrm{D} {g}(\chi_{1},\chi_{2},t)}{\mathrm{D} t}\equiv \dot{g}(\chi_{1},\chi_{2},t)\\
  =\cfrac{\partial {g}(\chi_{1},\chi_{2},t)}{\partial t}
  +\cfrac{\partial {g}(\chi_{1},\chi_{2},t)}{\partial {\chi}_{j}}\cdot w_{j}(X_{1},X_{2},t)
  \end{array}   
  
-

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} {g}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{g}(\boldsymbol\chi,t)
  =\cfrac{\partial {g}(\boldsymbol\chi,t)}{\partial t}
  +\cfrac{\partial {g}(\boldsymbol\chi,t)}{\partial {\chi}_{j}}\cdot w_{j}(\mathbf{X},t)
  \end{array}  
  
-

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}
  +\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial {\chi}_{j}}\cdot w_{j}(\mathbf{X},t)
  \end{array} 
  
-

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}
  +\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial {\chi}_{j}}\cdot w_{j}(\mathbf{X},t)
  \end{array}   
  
-

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}
  +\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial {\chi}_{j}}\cdot \cfrac{\partial \psi_{j}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}
  +\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial {\chi}_{j}}\cdot \cfrac{\partial \chi_{j}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  \end{array}

-

.. math:: 
   w_{i}(\mathbf{X},t)=\cfrac{\partial \psi_{i}(\mathbf{X},t)}{\partial t}= \cfrac{\partial \chi_{i}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\

  
the following expression for the material time derivative

.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}\Bigg |_{[\boldsymbol\chi]}
  +\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial {\chi}_{j}}\cdot \cfrac{\partial \chi_{i}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}\Bigg |_{[\boldsymbol\chi]}
  +\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial {\chi}_{j}}\cdot w_{j}(\mathbf{X},t)\\
  \end{array}
  
-

.. math::
  \cfrac{Df}{Dt} \equiv \dot{f} (\boldsymbol{\chi},t)
  = f_{,t}[\boldsymbol{\chi}] + \cfrac{\partial {f} }{\partial {\chi_{i}}}w_{i}
  
-

.. math::
  \begin{array}{l}
  \cfrac{Df}{Dt} \equiv \dot{f} (\boldsymbol{\chi},t)
  = f(\boldsymbol\chi,t)_{,t}[\boldsymbol{\chi}] + \cfrac{\partial {f}(\boldsymbol\chi,t) }{\partial {\chi_{i}}}w_{i}(\mathbf{X},t)
  \end{array}  
  
In the formulations to be given later, the ALE field variables are often treated as functions
of the material coordinates :math:`\mathbf{X}` and time. Hence, it is convenient to develop expressions for the
material time derivative in terms of the spatial gradient. 

.. math::
  \mathbf{x}={\Phi}(\mathbf{X}, t)=\hat{\Phi}({\Psi}(\mathbf{X}, t),t)=\hat{\Phi}\circ{\Psi}
  
-

.. math::
  \mathbf{v}(\mathbf{X},t)=\cfrac{\partial \mathbf{\Phi}(\mathbf{X},t)}{\partial x} 
  
-

.. math::
  \begin{align}
  {x}_{1}={x}_{1}(X_{1},X_{2},t)={\phi}_{1}(X_{1},X_{2},t)\\
  {x}_{2}={x}_{2}(X_{1},X_{2},t)={\phi}_{2}(X_{1},X_{2},t)
  \end{align}
  
  
-

.. math::
  \begin{align}
  {v}_{1}={v}_{1}(X_{1},X_{2},t)=\cfrac {\partial{\phi}_{1}(X_{1},X_{2},t)} {\partial t }\\
  {v}_{2}={v}_{2}(X_{1},X_{2},t)=\cfrac {\partial{\phi}_{2}(X_{1},X_{2},t)} {\partial t }
  \end{align}  
  
-  
  
.. math::
  \begin{align}
  {x}_{1}={x}_{1}(X_{1},X_{2},t)={\phi}_{1}(X_{1},X_{2},t)\\
  {x}_{2}={x}_{2}(X_{1},X_{2},t)={\phi}_{2}(X_{1},X_{2},t)
  \end{align}  
  
-  

.. math::
  \begin{align}
  {x}_{1}={x}_{1}(X_{1},X_{2},t)={\phi}_{1}(X_{1},X_{2},t)=\hat{\phi}_{1}(\chi_{1},\chi_{2},t)=\hat{\phi}_{1}(\chi_{1}(X_{1},X_{2},t),\chi_{2}(X_{1},X_{2},t),t)\\
  {x}_{2}={x}_{2}(X_{1},X_{2},t)={\phi}_{2}(X_{1},X_{2},t)=\hat{\phi}_{2}(\chi_{1},\chi_{2},t)=\hat{\phi}_{2}(\chi_{1}(X_{1},X_{2},t),\chi_{2}(X_{1},X_{2},t),t)
  \end{align} 
  
-  

.. math::
  \begin{align}
  \chi_{1}=\psi_{1}(X_{1},X_{2},t)\\
  \chi_{2}=\psi_{2}(X_{1},X_{2},t)
  \end{align}
  
-  

.. math::  
  \begin{align}
  {v}_{1}&={v}_{1}(X_{1},X_{2},t)=\cfrac{ \partial{\phi}_{1}(X_{1},X_{2},t)} {\partial t }\\
  &=\cfrac {\partial \hat {\phi}_{1}({\chi}_{1},{\chi}_{2},t)} {\partial t }+
  \cfrac {\partial \hat {\phi}_{1}({\chi}_{1},{\chi}_{2},t)} {\partial {\chi}_{1} }\cdot 
  \cfrac {\partial \psi_{1}(X_{1},X_{2},t)}{\partial {t} }+
  \cfrac {\partial \hat {\phi}_{1}({\chi}_{1},{\chi}_{2},t)} {\partial {\chi}_{2} }\cdot 
  \cfrac {\partial \psi_{2}(X_{1},X_{2},t)}{\partial {t} }\\
  {v}_{2}&={v}_{2}(X_{1},X_{2},t)=\cfrac {\partial{\phi}_{2}(X_{1},X_{2},t)} {\partial t }\\
  &=\cfrac {\partial \hat {\phi}_{2}({\chi}_{1},{\chi}_{2},t)} {\partial t }+
  \cfrac {\partial \hat {\phi}_{2}({\chi}_{1},{\chi}_{2},t)} {\partial {\chi}_{1} }\cdot 
  \cfrac {\partial \psi_{1}(X_{1},X_{2},t)}{\partial {t} }+
  \cfrac {\partial \hat {\phi}_{2}({\chi}_{1},{\chi}_{2},t)} {\partial {\chi}_{2} }\cdot 
  \cfrac {\partial \psi_{2}(X_{1},X_{2},t)}{\partial {t} }\\
  \end{align}

-  

.. math::  
  \begin{align}
  {v}_{j}&={v}_{j}(X_{1},X_{2},t)=\cfrac{ \partial{\phi}_{j}(X_{1},X_{2},t)} {\partial t }\\
  &=\cfrac {\partial \hat {\phi}_{j}({\chi}_{1},{\chi}_{2},t)} {\partial t }+
  \cfrac {\partial \hat {\phi}_{j}({\chi}_{1},{\chi}_{2},t)} {\partial {\chi}_{k} }\cdot 
  \cfrac {\partial \psi_{k}(X_{1},X_{2},t)}{\partial {t} }\\
  \end{align}
  
-  

.. math:: 
  \begin{align}
  {v}_{j}&={v}_{j}(\mathbf{X},t)=\cfrac{ \partial{\phi}_{j}(\mathbf{X},t)} {\partial t }
  =\cfrac {\partial \hat {\phi}_{j}(\boldsymbol{\chi},t)} {\partial t }+
  \cfrac {\partial \hat {\phi}_{j}(\boldsymbol{\chi},t)} {\partial {\chi}_{k} }\cdot 
  \cfrac {\partial \psi_{k}(\mathbf{X},t)}{\partial {t} }\\
  \end{align}
  
-  

.. math::   
  \begin{align}
  {v}_{j}&={v}_{j}(\mathbf{X},t)=
  \hat{v}_{j}(\boldsymbol{\chi},t)+
  \cfrac {\partial x_{j}(\boldsymbol{\chi},t)} {\partial {\chi}_{k}(\mathbf{X},t) }\cdot 
  \cfrac {\partial \psi_{k}(\mathbf{X},t)}{\partial {t} }\\
  \end{align}

-  

.. math::
  \begin{align}
  {v}_{j}&={v}_{j}(\mathbf{X},t)=
  \hat{v}_{j}(\boldsymbol{\chi},t)+
  \cfrac {\partial x_{j}(\boldsymbol{\chi},t)} {\partial {\chi}_{k}(\mathbf{X},t) }\cdot 
  \cfrac {\partial {\chi}_{k}(\mathbf{X},t)}{\partial {t} }\Bigg |_{[\mathbf{X}]} \\
  \end{align} 
  
-  

.. math::  
  w_{k}(\mathbf{X},t)=\cfrac{\partial \psi_{k}(\mathbf{X},t)}{\partial t}= \cfrac{\partial \chi_{k}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\  

-  

.. math:: 

 \begin{align}
  \cfrac {\partial x_{j}(\boldsymbol{\chi},t)} {\partial {\chi}_{k}(\mathbf{X},t) }\cdot 
  \cfrac {\partial {\chi}_{k}(\mathbf{X},t)}{\partial {t} }\Bigg |_{[\mathbf{X}]}=
  \cfrac {\partial x_{j}(\boldsymbol{\chi},t)} {\partial {\chi}_{k}(\mathbf{X},t) }\cdot 
  w_{k}(\mathbf{X},t)=  \cfrac {\partial x_{j}} {\partial {\chi}_{k} }\cdot w_{k}
  \end{align} 

 
Now we define the convective velocity, :math:`\mathbf{c}`, as the difference between the material and mesh
velocities:

.. math::
  \begin{align}
  c_{1} =v_{1}(\mathbf{X},t)-\hat{v}_{1}(\boldsymbol{\chi},t)\\
  c_{2} =v_{2}(\mathbf{X},t)-\hat{v}_{2}(\boldsymbol{\chi},t)\\
  c_{j} =v_{j}(\mathbf{X},t)-\hat{v}_{j}(\boldsymbol{\chi},t)
  \end{align}

-
  
.. math::
  \begin{align}
  c_{j} =v_{j}(\mathbf{X},t)-\hat{v}_{j}(\boldsymbol{\chi},t)=
  \cfrac {\partial x_{j}(\boldsymbol{\chi},t)} {\partial {\chi}_{k}(\mathbf{X},t) }\cdot 
  \cfrac {\partial {\chi}_{k}(\mathbf{X},t)}{\partial {t} }\Bigg |_{[\mathbf{X}]}=
  \cfrac {\partial x_{j}(\boldsymbol{\chi},t)} {\partial {\chi}_{k}(\mathbf{X},t) }\cdot 
  w_{k}(\mathbf{X},t)=  \cfrac {\partial x_{j}} {\partial {\chi}_{k} }\cdot w_{k}
  \end{align}   
  

This relationship between the convected velocity :math:`\mathbf{c}`, material velocity :math:`\mathbf{v}`, mesh velocity :math:`\mathbf{\hat{v}}` and
the referential velocity :math:`\mathbf{w}` will be used frequently in the ALE formulation.

.. math::
  \begin{align}
  x_{1} = x_{1}(\chi_{1},\chi_{2},t) = \hat{\phi}_{1}(\chi_{1},\chi_{2},t)\\
  x_{2} = x_{2}(\chi_{1},\chi_{2},t) = \hat{\phi}_{2}(\chi_{1},\chi_{2},t)
  \end{align}
  
-
  
.. math::
  \begin{align}
  f(x_{1},x_{2},t)=f(x_{1}(\chi_{1},\chi_{2},t),x_{2}(\chi_{1},\chi_{2},t),t)=
  f(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)=
  \hat{f}(\chi_{1},\chi_{2},t)
  \end{align}
  
-
  
.. math::
  \begin{align}
  \cfrac{\partial \hat{f}(\chi_{1},\chi_{2},t)}{\partial \chi_{1}} =\cfrac{\partial {f}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{1}}\cdot \cfrac{\partial \hat{\phi}_{1}(\chi_{1},\chi_{2},t)}{\partial \chi_{1}}\\
  +\cfrac{\partial {f}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{2}}\cdot \cfrac{\partial \hat{\phi}_{2}(\chi_{1},\chi_{2},t)}{\partial \chi_{1}}\\ 
  \end{align} 
  
-
  
.. math::
  \begin{align}
  \cfrac{\partial \hat{f}(\chi_{1},\chi_{2},t)}{\partial \chi_{2}} =\cfrac{\partial {f}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{1}}\cdot \cfrac{\partial \hat{\phi}_{1}(\chi_{1},\chi_{2},t)}{\partial \chi_{2}}\\
  +\cfrac{\partial {f}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{2}}\cdot \cfrac{\partial \hat{\phi}_{2}(\chi_{1},\chi_{2},t)}{\partial \chi_{2}}\\ 
  \end{align} 

-
  
.. math::
  \begin{align}
  \cfrac{\partial \hat{f}(\chi_{1},\chi_{2},t)}{\partial \chi_{j}} =\cfrac{\partial {f}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{1}}\cdot \cfrac{\partial \hat{\phi}_{1}(\chi_{1},\chi_{2},t)}{\partial \chi_{j}}\\
  +\cfrac{\partial {f}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{2}}\cdot \cfrac{\partial \hat{\phi}_{2}(\chi_{1},\chi_{2},t)}{\partial \chi_{j}}\\ 
  \end{align}
  
-
  
.. math::
  \begin{align}
  \cfrac{\partial \hat{f}(\chi_{1},\chi_{2},t)}{\partial \chi_{j}} =\cfrac{\partial {f}(\hat{\phi}_{1}(\chi_{1},\chi_{2},t),\hat{\phi}_{2}(\chi_{1},\chi_{2},t),t)}{\partial x_{k}}\cdot \cfrac{\partial \hat{\phi}_{k}(\chi_{1},\chi_{2},t)}{\partial \chi_{j}}\\
  \end{align}
  
-
  
.. math::  
  \begin{align}
  \cfrac{\partial \hat{f}(\boldsymbol\chi,t)}{\partial \chi_{j}} =\cfrac{\partial {f}(\mathbf{x},t)}{\partial x_{k}(\boldsymbol\chi,t)}\cdot \cfrac{\partial \hat{\phi}_{k}(\boldsymbol\chi,t)}{\partial \chi_{j}}\\
  \cfrac{\partial \hat{f}(\boldsymbol\chi,t)}{\partial \chi_{j}} =\cfrac{\partial {f}(\mathbf{x},t)}{\partial x_{k}(\boldsymbol\chi,t)}\cdot \cfrac{\partial {x}_{k}(\boldsymbol\chi,t)}{\partial \chi_{j}}\\
  \end{align} 
  
-
  
.. math:: 
  \begin{array}{l}
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}
  +\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial {\chi}_{j}}\cdot \cfrac{\partial \chi_{j}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  \end{array}
  
-
  
.. math:: 
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}
  +\cfrac{\partial \hat{f}(\mathbf{x},t)}{\partial x_{k}(\boldsymbol\chi,t)}\cdot \cfrac{\partial {x}_{k}(\boldsymbol\chi,t)}{\partial \chi_{j}(\mathbf{X},t)}\cdot \cfrac{\partial \chi_{j}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  
-
  
.. math:: 
  w_{k}(\mathbf{X},t)=\cfrac{\partial \psi_{k}(\mathbf{X},t)}{\partial t}= \cfrac{\partial \chi_{k}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\    
  
-
  
.. math:: 
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}
  +\cfrac{\partial \hat{f}(\mathbf{x},t)}{\partial x_{k}(\boldsymbol\chi,t)}\cdot \cfrac{\partial {x}_{k}(\boldsymbol\chi,t)}{\partial \chi_{j}(\mathbf{X},t)}\cdot { w_{j}(\mathbf{X},t)}  
  
-

.. math::
  \begin{align}
  c_{j} =v_{j}(\mathbf{X},t)-\hat{v}_{j}(\boldsymbol{\chi},t)=
  \cfrac {\partial x_{j}(\boldsymbol{\chi},t)} {\partial {\chi}_{k}(\mathbf{X},t) }\cdot 
  \cfrac {\partial {\chi}_{k}(\mathbf{X},t)}{\partial {t} }\Bigg |_{[\mathbf{X}]}=
  \cfrac {\partial x_{j}(\boldsymbol{\chi},t)} {\partial {\chi}_{k}(\mathbf{X},t) }\cdot 
  w_{k}(\mathbf{X},t)=  \cfrac {\partial x_{j}} {\partial {\chi}_{k} }\cdot w_{k}
  \end{align} 
  
-

.. math::  
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  =\cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}
  +\cfrac{\partial \hat{f}(\mathbf{x},t)}{\partial x_{k}(\boldsymbol\chi,t)}\cdot {c}_{k}
  
-

.. math:: 
  \begin{align}
  \hat{f}_{,k} \equiv \cfrac{\partial \hat{f}(\mathbf{x},t)}{\partial x_{k}(\boldsymbol\chi,t)}\\
  {f}_{,t[\boldsymbol{\chi}]} \equiv \cfrac{\partial {f}(\boldsymbol\chi,t)}{\partial t}
  \end{align} 
  
-

.. math::   
  \begin{align}
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  ={f}_{,t[\boldsymbol{\chi}]}
  +\hat{f}_{,k}\cdot {c}_{k}\\
  \end{align}
  
The above gives the material time derivative for the function in terms of the partial time derivative with
the ALE coordinates fixed and a spatial gradient. Note that a comma followed by an index
represents the spatial derivative with respect to an Eulerian coordinate, as in the rest of this
book. In vector notation, the above can be written as  

.. math::   
  \cfrac{\mathrm{D} {f}(\boldsymbol\chi,t)}{\mathrm{D} t}\equiv \dot{f}(\boldsymbol\chi,t)
  ={f}_{,t[\boldsymbol{\chi}]}+\mathbf{c}\cdot \mathrm{grad}{f}
  ={f}_{,t[\boldsymbol{\chi}]}+\mathbf{c}\cdot\nabla {f}\\

-

.. math::    
  \begin{align}
  \nabla & = \frac{\partial}{\partial x} \mathbf{i}+\frac{\partial}{\partial y} \mathbf{j}+\frac{\partial}{\partial z} \mathbf{k}   \\
  \nabla f& = \frac{\partial f}{\partial x} \mathbf{i}+\frac{\partial f}{\partial y} \mathbf{j}+\frac{\partial f}{\partial z} \mathbf{k}   
  \end{align}
  
-

.. math:: 
  \begin{align}
  \mathbf{x}=\begin{bmatrix}
   x\\y\\z
  \end{bmatrix}
  =\begin{bmatrix}
   x_{1}\\x_{2} \\x_{3}
  \end{bmatrix}
  \end{align}
  
-

.. math:: 
  \begin{align}
  \nabla f=\begin{bmatrix}
   \cfrac{\partial f}{\partial x} \\\cfrac{\partial f}{\partial y}\\\cfrac{\partial f}{\partial z}
  \end{bmatrix}
   =\begin{bmatrix}
   \cfrac{\partial f}{\partial x_{1}} \\\cfrac{\partial f}{\partial x_{2}}\\\cfrac{\partial f}{\partial x_{3}}
  \end{bmatrix}
  \end{align}
  
Relationship of ALE Description to Eulerian and Lagrangian Descriptions
---------------------------------------------------------------------------

It is worthwhile at this point to relate Lagrangian and Eulerian descriptions to the ALE description. We begin by letting :math:`\boldsymbol{\chi}=\mathbf{X}`, that is, by letting the ALE coordinates be coincident with the
material coordinates. The mesh motion is then given by

.. math:: 
  \begin{align}
  \boldsymbol{\chi} = \mathbf{X}\\
  \mathbf{x} = \hat{\mathbf{\Phi}}(\boldsymbol{\chi}, t)\longrightarrow \mathbf{x} & = \hat{\mathbf{\Phi}}(\mathbf{X}, t)
  \end{align} 

Since the mesh motion is now identical to the material motion, this indicates that the
mesh is now Lagrangian. This can also be seen by examining the :math:`\Psi` map, which becomes  

.. math:: 
  \boldsymbol{\chi}(\mathbf{X},t)=\mathbf{X}=\Psi(\mathbf{X},t)=I(\mathbf{X})
  
and as indicated above, the :math:`\Psi` map becomes the identity map, i.e. in this case the ALE coordinates are identical to the material coordinates. This does not really say anything new, since
this was our starting point. Nevertheless it is of interest because of the correspondence
which will emerge when we examine the reduction of an ALE formulation to an Eulerian
formulation.

When we let the ALE coordinates correspond to the Eulerian coordinate, i.e. :math:`\boldsymbol{\chi}=\mathbf{x}`, then the
mesh motion is given by

.. math:: 
  \mathbf{x} = \hat{\mathbf{\Phi}}(\boldsymbol{\chi}, t)=\hat{\mathbf{\Phi}}(\mathbf{x}, t)=I(\mathbf{x})\\
  
so the mesh motion is the identity map, i.e. the mesh is fixed in space. The motion for an
Eulerian description is given by

.. math:: 
  \mathbf{x}= {\mathbf{\Phi}}(\mathbf{X}, t)={\mathbf{\Phi}}(\Psi^{-1}(\mathbf{x}, t), t)=I(\mathbf{x})\\
  
-
  
.. math:: 
  {\mathbf{\Phi}}\circ \mathbf{\Psi}^{-1}=I(\mathbf{x})
  
So in the reduction of the ALE description to the Eulerian description,
  
.. math:: 
  {\mathbf{\Phi}}=  \mathbf{\Psi}
  
The reductions are thus duals of each other. In the reduction of ALE to the Lagrangian description, the :math:`\mathbf{\Psi}` map becomes the identity and the mesh motion becomes the material motion. In
the degeneration to the Eulerian description, the mesh motion becomes the identity map, and
the :math:`\mathbf{\Psi}` map becomes the material motion.  

It is also interesting to examine the Eulerian and Lagrangian forms of the material time
derivative which are embedded in the ALE form. Recall the material time derivative can be
expressed for the different descriptions as follows:

.. math:: 
  \begin{array}{l}
  \cfrac{Df(\boldsymbol\chi,t)}{Dt} \equiv \dot{f} (\boldsymbol{\chi},t)
  = f(\boldsymbol\chi,t)_{,t}[\boldsymbol{\chi}] + \cfrac{\partial {f}(\boldsymbol\chi,t) }{\partial {\chi_{i}}}w_{i}(\mathbf{X},t)
  \end{array}  

Lagrangian description :math:`(\mathbf{X},t)`

.. math:: 
  \begin{array}{l}
  \boldsymbol\chi=\mathbf{X}\\
  w_{k}(\mathbf{X},t)=\cfrac{\partial \psi_{k}(\mathbf{X},t)}{\partial t}= \cfrac{\partial \chi_{k}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  w_{k}(\mathbf{X},t)=0\\
  \cfrac{Df(\mathbf{X},t)}{Dt} \equiv \dot{f} (\mathbf{X},t)
  =\cfrac{\partial f(\mathbf{X},t)}{\partial t} + \cfrac{\partial {f}(\mathbf{X},t) }{\partial {X_{i}}}w_{i}(\mathbf{X},t)=\cfrac{\partial f(\mathbf{X},t)}{\partial t}
  \end{array}

Eulerian description :math:`(\mathbf{x},t)`

.. math:: 
  \begin{array}{l}
  \boldsymbol{\chi}=\mathbf{x}\\
  w_{k}(\mathbf{X},t)=\cfrac{\partial \psi_{k}(\mathbf{X},t)}{\partial t}= \cfrac{\partial \chi_{k}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}\\
  w_{k}(\mathbf{X},t)=\cfrac{\partial \psi_{k}(\mathbf{X},t)}{\partial t}= \cfrac{\partial x_{k}(\mathbf{X},t)}{\partial t}\Bigg |_{[\mathbf{X}]}=v_{k}(\mathbf{X},t)\\
  \cfrac{Df(\mathbf{x},t)}{Dt} \equiv \dot{f} (\mathbf{x},t)
  = \cfrac{\partial f(\mathbf{x},t)}{\partial t} + \cfrac{\partial {f}(\mathbf{x},t) }{\partial {x_{i}}}v_{i}(\mathbf{X},t)
  \end{array}  

ALE description :math:`(\boldsymbol{\chi},t)`

.. math:: 
  \begin{array}{l}
  \cfrac{Df(\boldsymbol\chi,t)}{Dt} \equiv \dot{f} (\boldsymbol{\chi},t)
  = f(\boldsymbol\chi,t)_{,t}[\boldsymbol{\chi}] + \cfrac{\partial {f}(\boldsymbol\chi,t) }{\partial {\chi_{i}}}w_{i}(\mathbf{X},t)
  \end{array}  
  
Motion

Lagrangian description


+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+
| Description      |            | General ALE                          |Lagrangian                                |  Eulerian                                |
+==================+============+======================================+==========================================+==========================================+
| Motion           | Material   |:math:`\mathbf{x}=                    | :math:`\mathbf{x}=                       | :math:`\mathbf{x}=                       |
|                  |            |{\mathbf{\Phi}}(\mathbf{X}, t)`       | {\mathbf{\Phi}}(\mathbf{X}, t)`          | {\mathbf{\Phi}}(\mathbf{X}, t)`          |
+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+
|                  | Mesh       |:math:`\mathbf{x}=\hat{\mathbf{\Phi}} | :math:`\mathbf{x}=                       | :math:`\mathbf{x}=                       |
|                  |            |(\boldsymbol{\chi}, t)`               | {\mathbf{\Phi}}(\mathbf{X}, t)`          | {\mathbf{I}}(\mathbf{x})`                |
+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+
|                  |            |                                      | :math:`(\boldsymbol{\chi}=\mathbf{X},    | :math:`(\boldsymbol{\chi}=\mathbf{x},    |
|                  |            |                                      | \hat{\mathbf{\Phi}}=\mathbf{\Phi})`      | \hat{\mathbf{\Phi}}=\mathbf{I})`         |
+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+
| Displacement     | Material   |:math:`\mathbf{u}=                    |:math:`\mathbf{u}=                        |:math:`\mathbf{u}=                        |
|                  |            |{\mathbf{x}}-{\mathbf{X}}`            |{\mathbf{x}}-{\mathbf{X}}`                |{\mathbf{x}}-{\mathbf{X}}`                |
+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+
|                  | Mesh       |:math:`\hat{\mathbf{u}}=              |:math:`\hat{\mathbf{u}}=                  |:math:`\hat{\mathbf{u}}=                  |
|                  |            |{\mathbf{x}}-{\boldsymbol{\chi}}`     |{\mathbf{x}}-{\mathbf{X}}=\mathbf{u}`     |{\mathbf{x}}-{\mathbf{x}}=\mathbf{0}`     |
+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+
| Velocity         | Material   |:math:`\mathbf{v}=\mathbf{u}_         |:math:`\mathbf{v}=\mathbf{u}_             |:math:`\mathbf{v}=\mathbf{u}_             | 
|                  |            |{,t[\mathbf{X}]}`                     |{,t[\mathbf{X}]}`                         |{,t[\mathbf{X}]}`                         |
+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+
|                  | Mesh       |:math:`\hat{\mathbf{v}}=\hat          |:math:`\hat{\mathbf{v}}=\hat{\mathbf{u}}  |:math:`\hat{\mathbf{v}}=\hat{\mathbf{u}}  |
|                  |            |{\mathbf{u}}_{,t[\boldsymbol{\chi}]}` |_{,t[\mathbf{X}]}={\mathbf{v}}`           |_{,t[\mathbf{X}]}={\mathbf{0}}`           |
+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+
| Acceleration     | Material   |:math:`\mathbf{a}=\mathbf{v}_         |:math:`\mathbf{a}=\mathbf{v}_             |:math:`\mathbf{a}=\mathbf{v}_             |
|                  |            |{,t[\mathbf{X}]}`                     |{,t[\mathbf{X}]}`                         |{,t[\mathbf{X}]}`                         |
+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+
|                  | Mesh       |:math:`\hat{\mathbf{a}}=\hat          |:math:`\hat{\mathbf{a}}=\hat{\mathbf{v}}  |:math:`\hat{\mathbf{a}}=\hat{\mathbf{v}}  |
|                  |            |{\mathbf{v}}_{,t[\boldsymbol{\chi}]}` |_{,t[\mathbf{X}]}={\mathbf{a}}`           |_{,t[\mathbf{X}]}={\mathbf{a}}`           | 
+------------------+------------+--------------------------------------+------------------------------------------+------------------------------------------+

  
Displacement, Velocity and Acceleration
----------------------------------------------------

The velocity :math:`\mathbf{v}(\mathbf{X}, t)` is the rate of change of the position vector for a material point, i.e. the
time derivative with :math:`\mathbf{X}` held constant. Time derivatives with :math:`\mathbf{X}` held constant are called **material
time derivatives**, or sometimes **material derivatives**. Material time derivatives are also called
**total derivatives**. The velocity can be written in various forms:

.. math::
  \mathbf{v}(\mathbf{X}, t)=\frac{\partial \mathbf{\Phi}(\mathbf{X}, t)}{\partial t}
  
The acceleration :math:`\mathbf{a}(\mathbf{X}, t)` is the rate of change of velocity of a material point, or in other
words the material time derivative of the velocity, and can be written in the forms

.. math::
  \mathbf{a}(\mathbf{X}, t)=\frac{\partial \mathbf{v}(\mathbf{X}, t)}{\partial t}
  
The above expression is called the material form of the acceleration.
When the velocity is expressed in terms of the spatial coordinates and the time, that is, in an
Eulerian description as in :math:`\mathbf{v}(\mathbf{x}, t)`, the material time derivative is obtained as follows. The
spatial coordinates in :math:`\mathbf{v}(\mathbf{x}, t)` are  expressed as a function of the material coordinates and
, giving :math:`\mathbf{v}(\mathbf{\Phi}(\mathbf{X}, t), t)`. The material time derivative is then obtained by the
chain rule:

.. math::
  \frac{D v_{i}(\mathbf{x}, t)}{D t}=\frac{\partial v_{i}(\mathbf{x}, t)}{\partial t}+\frac{\partial v_{i}(\mathbf{x}, t)}{\partial x_{j}} \frac{\partial \phi_{j}(\mathbf{X}, t)}{\partial t}=\frac{\partial v_{i}}{\partial t}+\frac{\partial v_{i}}{\partial x_{j}} v_{j}
  
-  

.. math::  
  \begin{align}
  {x}_{1}={x}_{1}(X_{1},X_{2},t)={\phi}_{1}(X_{1},X_{2},t)\\
  {x}_{2}={x}_{2}(X_{1},X_{2},t)={\phi}_{2}(X_{1},X_{2},t)
  \end{align}

-  

.. math::  
  \begin{align}
  {v}_{1}={f}_{1}(x_{1},x_{2},t)={f}_{1}({\phi}_{1}(X_{1},X_{2},t),{\phi}_{2}(X_{1},X_{2},t),t)=\hat{f}_{1}(X_{1},X_{2},t)\\
  {v}_{2}={f}_{2}(x_{1},x_{2},t)={f}_{2}({\phi}_{1}(X_{1},X_{2},t),{\phi}_{2}(X_{1},X_{2},t),t)=\hat{f}_{2}(X_{1},X_{2},t)
  \end{align}
  
-  

.. math::  
  \begin{align}
  \cfrac{\partial \hat{f}_{1}(X_{1},X_{2},t)}{\partial t} =\cfrac{\partial {f}_{1}({\phi}_{1}(X_{1},X_{2},t),{\phi}_{2}(X_{1},X_{2},t),t)}{\partial t}\\
  +\cfrac{\partial {f}_{1}({\phi}_{1}(X_{1},X_{2},t),{\phi}_{2}(X_{1},X_{2},t),t)}{\partial x_{1}}\cfrac{\partial{\phi}_{1}(X_{1},X_{2},t)}{\partial t}\\
  +\cfrac{\partial {f}_{1}({\phi}_{1}(X_{1},X_{2},t),{\phi}_{2}(X_{1},X_{2},t),t)}{\partial x_{2}}\cfrac{\partial{\phi}_{2}(X_{1},X_{2},t)}{\partial t}
  \end{align}  
  
-  

.. math::  
  \begin{align}
  \cfrac{\partial \hat{f}_{1}(X_{1},X_{2},t)}{\partial t} &=\cfrac{\partial {f}_{1}(x_{1},x_{2},t)}{\partial t}\\
  &+\cfrac{\partial {f}_{1}(x_{1},x_{2},t)}{\partial x_{1}}\cfrac{\partial{\phi}_{1}(X_{1},X_{2},t)}{\partial t}
  +\cfrac{\partial {f}_{1}(x_{1},x_{2},t)}{\partial x_{2}}\cfrac{\partial{\phi}_{2}(X_{1},X_{2},t)}{\partial t}
  \end{align}
  
-  

.. math:: 
  \begin{align}
  \frac{\mathrm{D} {v}_{1}(x_{1},x_{2},t)}{\mathrm{D} t} &= \cfrac{\partial \hat{v}_{1}(X_{1},X_{2},t)}{\partial t} \\&=\cfrac{\partial {v}_{1}(x_{1},x_{2},t)}{\partial t}\\
  &+\cfrac{\partial {v}_{1}(x_{1},x_{2},t)}{\partial x_{1}}\cfrac{\partial{\phi}_{1}(X_{1},X_{2},t)}{\partial t}
  +\cfrac{\partial {v}_{1}(x_{1},x_{2},t)}{\partial x_{2}}\cfrac{\partial{\phi}_{2}(X_{1},X_{2},t)}{\partial t}
  \end{align}
  
-  

.. math:: 
  \begin{align}
  \frac{\mathrm{D} {v}_{2}(x_{1},x_{2},t)}{\mathrm{D} t} &= \cfrac{\partial \hat{v}_{2}(X_{1},X_{2},t)}{\partial t} \\&=\cfrac{\partial {v}_{2}(x_{1},x_{2},t)}{\partial t}\\
  &+\cfrac{\partial {v}_{2}(x_{1},x_{2},t)}{\partial x_{1}}\cfrac{\partial{\phi}_{1}(X_{1},X_{2},t)}{\partial t}
  +\cfrac{\partial {v}_{2}(x_{1},x_{2},t)}{\partial x_{2}}\cfrac{\partial{\phi}_{2}(X_{1},X_{2},t)}{\partial t}
  \end{align}
  
-  

.. math::   
  \begin{align}
  \frac{\mathrm{D} {v}_{i}(x_{1},x_{2},t)}{\mathrm{D} t} &= \cfrac{\partial \hat{v}_{i}(X_{1},X_{2},t)}{\partial t} \\
  &=\cfrac{\partial {v}_{i}(x_{1},x_{2},t)}{\partial t}\\
  &+\cfrac{\partial {v}_{i}(x_{1},x_{2},t)}{\partial x_{1}}\cfrac{\partial{\phi}_{1}(X_{1},X_{2},t)}{\partial t}
  +\cfrac{\partial {v}_{i}(x_{1},x_{2},t)}{\partial x_{2}}\cfrac{\partial{\phi}_{2}(X_{1},X_{2},t)}{\partial t}
  \end{align}  
  
-  

.. math::    
  \begin{align}
 \frac{\mathrm{D} {v}_{i}(x_{1},x_{2},t)}{\mathrm{D} t} &= \cfrac{\partial \hat{v}_{i}(X_{1},X_{2},t)}{\partial t} \\
  &=\cfrac{\partial {v}_{i}(x_{1},x_{2},t)}{\partial t}+\cfrac{\partial {v}_{i}(x_{1},x_{2},t)}{\partial x_{j}}\cfrac{\partial{\phi}_{j}(X_{1},X_{2},t)}{\partial t}
  \end{align}  
  
-  

.. math::  
  \begin{align}
 \frac{\mathrm{D} {v}_{i}(\mathbf{x},t)}{\mathrm{D} t} &= \cfrac{\partial \hat{v}_{i}(\mathbf{X},t)}{\partial t}
  =\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial t}+\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial x_{j}}\cfrac{\partial{\phi}_{j}(\mathbf{X},t)}{\partial t}
  \end{align}
  
-  

.. math:: 
  \begin{align}
 \frac{\mathrm{D} {v}_{i}(\mathbf{x},t)}{\mathrm{D} t} = \cfrac{\partial \hat{v}_{i}(\mathbf{X},t)}{\partial t}
  &=\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial t}+\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial x_{j}}\hat{v}_{j}(\mathbf{X},t)
  \end{align}  
  
-  

.. math:: 
  \begin{align}
 \frac{\mathrm{D} {v}_{i}(\mathbf{x},t)}{\mathrm{D} t} = \cfrac{\partial \hat{v}_{i}(\mathbf{X},t)}{\partial t}
  &=\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial t}+\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial x_{j}}{v}_{j}(\mathbf{x},t)
  \end{align}  
  
-  

.. math::   
  \begin{align}
 \frac{\mathrm{D} {v}_{i}(\mathbf{x},t)}{\mathrm{D} t} = \cfrac{\partial \hat{v}_{i}(\mathbf{X},t)}{\partial t}
  &=\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial t}+\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial x_{j}}{v}_{j}(\mathbf{x},t)=\cfrac{\partial {v}_{i}}{\partial t}+\cfrac{\partial {v}_{i}}{\partial x_{j}}{v}_{j}
  \end{align}  
  
-  

.. math::   
  \begin{align}
 \frac{\mathrm{D} {v}_{i}(\mathbf{x},t)}{\mathrm{D} t} = \cfrac{\partial \hat{v}_{i}(\mathbf{X},t)}{\partial t}
  &=\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial t}+\cfrac{\partial {v}_{i}(\mathbf{x},t)}{\partial x_{j}}{v}_{j}(\mathbf{x},t)=\cfrac{\partial {v}_{i}}{\partial t}+\cfrac{\partial {v}_{i}}{\partial x_{j}}{v}_{j}
  \end{align}  

-  

.. math:: 
  \begin{align}
  \frac{\mathrm{D} \mathbf{v}(\mathbf{x},t)}{\mathrm{D} t} = \cfrac{\partial \hat{\mathbf{v}}(\mathbf{X},t)}{\partial t}
  &=\cfrac{\partial {\mathbf{v}}}{\partial t}+{\mathbf{v}}\cdot\nabla {\mathbf{v}}=\cfrac{\partial {\mathbf{v}}}{\partial t}+{\mathbf{v}}\cdot {\mathrm{grad}} \ {\mathbf{v}}\\
  \end{align}  
  
where :math:`\nabla {\mathbf{v}}` and :math:`{\mathrm{grad}} \ {\mathbf{v}}` are the left gradients of a vector field as defined in Malvern (1969: 58).
The matrix of the left gradient is given by 

.. math:: 
  \begin{align}
  \mathbf{v}=\begin{bmatrix}
   u\\v\\w
  \end{bmatrix}
  =\begin{bmatrix}
   v_{x}\\v_{y} \\v_{z}
  \end{bmatrix}
  \end{align} 

-
  
.. math:: 
  \nabla \mathbf{v}\equiv \mathrm{grad}\ \mathbf{v}=\begin{bmatrix}
  \cfrac{\partial u}{\partial x} &  \cfrac{\partial u}{\partial y}& \cfrac{\partial u}{\partial z}\\
  \cfrac{\partial v}{\partial x} &  \cfrac{\partial v}{\partial y}& \cfrac{\partial v}{\partial z}\\
  \cfrac{\partial w}{\partial x} &  \cfrac{\partial w}{\partial y}& \cfrac{\partial w}{\partial z}\\
  \end{bmatrix} 

-
  
.. math:: 
  \nabla \mathbf{v}\equiv \mathrm{grad}\ \mathbf{v}=\begin{bmatrix}
  \cfrac{\partial v_{x}}{\partial x} &  \cfrac{\partial v_{x}}{\partial y}& \cfrac{\partial v_{x}}{\partial z}\\
  \cfrac{\partial v_{y}}{\partial x} &  \cfrac{\partial v_{y}}{\partial y}& \cfrac{\partial v_{y}}{\partial z}\\
  \cfrac{\partial v_{z}}{\partial x} &  \cfrac{\partial v_{z}}{\partial y}& \cfrac{\partial v_{z}}{\partial z}\\
  \end{bmatrix}

-
  
.. math:: 
  \nabla \mathbf{v}\equiv \mathrm{grad}\ \mathbf{v}=\begin{bmatrix}
  \cfrac{\partial v_{1}}{\partial x_{1}} &  \cfrac{\partial v_{1}}{\partial x_{2}}& \cfrac{\partial v_{1}}{\partial x_{3}}\\
  \cfrac{\partial v_{2}}{\partial x_{1}} &  \cfrac{\partial v_{2}}{\partial x_{2}}& \cfrac{\partial v_{2}}{\partial x_{3}}\\
  \cfrac{\partial v_{3}}{\partial x_{1}} &  \cfrac{\partial v_{3}}{\partial x_{2}}& \cfrac{\partial v_{3}}{\partial x_{3}}\\
  \end{bmatrix}  
  
The material time derivative of any function of the spatial variables :math:`\mathbf{x}` and time :math:`t` can similarly
be obtained by the chain rule. Thus for a scalar function :math:`f(\mathbf{x},t)` and a tensor function :math:`\sigma_{i j}(\mathbf{x},t)`
the material time derivatives are given by

.. math:: 
  \begin{array}{c}
  \cfrac{\mathrm{D}  f}{\mathrm{D} t}=\cfrac{\partial f}{\partial t}+v_{i} \cfrac{\partial f}{\partial x_{i}}=\cfrac{\partial f}{\partial t}+\mathbf{v} \cdot \nabla f=\cfrac{\partial f}{\partial t}+\mathbf{v} \cdot \operatorname{grad} f \\
  \cfrac{\mathrm{D} \sigma_{i j}}{\mathrm{D} t}=\cfrac{\partial \sigma_{i j}}{\partial t}+v_{k} \cfrac{\partial \sigma_{i j}}{\partial x_{k}}=\cfrac{\partial \boldsymbol{\sigma}}{\partial t}+\mathbf{v} \cdot \nabla \boldsymbol{\sigma}=\cfrac{\partial \boldsymbol{\sigma}}{\partial t}+\mathbf{v} \cdot \operatorname{grad} \boldsymbol{\sigma}
  \end{array}
  
-

.. math::   
  \begin{array}{l}
  \cfrac{\mathrm{D} f(\mathbf{x},t)}{\mathrm{D} t}&=\cfrac{\partial f(\mathbf{x},t)}{\partial t}+v_{i}(\mathbf{x},t) \cfrac{\partial f(\mathbf{x},t)}{\partial x_{i}}\\
  &=\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\mathbf{v}(\mathbf{x},t) \cdot \nabla f(\mathbf{x},t)\\
  &=\cfrac{\partial f(\mathbf{x},t)}{\partial t}+\mathbf{v}(\mathbf{x},t) \cdot \operatorname{grad} f(\mathbf{x},t) \\
  \end{array}  
  
Motion
---------------------------------
The motion of the body is described by

.. math::   
  \mathbf{x}= {\mathbf{\Phi}}(\mathbf{X}, t)\quad or \quad  x_{i}=\phi_{i}(\mathbf{X}, t)
  
where :math:`\mathbf{x}` is the position of the material point :math:`\mathbf{X}` at time :math:`t`. The coordinates :math:`x_{i}` give the spatial
position, and are called spatial or Eulerian coordinates. The function :math:`{\mathbf{\Phi}}(\mathbf{X}, t)` maps the initial
configuration into the current configuration at time t, and is called a mapping or map from the
initial configuration to the current configuration.

Deformation Gradient
---------------------------------
The description of deformation and the measure of strain are essential parts of nonlinear
continuum mechanics. An important variable in the characterization of deformation is the
**deformation gradient**. The deformation gradient is defined by

.. math:: 
  F_{ij}=\cfrac{\partial \phi_{i}(\mathbf{X}, t)}{\partial X_{j}} \equiv \cfrac{\partial x_{i}(\mathbf{X}, t)}{\partial X_{j}} \\
  
-

.. math:: 
  F_{ij}=\cfrac{\partial \phi_{i}}{\partial X_{j}} \equiv \cfrac{\partial x_{i}}{\partial X_{j}} \\
  
- 
 
.. math:: 
  \mathbf{F}=\cfrac{\partial {\mathbf{\Phi}}(\mathbf{X}, t)}{\partial \mathbf{X}} \equiv \cfrac{\partial \mathbf{x}(\mathbf{X}, t)}{\partial\mathbf{X}}\equiv (\nabla_{0}  {\mathbf{\Phi}}) \\  
  
In the terminology of mathematics, the deformation gradient :math:`\mathbf{F}` is the Jacobian matrix of the
motion :math:`{\mathbf{\Phi}}(\mathbf{X}, t)`. Note in the above that the first index of :math:`F_{ij}` refers to the motion, the second to the
partial derivative. The operator :math:`\nabla_{0}` is the *left gradient with respect to the material coordinates*.

If we consider an infinitesimal line segment :math:`d\mathbf{X}` in the reference configuration, then
the corresponding line segment :math:`d\mathbf{x}` in the current configuration is given by  

.. math:: 
  d\mathbf{x}=\mathbf{F}\cdot d\mathbf{X} \quad or \quad dx_{i}=F_{ij}dX_{j}
  
 
In the above expression, the dot could have been omitted between the :math:`\mathbf{F}` and :math:`d\mathbf{X}`, since the
expression is also valid as a matrix expression. We have retained it to conform to our convention
of always explicitly indicating contractions in tensor expressions.
In two dimensions, the deformation gradient in a rectangular coordinate system is given by

.. math:: 
  \mathbf{F}=\begin{bmatrix}
  \cfrac{\partial x_{1}}{\partial X_{1}} & \cfrac{\partial x_{1}}{\partial X_{2}}\\
  \cfrac{\partial x_{2}}{\partial X_{1}} & \cfrac{\partial x_{2}}{\partial X_{2}}\\
  \end{bmatrix}
  =\begin{bmatrix}
  \cfrac{\partial x}{\partial X} & \cfrac{\partial x}{\partial Y}\\
  \cfrac{\partial y}{\partial X} & \cfrac{\partial y}{\partial Y}\\
  \end{bmatrix}  
  
As can be seen in the above, in writing a second-order tensor in matrix form, we use the first
index for the row number, and the second index for the column number. Note that :math:`\mathbf{F}` is the
transpose of the left-gradient.

The determinant of :math:`\mathbf{F}` is denoted by :math:`J` and called the *Jacobian determinant* or the determinant of the deformation gradient  

.. math::
  J=\mathrm{det}(\mathbf{F})
  
The Jacobian determinant can be used to relate integrals in the current and reference configurations by

.. math::
  \int\limits_{\Omega }^{}f(\mathbf{x},t)\mathrm{d}\Omega=\int\limits_{\Omega_{0} }^{}f(\mathbf{\Phi}(\mathbf{x},t),t)J\mathrm{d}\Omega_{0} \quad or \quad 
  \int\limits_{\Omega }^{}f\mathrm{d}\Omega=\int\limits_{\Omega_{0} }^{}fJ\mathrm{d}\Omega_{0}
  
or in two dimensions

.. math::
  \int\limits_{\Omega }^{}f(x,y)\mathrm{d}x\mathrm{d}y=\int\limits_{\Omega_{0} }f(X,Y)J\mathrm{d}X\mathrm{d}Y
  
The material derivative of the Jacobian determinant is given by

.. math:: 
  \cfrac{\mathrm{D}J}{\mathrm{D}t}\equiv j=J\mathrm{div}\ \mathbf{v}=J\nabla \mathbf{v}\equiv J\cfrac{\partial v_{i}}{\partial x_{i}}
  
  
The last mapping involves the relationship between the referential coordinates, :math:`\boldsymbol{\chi}`, and the
material coordinates :math:`\boldsymbol{X}`:  
  
  
.. figure:: images/ale4.png
   :width: 600
   :align: center
   
   One-dimensional example of Lagrangian, Eulerian and ALE mesh and particle motion.
   
.. figure:: images/ale5.png
   :width: 600
   :align: center
   
   Rotating rod example   
   
Mesh Descriptions
-----------------------

Spatial coordinates are denoted by :math:`\mathbf{x}` and are also called Eulerian coordinates. A spatial
coordinate specifies the location of a point in space. Material coordinates, also called
Lagrangian coordinates, are denoted by :math:`\mathbf{X}`. The material coordinate labels a material point:
each material point has a unique material coordinate, which is usually taken to be its spatial
coordinate in the initial configuration of the body, so at :math:`t = 0`, :math:`\mathbf{X}=\mathbf{x}`.

.. figure:: images/ale6.png
   :width: 600
   :align: center
   
   Space–time depiction of one-dimensional Lagrangian and Eulerian elements
   
The motion or deformation of a body is described by a function :math:`\boldsymbol{\phi}(\mathbf{X}, t)`, with the material
coordinates :math:`\mathbf{X}` and the time :math:`t` as the independent variables. This function gives the spatial
positions of the material points as a function of time through

.. math::
  :label: gs1

  \mathbf{x}=\boldsymbol{\phi}(\mathbf{X}, t)
  
This is also called a map between the initial and current configurations. The displacement :math:`\mathbf{u}` of
a material point is the difference between its current position and its original position:  
   
.. math::   
  :label: gs2
  
  \mathbf{u}(\mathbf{X}, t)=\boldsymbol{\phi}(\mathbf{X}, t)-\mathbf{X}   
  
To illustrate these definitions, consider the following motion in one dimension:

.. math::
  :label: gs3
  
  x={\phi}(X, t)=(1-X)t+\frac{1}{2}Xt^{2}+X
  
In these equations, the material and spatial coordinates have been changed to scalars since the
motion is one-dimensional. A motion is shown in the figure above; the motions
of several material points are plotted in space-time to exhibit their trajectories. The velocity of
a material point is the time derivative of the motion with the material coordinate fixed, that is,
the velocity is given by

.. math::
  :label: gs4
  
  v = v(X, t)=\cfrac{\partial {\phi}(X, t)}{\partial t}=(1-X)+2\cfrac{1}{2}Xt=1+X(t-1)
  
The mesh description depends on the choice of independent variables. For purposes of
illustration, let us consider the velocity field. We can describe the velocity field as a function
of the Lagrangian (material) coordinates, as in (:eq:`gs4`), or we can describe the velocity as a
function of the Eulerian (spatial) coordinates:

.. math::
  :label: gs5
  
  \overline{v} =\overline{v}(x,t)= v(X, t)=v({\phi}^{-1}(x,t), t)
  
In these expressions we have placed a bar over the velocity symbol to indicate that the velocity
field, when expressed in terms of the spatial coordinate :math:`x` and the time :math:`t`, will not be the same
function as that given in (:eq:`gs4`). We have also used an inverse map to express the material
coordinates in terms of the spatial coordinates:

.. math::
  :label: gs6
  
  X={\phi}^{-1}(x,t)
  
Such inverse mappings can generally not be expressed in closed form for arbitrary motions,
but they are an important conceptual device. For the simple motion given in (:eq:`gs3`), the inverse
map is given by

.. math::
  \begin{array}{c}
    x = {\phi}(X, t) = (1-X)t+\frac{1}{2}Xt^{2}+X\\
    x = t-Xt+\frac{1}{2}Xt^{2}+X=t+X(-t+\frac{1}{2}t^{2}+1)\\
    x-t=X(-t+\frac{1}{2}t^{2}+1)\\
    X=\cfrac{x-t}{-t+\cfrac{1}{2}t^{2}+1} \\
    X=\cfrac{x-t}{\cfrac{1}{2}t^{2}-t+1} \\
  \end{array}

.. math::
  :label: gs7
  
  \begin{array}{c}
    X=\cfrac{x-t}{\cfrac{1}{2}t^{2}-t+1} \\
  \end{array}


Substituting the (:eq:`gs7`) into (:eq:`gs4`) gives

.. math::
  :label: gs8
  
  \bar{v}(x, t)=1+\frac{(x-t)(t-1)}{\frac{1}{2} t^{2}-t+1}=\frac{1-x+x t-\frac{1}{2} t^{2}}{\frac{1}{2} t^{2}-t+1}

Equations (:eq:`gs4`) and (:eq:`gs8`) give the same physical velocity fields, but express them in terms of
different independent variables. Equation (:eq:`gs4`) is called a Lagrangian (material) description, for
it expresses the dependent variable in terms of the Lagrangian (material) coordinates. Equation
(:eq:`gs8`) is called an Eulerian (spatial) description, for it expresses the dependent variable as a
function of the Eulerian (spatial) coordinates. Mathematically, the velocities in the two
descriptions are different functions. Henceforth in this book, we will seldom use different symbols
for different functions when they pertain to the same field, but keep in mind that if a field variable
is expressed in terms of different independent variables, then the functions must be different. In
this book, a symbol for a dependent variable is associated with the field, not the function.

The differences between Lagrangian and Eulerian meshes are most clearly seen in the
behavior of the nodes. If the mesh is Eulerian, the Eulerian coordinates of nodes are fixed, that
is, the nodes are coincident with spatial points. If the mesh is Lagrangian, the Lagrangian
(material) coordinates of nodes are time invariant, that is, the nodes are coincident with
material points. This is illustrated in Figure 1.1. In the Eulerian mesh, the nodal trajectories
are vertical lines and material points pass across element interfaces. In the Lagrangian mesh,
nodal trajectories are coincident with material point trajectories, and no material passes
between elements. Furthermore, element quadrature points remain coincident with material
points in Lagrangian meshes, whereas in Eulerian meshes the material point at a given
quadrature point changes with time. We will see later that this complicates the treatment of
materials for which the stress is history-dependent.

The comparative advantages of Eulerian and Lagrangian meshes can be seen even in this
simple one-dimensional example. Since the nodes are coincident with material points in the
Lagrangian mesh, boundary nodes remain on the boundary throughout the evolution of the
problem. This simplifies the imposition of boundary conditions in Lagrangian meshes. In
Eulerian meshes, on the other hand, boundary nodes do not remain coincident with the
boundary. Therefore, boundary conditions must be imposed at points which are not nodes, and
this engenders significant complications in multi-dimensional problems. Similarly, if a node
is placed on an interface between two materials, it remains on the interface in a Lagrangian
mesh, but not in an Eulerian mesh.

In Lagrangian meshes, since the material points remain coincident with mesh points,
elements deform with the material. Therefore, elements in a Lagrangian mesh can become
severely distorted. This effect is apparent in a one-dimensional problem only in the element
lengths: in Eulerian meshes, element lengths are constant in time, whereas in Lagrangian
meshes, element lengths change with time. In multi-dimensional problems, these effects are
far more severe, and Lagrangian elements can get very distorted. Since element accuracy
degrades with distortion, the magnitude of deformation that can be simulated with a Lagrangian
mesh is limited. Eulerian elements, on the other hand, are unchanged by the deformation of
the material, so no degradation in accuracy occurs because of material deformation.

To illustrate the differences between Eulerian and Lagrangian mesh descriptions, a twodimensional example will be considered. The spatial coordinates are denoted by :math:`\mathbf{x} = [x, y]^{T}`
and the material coordinates by :math:`\mathbf{X} = [X, Y]^{T}`. The motion is given by

.. math::
  \mathbf{x} = \Phi(\mathbf{X},t)=\begin{bmatrix}
   x\\y
  \end{bmatrix}\\
  
where :math:`\Phi(\mathbf{X},t)` is a vector function, i.e. it gives a vector for every pair of the independent variables. Writing out the above expression gives

.. math::
  x=\phi_{1}(X, Y, t) \quad y=\phi_{2}(X, Y, t)
  
As an example of a motion, consider a pure shear

.. math::
  x=X+tY \quad y=Y
  
In a Lagrangian mesh, the nodes are coincident with material (Lagrangian) points, so for
Lagrangian nodes, :math:`\mathbf{X}_{I} =constant` in time

For an Eulerian mesh, the nodes are coincident with spatial (Eulerian) points, so for Eulerian
nodes, :math:`\mathbf{x}_{I} =constant` in time

Points on the edges of elements behave similarly to the nodes: in two-dimensional
Lagrangian meshes, element edges remain coincident with material lines, whereas in Eulerian
meshes, the element edges remain fixed in space.