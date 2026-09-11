Generalized Curvilinear Coordinate System(3D)
==============================================

- The Geometric Conservation Law - A link between finite-difference and finite-volume methods of flow computation on moving grids
- Development of CFD Algorithms for Transient and Steady Aerodynamics- Fergal Boyle [Thesis] 


Applying the Curvilinear Transformation
-------------------------------------------
The general curvilinear axes for a time-dependent curvilinear coordinate system are:

.. math::  
  \begin{align}
  \xi & = \xi(x,y,z,t)\\
  \eta & = \eta(x,y,z,t)\\
  \zeta & = \zeta (x,y,z,t)\\
  \tau&=\tau(x,y,z,t)= t
  \end{align} 
  
Using the chain rule for a function of multiple variables, the Cartesian derivatives can be written in terms of the curvilinear derivatives as:

.. math::  
  \begin{align}
  \cfrac{\partial }{\partial x}
  & = \xi_{x} \cfrac{\partial}{\partial \xi}
   + \eta_{x} \cfrac{\partial}{\partial \eta}
   + \zeta_{x} \cfrac{\partial}{\partial \zeta}
   + \tau_{x} \cfrac{\partial}{\partial \tau}\\
  \cfrac{\partial }{\partial y}
  & = \xi_{y} \cfrac{\partial}{\partial \xi}
   + \eta_{y} \cfrac{\partial}{\partial \eta}
   + \zeta_{y} \cfrac{\partial}{\partial \zeta}
   + \tau_{y} \cfrac{\partial}{\partial \tau}\\
  \cfrac{\partial }{\partial z}
  & = \xi_{z} \cfrac{\partial}{\partial \xi}
   + \eta_{z} \cfrac{\partial}{\partial \eta}
   + \zeta_{z} \cfrac{\partial}{\partial \zeta}
   + \tau_{z} \cfrac{\partial}{\partial \tau}\\
  \cfrac{\partial }{\partial t}
  & = \xi_{t} \cfrac{\partial}{\partial \xi}
   + \eta_{t} \cfrac{\partial}{\partial \eta}
   + \zeta_{t} \cfrac{\partial}{\partial \zeta}
   + \tau_{t} \cfrac{\partial}{\partial \tau}
  \end{align}
  
where :math:`\xi_{x},\xi_{y},\xi_{z},\xi_{t},\eta_{x},\eta_{y},\eta_{z},\eta_{t},\zeta_{x},\zeta_{y},\zeta_{z},\zeta_{t}`  are the metrics of the
transformation. (noting that :math:`\tau_{x}=\tau_{y}=\tau_{z}=0` and :math:`\tau_{t}=1`). The above equation can be written in the following matrix form:

.. math::
  \left[\begin{array}{c}
  \frac{\partial}{\partial x} \\
  \frac{\partial}{\partial y} \\
  \frac{\partial}{\partial z} \\
  \frac{\partial}{\partial t}
  \end{array}\right]=\left[\begin{array}{llll}
  \xi_{x} & \eta_{x} & \zeta_{x} & 0 \\
  \xi_{y} & \eta_{y} & \zeta_{y} & 0 \\
  \xi_{z} & \eta_{z} & \zeta_{z} & 0 \\
  \xi_{t} & \eta_{t} & \zeta_{t} & 1
  \end{array}\right]\left[\begin{array}{c}
  \frac{\partial}{\partial \xi} \\
  \frac{\partial}{\partial \eta} \\
  \frac{\partial}{\partial \zeta} \\
  \frac{\partial}{\partial \tau}
  \end{array}\right] 
  
It is also possible to expand the curvilinear derivatives in terms of the Cartesian derivatives with the aid of the chain rule:

.. math::
  \begin{aligned}
  \frac{\partial}{\partial \xi} & =x_{\xi} \frac{\partial}{\partial x}+y_{\xi} \frac{\partial}{\partial y}+z_{\xi} \frac{\partial}{\partial z}+t_{\xi} \frac{\partial}{\partial t} \\
  \frac{\partial}{\partial \eta} & =x_{\eta} \frac{\partial}{\partial x}+y_{\eta} \frac{\partial}{\partial y}+z_{\eta} \frac{\partial}{\partial z}+t_{\eta} \frac{\partial}{\partial t} \\
  \frac{\partial}{\partial \zeta} & =x_{\zeta} \frac{\partial}{\partial x}+y_{\zeta} \frac{\partial}{\partial y}+z_{\zeta} \frac{\partial}{\partial z}+t_{\zeta} \frac{\partial}{\partial t} \\
  \frac{\partial}{\partial \tau} & =x_{\tau} \frac{\partial}{\partial x}+y_{\tau} \frac{\partial}{\partial y}+z_{\tau} \frac{\partial}{\partial z}+t_{\tau} \frac{\partial}{\partial t}
  \end{aligned}
  
Note again that :math:`t_{\xi}=t_{\eta}=t_{\zeta}=0` and :math:`t_{\tau}=1`. The equation set above can also be written in matrix form:

.. math::
  \left[\begin{array}{c}
  \frac{\partial}{\partial \xi} \\
  \frac{\partial}{\partial \eta} \\
  \frac{\partial}{\partial \zeta} \\
  \frac{\partial}{\partial \tau}
  \end{array}\right]
  =\left[\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}& {z}_{\xi} & 0 \\
  {x}_{\eta} & {y}_{\eta}& {z}_{\eta} & 0 \\
  {x}_{\zeta} & {y}_{\zeta}& {z}_{\zeta} & 0 \\
  x_{\tau} & y_{\tau} & z_{\tau} & 1
  \end{array}\right]=\left[\begin{array}{c}
  \frac{\partial}{\partial x} \\
  \frac{\partial}{\partial y} \\
  \frac{\partial}{\partial z} \\
  \frac{\partial}{\partial t}
  \end{array}\right]
  
Comparing the above equation sets, it is evident that the following holds true:

.. math::
  \left[\begin{array}{llll}
  \xi_{x} & \eta_{x} & \zeta_{x} & 0 \\
  \xi_{y} & \eta_{y} & \zeta_{y} & 0 \\
  \xi_{z} & \eta_{z} & \zeta_{z} & 0 \\
  \xi_{t} & \eta_{t} & \zeta_{t} & 1
  \end{array}\right]=
  \left[\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}& {z}_{\xi} & 0 \\
  {x}_{\eta} & {y}_{\eta}& {z}_{\eta} & 0 \\
  {x}_{\zeta} & {y}_{\zeta}& {z}_{\zeta} & 0 \\
  x_{\tau} & y_{\tau} & z_{\tau} & 1
  \end{array}\right]^{-1}
  
If the 4x4 matrix on the right-hand side is inverted, it is then possible to solve for the metrics of the transformation.
The inverse of a matrix :math:`\mathbf{D}` can be obtained using the following expression:

.. math::
  \mathbf{D}=\frac{1}{determinant(\mathbf{D})}adjoint(\mathbf{D})
  
Let

.. math::
  A=\left[\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}& {z}_{\xi} & 0 \\
  {x}_{\eta} & {y}_{\eta}& {z}_{\eta} & 0 \\
  {x}_{\zeta} & {y}_{\zeta}& {z}_{\zeta} & 0 \\
  x_{\tau} & y_{\tau} & z_{\tau} & 1
  \end{array}\right]
  
 
The determinant of the matrix :math:`A` is:

.. math::
  \text{det}{A}=
  \left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}& {z}_{\xi} & 0 \\
  {x}_{\eta} & {y}_{\eta}& {z}_{\eta} & 0 \\
  {x}_{\zeta} & {y}_{\zeta}& {z}_{\zeta} & 0 \\
  x_{\tau} & y_{\tau} & z_{\tau} & 1
  \end{array}\right|
  ={x}_{\xi}({y}_{\eta}{z}_{\zeta}-{z}_{\eta}{y}_{\zeta})
  -{y}_{\xi}({x}_{\eta}{z}_{\zeta}-{z}_{\eta}{x}_{\zeta})
  +{z}_{\xi}({x}_{\eta}{y}_{\zeta}-{y}_{\eta}{x}_{\zeta}) 

The Jacobian of the inverse transformation is defined as:

.. math::
  J^{-1}=\cfrac{\partial (x,y,z,t)}{\partial (\xi,\eta,\zeta,\tau)}=
  \left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}& {z}_{\xi} & 0 \\
  {x}_{\eta} & {y}_{\eta}& {z}_{\eta} & 0 \\
  {x}_{\zeta} & {y}_{\zeta}& {z}_{\zeta} & 0 \\
  x_{\tau} & y_{\tau} & z_{\tau} & 1
  \end{array}\right|

-
  
.. math::
  J^{-1}=\text{det}({A})
  
Following the matrix inversion the metrics of the transformation are evaluated as:

.. math::
  A^{-1}=\cfrac{1}{\text{det}{A}}A^{*}

-

.. math::
  \mathbf{A}^{*}=\begin{bmatrix}
  A_{11}&A_{21}  &\cdots   & A_{n1}\\
  A_{12}&A_{22}  &\cdots   & A_{n2}\\
  \vdots& \vdots &  &\vdots \\
  A_{1n}&A_{2n}  &\cdots   & A_{nn}\\
  \end{bmatrix}  

-
  
.. math::
  A_{ij}=(-1)^{i+j}M_{ij}
  
then

.. math::
  A^{-1}=J\begin{bmatrix}
  A_{11}& A_{21} & A_{31} & A_{41}\\
  A_{12}& A_{22} & A_{32} & A_{42}\\
  A_{13}& A_{23} & A_{33} & A_{43}\\
  A_{14}& A_{24} & A_{34} & A_{44}\\
  \end{bmatrix}

-
  
.. math::
  \left[\begin{array}{llll}
  \xi_{x} & \eta_{x} & \zeta_{x} & 0 \\
  \xi_{y} & \eta_{y} & \zeta_{y} & 0 \\
  \xi_{z} & \eta_{z} & \zeta_{z} & 0 \\
  \xi_{t} & \eta_{t} & \zeta_{t} & 1
  \end{array}\right]=
  A^{-1}=J\begin{bmatrix}
  A_{11}& A_{21} & A_{31} & A_{41}\\
  A_{12}& A_{22} & A_{32} & A_{42}\\
  A_{13}& A_{23} & A_{33} & A_{43}\\
  A_{14}& A_{24} & A_{34} & A_{44}\\
  \end{bmatrix}  

specifically:
  
-
  
.. math::
  A_{11}=
  +\left|\begin{array}{llll}
  {y}_{\eta}& {z}_{\eta} & 0 \\
  {y}_{\zeta}& {z}_{\zeta} & 0 \\
   y_{\tau} & z_{\tau} & 1
  \end{array}\right| =+({y}_{\eta}{z}_{\zeta}-{z}_{\eta}{y}_{\zeta})

-
  
.. math::
  A_{12}=
  -\left|\begin{array}{llll}
  {x}_{\eta} & {z}_{\eta} & 0 \\
  {x}_{\zeta} & {z}_{\zeta} & 0 \\
  x_{\tau} &  z_{\tau} & 1
  \end{array}\right| =-({x}_{\eta}{z}_{\zeta}-{z}_{\eta}{x}_{\zeta}) 
  
-
  
.. math::
  A_{13}
  =+  \left|\begin{array}{llll}
  {x}_{\eta} & {y}_{\eta} & 0 \\
  {x}_{\zeta} & {y}_{\zeta} & 0 \\
  x_{\tau} & y_{\tau}  & 1
  \end{array}\right| =+({x}_{\eta}{y}_{\zeta}-{y}_{\eta}{x}_{\zeta}) \\
  
-
  
.. math::
  \begin{align}
  A_{14} & = - \left|\begin{array}{lll}
  {x}_{\eta} & {y}_{\eta}& {z}_{\eta}\\
  {x}_{\zeta} & {y}_{\zeta}& {z}_{\zeta} \\
  x_{\tau} & y_{\tau} & z_{\tau}
  \end{array}\right| = -x_{\tau}\left|\begin{array}{ll}
  {y}_{\eta}& {z}_{\eta}\\
  {y}_{\zeta}& {z}_{\zeta} \\
  \end{array}\right|
  +y_{\tau}\left|\begin{array}{ll}
  {x}_{\eta} & {z}_{\eta}\\
  {x}_{\zeta} & {z}_{\zeta} \\
 \end{array}\right|
 -z_{\tau}\left|\begin{array}{ll}
  {x}_{\eta} & {y}_{\eta}\\
  {x}_{\zeta} & {y}_{\zeta} \\
 \end{array}\right|\\
  &=-x_{\tau}({y}_{\eta}{z}_{\zeta}-{y}_{\zeta}{z}_{\eta})
  -y_{\tau}({z}_{\eta}{x}_{\zeta}-{z}_{\zeta}{x}_{\eta})
  -z_{\tau}({x}_{\eta}{y}_{\zeta}-{x}_{\zeta}{y}_{\eta})
  \end{align}  
  
-
  
.. math::
  A_{21}
  =-\left|\begin{array}{llll}
  {y}_{\xi}& {z}_{\xi} & 0 \\
  {y}_{\zeta}& {z}_{\zeta} & 0 \\
  y_{\tau} & z_{\tau} & 1
  \end{array}\right| =-({y}_{\xi}{z}_{\zeta}-{z}_{\xi}{y}_{\zeta})
  
-
  
.. math::
  A_{22}
  = +\left|\begin{array}{llll}
  {x}_{\xi} & {z}_{\xi} & 0 \\
  {x}_{\zeta} & {z}_{\zeta} & 0 \\
  x_{\tau}  & z_{\tau} & 1
  \end{array}\right| =({x}_{\xi}{z}_{\zeta}-{z}_{\xi}{x}_{\zeta})\\
  
-
  
.. math::  
  A_{23}
  =- \left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi} & 0 \\
  {x}_{\zeta} & {y}_{\zeta} & 0 \\
  x_{\tau} & y_{\tau}  & 1
  \end{array}\right|=- ({x}_{\xi}{y}_{\zeta}-{y}_{\xi}{x}_{\zeta})\\

-
  
.. math::
  \begin{align}
  A_{24} & = + \left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}& {z}_{\xi}\\
  {x}_{\zeta} & {y}_{\zeta}& {z}_{\zeta} \\
  x_{\tau} & y_{\tau} & z_{\tau}
  \end{array}\right|  = x_{\tau}\left|\begin{array}{llll}
  {y}_{\xi}& {z}_{\xi}\\
  {y}_{\zeta}& {z}_{\zeta} \\
  \end{array}\right|
  -y_{\tau}\left|\begin{array}{llll}
  {x}_{\xi} & {z}_{\xi}\\
  {x}_{\zeta} & {z}_{\zeta} \\
  \end{array}\right|
  +z_{\tau}\left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}\\
  {x}_{\zeta} & {y}_{\zeta} \\
  \end{array}\right|\\
  &=x_{\tau}({y}_{\xi}{z}_{\zeta}-{z}_{\xi}{y}_{\zeta})
  -y_{\tau}({x}_{\xi}{z}_{\zeta}-{z}_{\xi}{x}_{\zeta})
  +z_{\tau}({x}_{\xi}{y}_{\zeta}-{y}_{\xi}{x}_{\zeta})
  \end{align}

  
-
  
.. math::  
  A_{31}
  = + \left|\begin{array}{llll}
   {y}_{\xi}& {z}_{\xi} & 0 \\
   {y}_{\eta}& {z}_{\eta} & 0 \\
   y_{\tau} & z_{\tau} & 1
  \end{array}\right|=+({y}_{\xi}{z}_{\eta}-{z}_{\xi}{y}_{\eta})
  
-
  
.. math::
  A_{32}
  =- \left|\begin{array}{llll}
  {x}_{\xi} &  {z}_{\xi} & 0 \\
  {x}_{\eta} &  {z}_{\eta} & 0 \\
  x_{\tau} & z_{\tau} & 1
  \end{array}\right| =-({x}_{\xi}{z}_{\eta}-{z}_{\xi}{x}_{\eta})\\
  
-
  
.. math::
  A_{33}
  = +\left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi} & 0 \\
  {x}_{\eta} & {y}_{\eta} & 0 \\
  x_{\tau} & y_{\tau}  & 1
  \end{array}\right| =+({x}_{\xi}{y}_{\eta}-{y}_{\xi}{x}_{\eta})\\
  
-
  
.. math::
  \begin{align}
  A_{34}
  &=-\left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}& {z}_{\xi}  \\
  {x}_{\eta} & {y}_{\eta}& {z}_{\eta}  \\
  x_{\tau} & y_{\tau} & z_{\tau} 
  \end{array}\right|
  =-x_{\tau}\left|\begin{array}{llll}
   {y}_{\xi}& {z}_{\xi}  \\
   {y}_{\eta}& {z}_{\eta}  \\
  \end{array}\right|
  +y_{\tau}\left|\begin{array}{llll}
  {x}_{\xi} &  {z}_{\xi}  \\
  {x}_{\eta} & {z}_{\eta}  \\
  \end{array}\right|
  -z_{\tau}\left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}  \\
  {x}_{\eta} & {y}_{\eta} \\
  \end{array}\right|\\
  &=-x_{\tau}({y}_{\xi}{z}_{\eta}-{z}_{\xi}{y}_{\eta})
  +y_{\tau}({x}_{\xi}{z}_{\eta}-{z}_{\xi}{x}_{\eta})
  -z_{\tau}({x}_{\xi}{y}_{\eta}-{y}_{\xi}{x}_{\eta})\\
  \end{align}

  
-
  
.. math::
  A_{41}
  =-\left|\begin{array}{llll}
   {y}_{\xi}& {z}_{\xi} & 0 \\
   {y}_{\eta}& {z}_{\eta} & 0 \\
   {y}_{\zeta}& {z}_{\zeta} & 0 \\
  \end{array}\right| =0 
  
-
  
.. math::  
  A_{42}
  =+\left|\begin{array}{llll}
  {x}_{\xi} &  {z}_{\xi} & 0 \\
  {x}_{\eta} & {z}_{\eta} & 0 \\
  {x}_{\zeta} & {z}_{\zeta} & 0 \\
  \end{array}\right|=0
  
-
  
.. math:: 
  A_{43}
  =-\left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi} & 0 \\
  {x}_{\eta} & {y}_{\eta} & 0 \\
  {x}_{\zeta} & {y}_{\zeta} & 0 \\
  \end{array}\right|=0
  
-
  
.. math::
  A_{44}
  =+ \left|\begin{array}{llll}
  {x}_{\xi} & {y}_{\xi}& {z}_{\xi}\\
  {x}_{\eta} & {y}_{\eta}& {z}_{\eta} \\
  {x}_{\zeta} & {y}_{\zeta}& {z}_{\zeta} \\
  \end{array}\right|=J^{-1}  
  
Finally:

.. math::
  \begin{align}
  \xi_{x} = JA_{11} & = +J({y}_{\eta}{z}_{\zeta}-{z}_{\eta}{y}_{\zeta})\\
  \xi_{y} = JA_{12} & = -J({x}_{\eta}{z}_{\zeta}-{z}_{\eta}{x}_{\zeta})\\
  \xi_{z} = JA_{13} & = +J({x}_{\eta}{y}_{\zeta}-{y}_{\eta}{x}_{\zeta})\\
  \xi_{t} = JA_{14} & = +J(-x_{\tau}({y}_{\eta}{z}_{\zeta}-{y}_{\zeta}{z}_{\eta})
  -y_{\tau}({z}_{\eta}{x}_{\zeta}-{z}_{\zeta}{x}_{\eta})
  -z_{\tau}({x}_{\eta}{y}_{\zeta}-{x}_{\zeta}{y}_{\eta}))\\
  & = -x_{\tau}\xi_{x}-y_{\tau}\xi_{y}-z_{\tau}\xi_{z}
  \end{align} 

-

.. math::
  \begin{align}
  \eta_{x} = JA_{21} & = -J({y}_{\xi}{z}_{\zeta}-{z}_{\xi}{y}_{\zeta})\\
  \eta_{y} = JA_{22} & = +J({x}_{\xi}{z}_{\zeta}-{z}_{\xi}{x}_{\zeta})\\
  \eta_{z} = JA_{23} & = -J({x}_{\xi}{y}_{\zeta}-{y}_{\xi}{x}_{\zeta})\\
  \eta_{t} = JA_{24} & = +J(x_{\tau}({y}_{\xi}{z}_{\zeta}-{z}_{\xi}{y}_{\zeta})
  -y_{\tau}({x}_{\xi}{z}_{\zeta}-{z}_{\xi}{x}_{\zeta})
  +z_{\tau}({x}_{\xi}{y}_{\zeta}-{y}_{\xi}{x}_{\zeta}))\\
  & = -x_{\tau}\eta_{x}-y_{\tau}\eta_{y}-z_{\tau}\eta_{z}
  \end{align} 
  
-

.. math::
  \begin{align}
  \zeta_{x} = JA_{31} & = +J({y}_{\xi}{z}_{\eta}-{z}_{\xi}{y}_{\eta})\\
  \zeta_{y} = JA_{32} & = -J({x}_{\xi}{z}_{\eta}-{z}_{\xi}{x}_{\eta})\\
  \zeta_{z} = JA_{33} & = +J({x}_{\xi}{y}_{\eta}-{y}_{\xi}{x}_{\eta})\\
  \zeta_{t} = JA_{34} & = +J(-x_{\tau}({y}_{\xi}{z}_{\eta}-{z}_{\xi}{y}_{\eta})
  +y_{\tau}({x}_{\xi}{z}_{\eta}-{z}_{\xi}{x}_{\eta})
  -z_{\tau}({x}_{\xi}{y}_{\eta}-{y}_{\xi}{x}_{\eta}))\\
  & = -x_{\tau}\zeta_{x}-y_{\tau}\zeta_{y}-z_{\tau}\zeta_{z}
  \end{align}  

Euler Equations in Cartesian Coordinates
------------------------------------------------
The partial differential equation form of the non-dimensional, three-dimensional, Euler equations in Cartesian coordinates in an inertial
reference frame, neglecting volumetric heat addition and body forces, is:

.. math::
  \cfrac{\partial \mathbf{q}}{\partial \text{t}}+
  \cfrac{\partial \mathbf{f}}{\partial \text{x}}+
  \cfrac{\partial \mathbf{g}}{\partial \text{y}}+
  \cfrac{\partial \mathbf{h}}{\partial \text{z}}=0
  
where the vector of sonserved variables, :math:`\mathbf{q}`  , and the vectors of the inviscid flux terms,
:math:`\mathbf{f}`, :math:`\mathbf{g}`, and :math:`\mathbf{h}`, are:


.. math::
  \begin{array}{l}
  \mathbf{q}=\begin{bmatrix}
   \rho\\ \rho u\\ \rho v\\ \rho w \\ \rho E\\
  \end{bmatrix} \quad
  \mathbf{f}=\begin{bmatrix}
   \rho u\\ \rho uu+p\\ \rho vu\\ \rho wu \\ \rho Hu\\
  \end{bmatrix} \quad
  \mathbf{g}=\begin{bmatrix}
   \rho v\\ \rho uv\\ \rho vv+p\\ \rho wv \\ \rho Hv\\
  \end{bmatrix} \quad
  \mathbf{h}=\begin{bmatrix}
   \rho w\\ \rho uw\\ \rho vw\\ \rho ww+p \\ \rho Hw\\
  \end{bmatrix} \quad
  \end{array}
  
where :math:`\rho` is the density, :math:`p` is the static pressure, :math:`u`, :math:`v` and :math:`w` are the Cartesian 
velocity components in the :math:`x`, :math:`y` and :math:`z` directions respectively, :math:`E` is the total energy perunit mass
and :math:`H` is the total enthalpy per unit mass.

.. math::
  \cfrac{\partial\mathbf{q} }{\partial \text{t}}
  = \xi_{t} \cfrac{\partial \mathbf{\hat{q}}}{\partial \xi}
  + \eta_{t} \cfrac{\partial\mathbf{\hat{q}}}{\partial \eta}
  + \zeta_{t} \cfrac{\partial\mathbf{\hat{q}}}{\partial \zeta}
  + \tau_{t} \cfrac{\partial\mathbf{\hat{q}}}{\partial \tau}
  
- 
 
.. math::
  \cfrac{\partial \mathbf{f}}{\partial \text{x}}
  & = \xi_{x} \cfrac{\partial\mathbf{\hat{f}}}{\partial \xi}
  + \eta_{x} \cfrac{\partial\mathbf{\hat{f}}}{\partial \eta}
  + \zeta_{x} \cfrac{\partial\mathbf{\hat{f}}}{\partial \zeta}
  + \tau_{x} \cfrac{\partial\mathbf{\hat{f}}}{\partial \tau}\\  
  
- 
 
.. math::
  \cfrac{\partial \mathbf{g}}{\partial \text{y}}
  = \xi_{y} \cfrac{\partial\mathbf{\hat{g}}}{\partial \xi}
  + \eta_{y} \cfrac{\partial\mathbf{\hat{g}}}{\partial \eta}
  + \zeta_{y} \cfrac{\partial\mathbf{\hat{g}}}{\partial \zeta}
  + \tau_{y} \cfrac{\partial\mathbf{\hat{g}}}{\partial \tau}\\  
  
- 
 
.. math::
  \cfrac{\partial \mathbf{h}}{\partial \text{z}}
  = \xi_{z} \cfrac{\partial\mathbf{\hat{h}}}{\partial \xi}
  + \eta_{z} \cfrac{\partial\mathbf{\hat{h}}}{\partial \eta}
  + \zeta_{z} \cfrac{\partial\mathbf{\hat{h}}}{\partial \zeta}
  + \tau_{z} \cfrac{\partial\mathbf{\hat{h}}}{\partial \tau}\\ 
  

multiply by :math:`J^{-1}` to get:
 
.. math::
  J^{-1}(\cfrac{\partial \mathbf{q}}{\partial \text{t}}+
  \cfrac{\partial \mathbf{f}}{\partial \text{x}}+
  \cfrac{\partial \mathbf{g}}{\partial \text{y}}+
  \cfrac{\partial \mathbf{h}}{\partial \text{z}})=0\\
  
-
  
.. math::
  \begin{align}
  0 & = \hat\xi_{t} \cfrac{\partial \mathbf{\hat{q}}}{\partial \xi}
  + \hat\eta_{t} \cfrac{\partial\mathbf{\hat{q}}}{\partial \eta}
  + \hat\zeta_{t} \cfrac{\partial\mathbf{\hat{q}}}{\partial \zeta}
  + \hat\tau_{t} \cfrac{\partial\mathbf{\hat{q}}}{\partial \tau}\\
  &+  \hat\xi_{x} \cfrac{\partial\mathbf{\hat{f}}}{\partial \xi}
  + \hat\eta_{x} \cfrac{\partial\mathbf{\hat{f}}}{\partial \eta}
  + \hat\zeta_{x} \cfrac{\partial\mathbf{\hat{f}}}{\partial \zeta}
  + \hat\tau_{x} \cfrac{\partial\mathbf{\hat{f}}}{\partial \tau}\\ 
  &+ \hat\xi_{y} \cfrac{\partial\mathbf{\hat{g}}}{\partial \xi}
  + \hat\eta_{y} \cfrac{\partial\mathbf{\hat{g}}}{\partial \eta}
  + \hat\zeta_{y} \cfrac{\partial\mathbf{\hat{g}}}{\partial \zeta}
  + \hat\tau_{y} \cfrac{\partial\mathbf{\hat{g}}}{\partial \tau}\\
  &+  \hat\xi_{z} \cfrac{\partial\mathbf{\hat{h}}}{\partial \xi}
  + \hat\eta_{z} \cfrac{\partial\mathbf{\hat{h}}}{\partial \eta}
  + \hat\zeta_{z} \cfrac{\partial\mathbf{\hat{h}}}{\partial \zeta}
  + \hat\tau_{z} \cfrac{\partial\mathbf{\hat{h}}}{\partial \tau}\\ 
  \end{align}  
  
where

.. math::
  \begin{align}
  (\hat\xi_{x},\hat\xi_{y},\hat\xi_{z},\hat\xi_{t}) & = J^{-1}(\xi_{x},\xi_{y},\xi_{z},\xi_{t})\\
  (\hat\eta_{x},\hat\eta_{y},\hat\eta_{z},\hat\eta_{t}) & = J^{-1}(\eta_{x},\eta_{y},\eta_{z},\eta_{t})\\
  (\hat\zeta_{x},\hat\zeta_{y},\hat\zeta_{z},\hat\zeta_{t}) & = J^{-1}(\zeta_{x},\zeta_{y},\zeta_{z},\zeta_{t})\\
  (\hat\tau_{x},\hat\tau_{y},\hat\tau_{z},\hat\tau_{t}) & = J^{-1}(\tau_{x},\tau_{y},\tau_{z},\tau_{t})  = J^{-1}(0,0,0,1)\\
  \end{align}
  
then
  
.. math::
  \begin{align}
  0 & = \cfrac{\partial (\hat\xi_{t}\mathbf{\hat{q}})}{\partial \xi}
  + \cfrac{\partial(\hat\eta_{t}\mathbf{\hat{q}})}{\partial \eta}
  +\cfrac{\partial(\hat\zeta_{t}\mathbf{\hat{q}})}{\partial \zeta}
  + \cfrac{\partial(\hat\tau_{t}\mathbf{\hat{q}})}{\partial \tau}\\
  &-\mathbf{\hat{q}}
  (\cfrac{\partial (\hat\tau_{t}\equiv J^{-1})}{\partial \tau}
  +\cfrac{\partial (\hat\xi_{t})}{\partial \xi}
  +\cfrac{\partial(\hat\eta_{t})}{\partial \eta}
  +\cfrac{\partial(\hat\zeta_{t})}{\partial \zeta}
  )\\
  &+  \cfrac{\partial(\hat\xi_{x}\mathbf{\hat{f}})}{\partial \xi}
  + \cfrac{\partial(\hat\eta_{x}\mathbf{\hat{f}})}{\partial \eta}
  + \cfrac{\partial(\hat\zeta_{x}\mathbf{\hat{f}})}{\partial \zeta}
  + \cfrac{\partial(\hat\tau_{x}\mathbf{\hat{f}}\equiv 0)}{\partial \tau}\\
  &- \mathbf{\hat{f}}(\cfrac{\partial(\hat\xi_{x})}{\partial \xi}
  +\cfrac{\partial(\hat\eta_{x})}{\partial \eta}
  +\cfrac{\partial(\hat\zeta_{x})}{\partial \zeta}
  )\\
  &+ \cfrac{\partial(\hat\xi_{y} \mathbf{\hat{g}})}{\partial \xi}
  + \cfrac{\partial(\hat\eta_{y}\mathbf{\hat{g}})}{\partial \eta}
  + \cfrac{\partial(\hat\zeta_{y}\mathbf{\hat{g}})}{\partial \zeta}
  + \cfrac{\partial(\hat\tau_{y}\mathbf{\hat{g}}\equiv 0)}{\partial \tau}\\
  &-\mathbf{\hat{g}}( \cfrac{\partial(\hat\xi_{y})}{\partial \xi}
  + \cfrac{\partial(\hat\eta_{y})}{\partial \eta}
  + \cfrac{\partial(\hat\zeta_{y})}{\partial \zeta})\\
  &+  \cfrac{\partial(\hat\xi_{z}\mathbf{\hat{h}})}{\partial \xi}
  + \cfrac{\partial(\hat\eta_{z}\mathbf{\hat{h}})}{\partial \eta}
  + \cfrac{\partial(\hat\zeta_{z}\mathbf{\hat{h}})}{\partial \zeta}
  + \cfrac{\partial(\hat\tau_{z}\mathbf{\hat{h}}\equiv 0)}{\partial \tau}\\ 
  &-\mathbf{\hat{h}}(\cfrac{\partial(\hat\xi_{z})}{\partial \xi}
  + \cfrac{\partial(\hat\eta_{z})}{\partial \eta}
  + \cfrac{\partial(\hat\zeta_{z})}{\partial \zeta})
  \end{align}
  
-
  
.. math::
  \begin{align}
  &\cfrac{\partial(\hat\xi_{x})}{\partial \xi}
  +\cfrac{\partial(\hat\eta_{x})}{\partial \eta}
  +\cfrac{\partial(\hat\zeta_{x})}{\partial \zeta}\\
  & = \cfrac{\partial({y}_{\eta}{z}_{\zeta}-{z}_{\eta}{y}_{\zeta})}{\partial \xi}
  +\cfrac{\partial({z}_{\xi}{y}_{\zeta}-{y}_{\xi}{z}_{\zeta})}{\partial \eta}
  +\cfrac{\partial({y}_{\xi}{z}_{\eta}-{z}_{\xi}{y}_{\eta})}{\partial \zeta}\\
  &=({y}_{\xi\eta}{z}_{\zeta}+{y}_{\eta}{z}_{\xi\zeta}-{z}_{\xi\eta}{y}_{\zeta}-{z}_{\eta}{y}_{\xi\zeta})\\
  &+({z}_{\xi\eta}{y}_{\zeta}+{z}_{\xi}{y}_{\eta\zeta}-{y}_{\xi\eta}{z}_{\zeta}-{y}_{\xi}{z}_{\eta\zeta})\\
  &+({y}_{\xi\zeta}{z}_{\eta}+{y}_{\xi}{z}_{\eta\zeta}-{z}_{\xi\zeta}{y}_{\eta}-{z}_{\xi}{y}_{\eta\zeta})\\
  &={y}_{\xi}({z}_{\eta\zeta}-{z}_{\eta\zeta})+{y}_{\eta}({z}_{\xi\zeta}-{z}_{\xi\zeta})+{y}_{\zeta}({z}_{\xi\eta}-{z}_{\xi\eta})\\
  &+{z}_{\xi}({y}_{\eta\zeta}-{y}_{\eta\zeta})+{z}_{\eta}(-{y}_{\xi\zeta}+{y}_{\xi\zeta})+{z}_{\zeta}({y}_{\xi\eta}-{y}_{\xi\eta})\\
  &=0
  \end{align}

-
  
.. math::
  \begin{align}
  &\cfrac{\partial(\hat\xi_{y})}{\partial \xi}
  + \cfrac{\partial(\hat\eta_{y})}{\partial \eta}
  + \cfrac{\partial(\hat\zeta_{y})}{\partial \zeta}\\
  &=\cfrac{\partial({z}_{\eta}{x}_{\zeta}-{x}_{\eta}{z}_{\zeta})}{\partial \xi}
  +\cfrac{\partial({x}_{\xi}{z}_{\zeta}-{z}_{\xi}{x}_{\zeta})}{\partial \eta}
  +\cfrac{\partial({z}_{\xi}{x}_{\eta}-{x}_{\xi}{z}_{\eta})}{\partial \zeta}\\
  &={z}_{\xi\eta}{x}_{\zeta}+{z}_{\eta}{x}_{\xi\zeta}-{x}_{\xi\eta}{z}_{\zeta}-{x}_{\eta}{z}_{\xi\zeta}\\
  &+{x}_{\xi\eta}{z}_{\zeta}+{x}_{\xi}{z}_{\eta\zeta}-{z}_{\xi\eta}{x}_{\zeta}-{z}_{\xi}{x}_{\eta\zeta}\\
  &+{z}_{\xi\zeta}{x}_{\eta}+{z}_{\xi}{x}_{\eta\zeta}-{x}_{\xi\zeta}{z}_{\eta}-{x}_{\xi}{z}_{\eta\zeta}\\
  &={x}_{\xi}({z}_{\eta\zeta}-{z}_{\eta\zeta})+{x}_{\eta}(-{z}_{\xi\zeta}+{z}_{\xi\zeta})+{x}_{\zeta}({z}_{\xi\eta}-{z}_{\xi\eta})\\
  &+{z}_{\xi}(-{x}_{\eta\zeta}+{x}_{\eta\zeta})+{z}_{\eta}({x}_{\xi\zeta}-{x}_{\xi\zeta})+{z}_{\zeta}(-{x}_{\xi\eta}+{x}_{\xi\eta})\\
  &=0
  \end{align}  
  
-
  
.. math::  
  \begin{align}
  &\cfrac{\partial(\hat\xi_{z})}{\partial \xi}
    + \cfrac{\partial(\hat\eta_{z})}{\partial \eta}
    + \cfrac{\partial(\hat\zeta_{z})}{\partial \zeta}\\
  &=\cfrac{\partial({x}_{\eta}{y}_{\zeta}-{y}_{\eta}{x}_{\zeta})}{\partial \xi}
    + \cfrac{\partial({y}_{\xi}{x}_{\zeta}-{x}_{\xi}{y}_{\zeta})}{\partial \eta}
    + \cfrac{\partial({x}_{\xi}{y}_{\eta}-{y}_{\xi}{x}_{\eta})}{\partial \zeta}\\
  &={x}_{\xi\eta}{y}_{\zeta}+{x}_{\eta}{y}_{\xi\zeta}-{y}_{\xi\eta}{x}_{\zeta}-{y}_{\eta}{x}_{\xi\zeta}\\
  &+{y}_{\xi\eta}{x}_{\zeta}+{y}_{\xi}{x}_{\eta\zeta}-{x}_{\xi\eta}{y}_{\zeta}-{x}_{\xi}{y}_{\eta\zeta}\\
  &+{x}_{\xi\zeta}{y}_{\eta}+{x}_{\xi}{y}_{\eta\zeta}-{y}_{\xi\zeta}{x}_{\eta}-{y}_{\xi}{x}_{\eta\zeta}\\
  &={x}_{\xi}(-{y}_{\eta\zeta}+{y}_{\eta\zeta})+{x}_{\eta}({y}_{\xi\zeta}-{y}_{\xi\zeta})+{x}_{\zeta}(-{y}_{\xi\eta}+{y}_{\xi\eta})\\
  &+{y}_{\xi}({x}_{\eta\zeta}-{x}_{\eta\zeta})+{y}_{\eta}(-{x}_{\xi\zeta}+{x}_{\xi\zeta})+{y}_{\zeta}({x}_{\xi\eta}-{x}_{\xi\eta})\\
  &=0
  \end{align}  
  
-
  
.. math::  
  \begin{align}
  &\cfrac{\partial (\hat\xi_{t})}{\partial \xi}
   +\cfrac{\partial(\hat\eta_{t})}{\partial \eta}
   +\cfrac{\partial(\hat\zeta_{t})}{\partial \zeta}\\
  &=\cfrac{\partial (-x_{\tau}\hat\xi_{x}-y_{\tau}\hat\xi_{y}-z_{\tau}\hat\xi_{z})}{\partial \xi}
   +\cfrac{\partial(-x_{\tau}\hat\eta_{x}-y_{\tau}\hat\eta_{y}-z_{\tau}\hat\eta_{z})}{\partial \eta}
   +\cfrac{\partial(-x_{\tau}\hat\zeta_{x}-y_{\tau}\hat\zeta_{y}-z_{\tau}\hat\zeta_{z})}{\partial \zeta}\\
  &=-x_{\tau}(\cfrac{\partial\hat\xi_{x}}{\partial \xi}
             +\cfrac{\partial\hat\eta_{x}}{\partial \eta}
             +\cfrac{\partial\hat\zeta_{x}}{\partial \zeta})
  -y_{\tau}(\cfrac{\partial\hat\xi_{y}}{\partial \xi}
             +\cfrac{\partial\hat\eta_{y}}{\partial \eta}
             +\cfrac{\partial\hat\zeta_{y}}{\partial \zeta})
  -z_{\tau}(\cfrac{\partial\hat\xi_{z}}{\partial \xi}
             +\cfrac{\partial\hat\eta_{z}}{\partial \eta}
             +\cfrac{\partial\hat\zeta_{z}}{\partial \zeta})\\
  &(-\hat\xi_{x}x_{\xi\tau}-\hat\xi_{y}y_{\xi\tau}-\hat\xi_{z}z_{\xi\tau})
  +(-\hat\eta_{x}x_{\eta\tau}-\hat\eta_{y}y_{\eta\tau}-\hat\eta_{z}z_{\eta\tau})
  +(-\hat\zeta_{x}x_{\zeta\tau}-\hat\zeta_{y}y_{\zeta\tau}-\hat\zeta_{z}z_{\zeta\tau})\\
  &=(-\hat\xi_{x}x_{\xi\tau}-\hat\xi_{y}y_{\xi\tau}-\hat\xi_{z}z_{\xi\tau})
  +(-\hat\eta_{x}x_{\eta\tau}-\hat\eta_{y}y_{\eta\tau}-\hat\eta_{z}z_{\eta\tau})
  +(-\hat\zeta_{x}x_{\zeta\tau}-\hat\zeta_{y}y_{\zeta\tau}-\hat\zeta_{z}z_{\zeta\tau})
  \end{align}
  
-
  
.. math::   
  \begin{align}
  J^{-1} 
  & = {x}_{\xi}({y}_{\eta}{z}_{\zeta}-{z}_{\eta}{y}_{\zeta})
    - {y}_{\xi}({x}_{\eta}{z}_{\zeta}-{z}_{\eta}{x}_{\zeta})
    + {z}_{\xi}({x}_{\eta}{y}_{\zeta}-{y}_{\eta}{x}_{\zeta}) \\
  & = {x}_{\eta}({y}_{\zeta}{z}_{\xi}-{y}_{\xi}{z}_{\zeta})
    + {y}_{\eta}({x}_{\xi}{z}_{\zeta}-{z}_{\xi}{x}_{\zeta})
    + {z}_{\eta}({x}_{\zeta}{y}_{\xi}-{x}_{\xi}{y}_{\zeta}) \\
  & = {x}_{\zeta}({y}_{\xi}{z}_{\eta}-{z}_{\xi}{y}_{\zeta})
    + {y}_{\zeta}({z}_{\xi}{x}_{\eta}-{x}_{\xi}{z}_{\eta})
    + {z}_{\zeta}({x}_{\xi}{y}_{\eta}-{y}_{\xi}{x}_{\eta}) \\
  \end{align}
  
-
  
.. math::   
  \begin{align}
  J^{-1} 
  & = {x}_{\xi}({y}_{\eta}{z}_{\zeta}-{z}_{\eta}{y}_{\zeta})
    + {y}_{\xi}({z}_{\eta}{x}_{\zeta}-{x}_{\eta}{z}_{\zeta})
    + {z}_{\xi}({x}_{\eta}{y}_{\zeta}-{y}_{\eta}{x}_{\zeta}) \\
  & = {x}_{\eta}({y}_{\zeta}{z}_{\xi}-{z}_{\zeta}{y}_{\xi})
    + {y}_{\eta}({z}_{\zeta}{x}_{\xi}-{x}_{\zeta}{z}_{\xi})
    + {z}_{\eta}({x}_{\zeta}{y}_{\xi}-{y}_{\zeta}{x}_{\xi}) \\
  & = {x}_{\zeta}({y}_{\xi}{z}_{\eta}-{z}_{\xi}{y}_{\zeta})
    + {y}_{\zeta}({z}_{\xi}{x}_{\eta}-{x}_{\xi}{z}_{\eta})
    + {z}_{\zeta}({x}_{\xi}{y}_{\eta}-{y}_{\xi}{x}_{\eta}) \\
  \end{align}  
  
-
  
.. math:: 
  \begin{align}
  J^{-1} & = x_{\xi}\hat\xi_{x}+y_{\xi}\hat\xi_{y}+z_{\xi}\hat\xi_{z}\\
         & = x_{\eta}\hat\eta_{x}+y_{\eta}\hat\eta_{y}+z_{\eta}\hat\eta_{z}\\
         & = x_{\zeta}\hat\zeta_{x}+y_{\zeta}\hat\zeta_{y}+z_{\zeta}\hat\zeta_{z}\\
  \end{align}
  
-
  
.. math:: 
  \begin{align}
  1 & = x_{\xi}\xi_{x}+y_{\xi}\xi_{y}+z_{\xi}\xi_{z}\\
         & = x_{\eta}\eta_{x}+y_{\eta}\eta_{y}+z_{\eta}\eta_{z}\\
         & = x_{\zeta}\zeta_{x}+y_{\zeta}\zeta_{y}+z_{\zeta}\zeta_{z}\\
  \end{align} 
  
-
  
.. math:: 
  \begin{align}
  \xi_{t} &= -x_{\tau}\xi_{x}-y_{\tau}\xi_{y}-z_{\tau}\xi_{z}\\
  \eta_{t} &= -x_{\tau}\eta_{x}-y_{\tau}\eta_{y}-z_{\tau}\eta_{z}\\
  \zeta_{t} &= -x_{\tau}\zeta_{x}-y_{\tau}\zeta_{y}-z_{\tau}\zeta_{z}\\
  \end{align} 
  
-
  
.. math:: 
  \begin{align}
  \hat\xi_{t} &= -x_{\tau}\hat\xi_{x}-y_{\tau}\hat\xi_{y}-z_{\tau}\hat\xi_{z}\\
  \hat\eta_{t} &= -x_{\tau}\hat\eta_{x}-y_{\tau}\hat\eta_{y}-z_{\tau}\hat\eta_{z}\\
  \hat\zeta_{t} &= -x_{\tau}\hat\zeta_{x}-y_{\tau}\hat\zeta_{y}-z_{\tau}\hat\zeta_{z}\\
  \end{align}  
  
-
  
.. math:: 
  \begin{align}
  &\cfrac{\partial (\hat\xi_{t})}{\partial \xi}
   +\cfrac{\partial(\hat\eta_{t})}{\partial \eta}
   +\cfrac{\partial(\hat\zeta_{t})}{\partial \zeta}\\
  &=(-\hat\xi_{x}x_{\xi\tau}-\hat\xi_{y}y_{\xi\tau}-\hat\xi_{z}z_{\xi\tau})
  +(-\hat\eta_{x}x_{\eta\tau}-\hat\eta_{y}y_{\eta\tau}-\hat\eta_{z}z_{\eta\tau})
  +(-\hat\zeta_{x}x_{\zeta\tau}-\hat\zeta_{y}y_{\zeta\tau}-\hat\zeta_{z}z_{\zeta\tau})
  \end{align}
  
-
  
.. math::
  \begin{array}{c}
  \hat{\xi}_{x}=({y}_{\eta}{z}_{\zeta}-{z}_{\eta}{y}_{\zeta})\\
  \hat{\xi}_{y}=({z}_{\eta}{x}_{\zeta}-{x}_{\eta}{z}_{\zeta})\\
  \hat{\xi}_{z}=({x}_{\eta}{y}_{\zeta}-{y}_{\eta}{x}_{\zeta})\\
  \hat{\xi}_{x\tau}=({y}_{\eta\tau}{z}_{\zeta}+{y}_{\eta}{z}_{\zeta\tau}-{z}_{\eta\tau}{y}_{\zeta}-{z}_{\eta}{y}_{\zeta\tau})\\
  \hat{\xi}_{y\tau}=({z}_{\eta\tau}{x}_{\zeta}+{z}_{\eta}{x}_{\zeta\tau}-{x}_{\eta\tau}{z}_{\zeta}-{x}_{\eta}{z}_{\zeta\tau})\\
  \hat{\xi}_{z\tau}=({x}_{\eta\tau}{y}_{\zeta}+{x}_{\eta}{y}_{\zeta\tau}-{y}_{\eta\tau}{x}_{\zeta}-{y}_{\eta}{x}_{\zeta\tau})\\
  x_{\xi}\hat{\xi}_{x\tau}=x_{\xi}({y}_{\eta\tau}{z}_{\zeta}+{y}_{\eta}{z}_{\zeta\tau}-{z}_{\eta\tau}{y}_{\zeta}-{z}_{\eta}{y}_{\zeta\tau})\\
  y_{\xi}\hat{\xi}_{y\tau}=y_{\xi}({z}_{\eta\tau}{x}_{\zeta}+{z}_{\eta}{x}_{\zeta\tau}-{x}_{\eta\tau}{z}_{\zeta}-{x}_{\eta}{z}_{\zeta\tau})\\
  z_{\xi}\hat{\xi}_{z\tau}=z_{\xi}({x}_{\eta\tau}{y}_{\zeta}+{x}_{\eta}{y}_{\zeta\tau}-{y}_{\eta\tau}{x}_{\zeta}-{y}_{\eta}{x}_{\zeta\tau})\\
  \end{array} 
  
-
  
.. math::
  \begin{align}
  x_{\xi}\hat{\xi}_{x\tau}+y_{\xi}\hat{\xi}_{y\tau}+z_{\xi}\hat{\xi}_{z\tau}
  & = x_{\xi}({y}_{\eta\tau}{z}_{\zeta}+{y}_{\eta}{z}_{\zeta\tau}-{z}_{\eta\tau}{y}_{\zeta}-{z}_{\eta}{y}_{\zeta\tau})\\
  & + y_{\xi}({z}_{\eta\tau}{x}_{\zeta}+{z}_{\eta}{x}_{\zeta\tau}-{x}_{\eta\tau}{z}_{\zeta}-{x}_{\eta}{z}_{\zeta\tau})\\
  & + z_{\xi}({x}_{\eta\tau}{y}_{\zeta}+{x}_{\eta}{y}_{\zeta\tau}-{y}_{\eta\tau}{x}_{\zeta}-{y}_{\eta}{x}_{\zeta\tau})\\
  &={x}_{\eta\tau}({y}_{\zeta}z_{\xi}-{z}_{\zeta}y_{\xi})
  +{y}_{\eta\tau}(x_{\xi}{z}_{\zeta}-z_{\xi}{x}_{\zeta})
  +{z}_{\eta\tau}(y_{\xi}{x}_{\zeta}-x_{\xi}{y}_{\zeta})\\
  &+{x}_{\zeta\tau}(y_{\xi}{z}_{\eta}-z_{\xi}{y}_{\eta})
  +{y}_{\zeta\tau}(z_{\xi}{x}_{\eta}-x_{\xi}{z}_{\eta})
  +{z}_{\zeta\tau}(x_{\xi}{y}_{\eta}-y_{\xi}{x}_{\eta})
  \end{align}
  
-
  
.. math::
  \begin{align}
  \hat\eta_{x} & = +({z}_{\xi}{y}_{\zeta}-{y}_{\xi}{z}_{\zeta})\\
  \hat\eta_{y} & = +({x}_{\xi}{z}_{\zeta}-{z}_{\xi}{x}_{\zeta})\\
  \hat\eta_{z} & = +({y}_{\xi}{x}_{\zeta}-{x}_{\xi}{y}_{\zeta})\\
  \end{align} 

-
  
.. math::
  \begin{align}
  \hat\zeta_{x} & = +({y}_{\xi}{z}_{\eta}-{z}_{\xi}{y}_{\eta})\\
  \hat\zeta_{y} & = +({z}_{\xi}{x}_{\eta}-{x}_{\xi}{z}_{\eta})\\
  \hat\zeta_{z} & = +({x}_{\xi}{y}_{\eta}-{y}_{\xi}{x}_{\eta})\\
  \end{align}  
  
-
  
.. math::  
  \begin{align}
  x_{\xi}\hat{\xi}_{x\tau}+y_{\xi}\hat{\xi}_{y\tau}+z_{\xi}\hat{\xi}_{z\tau}
  &={x}_{\eta\tau}({y}_{\zeta}z_{\xi}-{z}_{\zeta}y_{\xi})
  +{y}_{\eta\tau}(x_{\xi}{z}_{\zeta}-z_{\xi}{x}_{\zeta})
  +{z}_{\eta\tau}(y_{\xi}{x}_{\zeta}-x_{\xi}{y}_{\zeta})\\
  &+{x}_{\zeta\tau}(y_{\xi}{z}_{\eta}-z_{\xi}{y}_{\eta})
  +{y}_{\zeta\tau}(z_{\xi}{x}_{\eta}-x_{\xi}{z}_{\eta})
  +{z}_{\zeta\tau}(x_{\xi}{y}_{\eta}-y_{\xi}{x}_{\eta})\\
  &={x}_{\eta\tau}(\hat\eta_{x})
  +{y}_{\eta\tau}(\hat\eta_{y})
  +{z}_{\eta\tau}(\hat\eta_{z})\\
  &+{x}_{\zeta\tau}(\hat\zeta_{x})
  +{y}_{\zeta\tau}(\hat\zeta_{y})
  +{z}_{\zeta\tau}(\hat\zeta_{z})\\
  \end{align}  
  
-
  
.. math::
  \begin{align}
  \cfrac{\partial J^{-1}}{\partial \tau}=\cfrac{\partial (x_{\xi}\hat\xi_{x}+y_{\xi}\hat\xi_{y}+z_{\xi}\hat\xi_{z})}{\partial \tau}
  =(x_{\xi\tau}\hat\xi_{x}+y_{\xi\tau}\hat\xi_{y}+z_{\xi\tau}\hat\xi_{z})
  +(x_{\xi}\hat\xi_{x\tau}+y_{\xi}\hat\xi_{y\tau}+z_{\xi}\hat\xi_{z\tau})\\
  \end{align} 
  
-
  
.. math::
  \begin{align}
  &\cfrac{\partial J^{-1}}{\partial \tau}+  \cfrac{\partial (\hat\xi_{t})}{\partial \xi}
   +\cfrac{\partial(\hat\eta_{t})}{\partial \eta}
   +\cfrac{\partial(\hat\zeta_{t})}{\partial \zeta}\\
  & = (x_{\xi\tau}\hat\xi_{x}+y_{\xi\tau}\hat\xi_{y}+z_{\xi\tau}\hat\xi_{z})
  +(x_{\xi}\hat\xi_{x\tau}+y_{\xi}\hat\xi_{y\tau}+z_{\xi}\hat\xi_{z\tau})\\
  & + (-\hat\xi_{x}x_{\xi\tau}-\hat\xi_{y}y_{\xi\tau}-\hat\xi_{z}z_{\xi\tau})
  +(-\hat\eta_{x}x_{\eta\tau}-\hat\eta_{y}y_{\eta\tau}-\hat\eta_{z}z_{\eta\tau})
  +(-\hat\zeta_{x}x_{\zeta\tau}-\hat\zeta_{y}y_{\zeta\tau}-\hat\zeta_{z}z_{\zeta\tau})\\
  &=(x_{\xi}\hat\xi_{x\tau}+y_{\xi}\hat\xi_{y\tau}+z_{\xi}\hat\xi_{z\tau})+(-\hat\eta_{x}x_{\eta\tau}-\hat\eta_{y}y_{\eta\tau}-\hat\eta_{z}z_{\eta\tau})
  +(-\hat\zeta_{x}x_{\zeta\tau}-\hat\zeta_{y}y_{\zeta\tau}-\hat\zeta_{z}z_{\zeta\tau})\\
  &=0
  \end{align}
  
Therefor the general curvilinear equation can now be expressed:
  
.. math::
  \begin{align}
  \cfrac{\partial(J^{-1}\hat{\mathbf{q}})}{\partial t}
  &+\cfrac{\partial}{\partial \xi}[\hat\xi_{t}\hat{\mathbf{q}}+\hat\xi_{x}\hat{\mathbf{f}}+\hat\xi_{y}\hat{\mathbf{g}}+\hat\xi_{z}\hat{\mathbf{h}}]\\
  &+\cfrac{\partial}{\partial \eta}[\hat\eta_{t}\hat{\mathbf{q}}+\hat\eta_{x}\hat{\mathbf{f}}+\hat\eta_{y}\hat{\mathbf{g}}+\hat\eta_{z}\hat{\mathbf{h}}]\\
  &+\cfrac{\partial}{\partial \zeta}[\hat\zeta_{t}\hat{\mathbf{q}}+\hat\zeta_{x}\hat{\mathbf{f}}+\hat\zeta_{y}\hat{\mathbf{g}}+\hat\zeta_{z}\hat{\mathbf{h}}]\\
  &=0
  \end{align}
  
or in the more compact form:

.. math::
  \begin{align}
  \cfrac{\partial(J^{-1}{\mathbf{Q}})}{\partial \tau}
  +\cfrac{\partial{\mathbf{F}}}{\partial \xi}
  +\cfrac{\partial{\mathbf{G}}}{\partial \eta}
  +\cfrac{\partial{\mathbf{H}}}{\partial \zeta}
  =0
  \end{align}  
  
where the vector of conserved variables, :math:`\mathbf{Q}`, and the vector of inviscid flux terms, :math:`\mathbf{F}`,
:math:`\mathbf{G}` and :math:`\mathbf{H}`, are:

.. math::
  \mathbf{Q}=\mathbf{\hat{q}}=\begin{bmatrix}
  \rho\\ \rho u \\ \rho v \\ \rho w\\ \rho E
  \end{bmatrix}
  
-
  
.. math::
  \mathbf{F}=\hat\xi_{t}\hat{\mathbf{q}}+\hat\xi_{x}\hat{\mathbf{f}}+\hat\xi_{y}\hat{\mathbf{g}}+\hat\xi_{z}\hat{\mathbf{h}}=\begin{bmatrix}
  \rho U\\
  \rho u U + \hat\xi_{x}p\\
  \rho v U + \hat\xi_{y}p\\
  \rho w U + \hat\xi_{z}p\\
  \rho H U - \hat\xi_{t}p\\
  \end{bmatrix}
  
-
  
.. math::
  \mathbf{G}=\hat\eta_{t}\hat{\mathbf{q}}+\hat\eta_{x}\hat{\mathbf{f}}+\hat\eta_{y}\hat{\mathbf{g}}+\hat\eta_{z}\hat{\mathbf{h}}=\begin{bmatrix}
  \rho V\\
  \rho u V + \hat\eta_{x}p\\
  \rho v V + \hat\eta_{y}p\\
  \rho w V + \hat\eta_{z}p\\
  \rho H V - \hat\eta_{t}p\\
  \end{bmatrix}  
  
-
  
.. math::
  \mathbf{H}=\hat\zeta_{t}\hat{\mathbf{q}}+\hat\zeta_{x}\hat{\mathbf{f}}+\hat\zeta_{y}\hat{\mathbf{g}}+\hat\zeta_{z}\hat{\mathbf{h}}=\begin{bmatrix}
  \rho W\\
  \rho u W + \hat\zeta_{x}p\\
  \rho v W + \hat\zeta_{y}p\\
  \rho w W + \hat\zeta_{z}p\\
  \rho H W - \hat\zeta_{t}p\\
  \end{bmatrix}  
  
:math:`U`, :math:`V` and :math:`W` are the contravariant velocity components in the :math:`\xi`, :math:`\eta` and :math:`\zeta`
directions and are defined as:

.. math::
  \begin{align}
  U & = \hat\xi_{x}u+\hat\xi_{y}v+\hat\xi_{z}w+\hat\xi_{t}\\
  V & = \hat\eta_{x}u+\hat\eta_{y}v+\hat\eta_{z}w+\hat\eta_{t}\\
  W & = \hat\zeta_{x}u+\hat\zeta_{y}v+\hat\zeta_{z}w+\hat\zeta_{t}\\
  \end{align}
