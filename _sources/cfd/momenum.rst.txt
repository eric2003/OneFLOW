Momentum Analysis Of Flow Systems
==================================

- Fluid Mechanics Fundamentals and Applications, Fourth Edition, Yunus Cengel, John Cimbala

pp-250

Forces Acting On A Control Volume
``````````````````````````````````````````

Total force acting on control volume:

.. math::
  \sum \vec{F}=\sum \vec{F}_{\text {body }}+\sum \vec{F}_{\text {surface }}
  =\int_{\mathrm{CV}} \rho \vec{g}  d V+\int_{\mathrm{CS}} \sigma_{ij}\cdot \vec{n} d A
  
  
Total force:

.. math::
  \underbrace{\sum \vec{F}}_{\text {total force }}=\underbrace{\sum \vec{F}_{\text {gravity }}}_{\text {body force }}+\underbrace{\sum \vec{F}_{\text {pressure }}+\sum \vec{F}_{\text {viscous }}+\sum \vec{F}_{\text {other }}}_{\text {surface forces }}  
  
The Linear Momentum Equation:

The general form of the linear momentum equation that applies to
fixed, moving, or deforming control volumes is

.. math::
  \frac{d ({m}  \vec{V})_{\mathrm{sys}}}{d t}=\frac{d}{d t} \int_{\mathrm{CV}} \rho \vec{V} d V+\int_{\mathrm{CS}} (\rho\vec{V}) (\vec{V}_{r} \cdot \vec{n}) d A  

-

.. math::
  \frac{d ({m}  \vec{V})_{\mathrm{sys}}}{d t}=\frac{d}{d t} \int_{\mathrm{CV}} \rho \vec{V} d V+\int_{\mathrm{CS}} (\rho\vec{V}) ((\vec{V}-\vec{V}_{\text{CS}}) \cdot \vec{n}) d A
  
General:

.. math::
  \sum \vec{F}=\frac{d}{d t} \int_{\mathrm{CV}} \rho \vec{V} d V+\int_{\mathrm{CS}} (\rho\vec{V}) ((\vec{V}-\vec{V}_{\text{CS}}) \cdot \vec{n}) d A  

which is stated in words as

.. math::
  \begin{pmatrix}
  \text{The sum of all}\\
  \text{external forces}\\
  \text{acting on a CV}\\
  \end{pmatrix}
  =\begin{pmatrix}
  \text{The time rate of change}\\
  \text{of the linear momentum}\\
  \text{of the contents of the CV}\\
  \end{pmatrix}
  +\begin{pmatrix}
  \text{The net flow rate of}\\
  \text{the linear momentum out of the}\\
  \text{control surface by mass flow}\\
  \end{pmatrix}
  
Here :math:`\vec{V}_{\text{r}}=\vec{V}-\vec{V}_{\text{CS}}` is the fluid velocity relative to the control surface (for
use in mass flow rate calculations at all locations where the fluid crosses the
control surface), and :math:`\vec{V}` is the fluid velocity as viewed from an inertial reference frame. The product 
:math:`\rho(\vec{V}_{r}\cdot\vec{n})dA` represents the mass flow rate through
area element :math:`dA` into or out of the control volume.  

For a fixed control volume (no motion or deformation of the control volume),
:math:`\vec{V}_{\text{r}}=\vec{V}` and the linear momentum equation becomes

Fixed CV:

.. math::
  \sum \vec{F}=\frac{d}{d t} \int_{\mathrm{CV}} \rho \vec{V} d V+\int_{\mathrm{CS}} (\rho\vec{V}) (\vec{V} \cdot \vec{n}) d A    
  
The Differential Linear Momentum Equation-Cauchy's Equation
`````````````````````````````````````````````````````````````

.. math::
  \sum \vec{F}=\int_{\mathrm{CV}} \rho \vec{g}  d V+\int_{\mathrm{CS}} \sigma_{ij}\cdot \vec{n} d A
  =\int_{\mathrm{CV}}\frac{\partial }{\partial t} (\rho \vec{V}) d V+\int_{\mathrm{CS}} (\rho\vec{V}) (\vec{V} \cdot \vec{n}) d A  
  
The Divergence Theorem-tensors  
------------------------------------
Consider an arbitrary differentiable tensor field :math:`T_{ij\cdots k}(\mathbf{x},t)` defined in some finite region of
physical space. Let :math:`S` be a closed surface bounding a volume :math:`V` in this space, and let the
outward normal to :math:`S` be :math:`\mathbf{n}`. The divergence theorem of Gauss then states that

.. math::
  \int\limits_{S}T_{ij\cdots k}\cdot{n}_{k}\text{d}S=\int\limits_{V}\cfrac{\partial T_{ij\cdots k}}{\partial {x}_{k}}\text{d}V\\  
  
For a second order tensor,

.. math::
  \int\limits_{S}\mathbf{T}\cdot\mathbf{n}\text{d}S=\int\limits_{V}\text{div }\mathbf{T}\text{d}V\\

-
  
.. math::
  \int\limits_{S}{T}_{ij}\cdot{n}_{j}\text{d}S=\int\limits_{V}\cfrac{\partial {T}_{ij}}{\partial {x}_{j}}\text{d}V\\
  
Derivation Using the Divergence Theorem  

.. math::
  \int_{\mathrm{CS}} ((\rho v_{1}\vec{V}) \cdot \vec{n}) d A 
  =\int_{\mathrm{CV}}  \nabla\cdot(\rho v_{1}\vec{V}) d V
  
-

.. math::
  \int_{\mathrm{CS}} ((\rho v_{2}\vec{V}) \cdot \vec{n}) d A 
  =\int_{\mathrm{CV}}  \nabla\cdot(\rho v_{2}\vec{V}) d V
  
-

.. math::
  \int_{\mathrm{CS}} ((\rho v_{3}\vec{V}) \cdot \vec{n}) d A 
  =\int_{\mathrm{CV}}  \nabla\cdot(\rho v_{3}\vec{V}) d V  
  
-
  
.. math::
  \int_{\mathrm{CS}} (\rho v_{i}v_{j}n_{j}) d A 
  =\int_{\mathrm{CV}}  \cfrac{\partial (\rho v_{i}v_{j})}{\partial x_{j}} d V
  
-
  
.. math::
  \int_{\mathrm{CS}} (\rho \vec{V}\otimes\vec{V})\cdot\vec{n} d A 
  =\int_{\mathrm{CV}}\text{div} (\rho \vec{V}\otimes\vec{V}) d V  
  
-
  
.. math::
  \text{div }\mathbf{T}=\nabla\cdot [\mathbf{T}^{\text{T}}]  
  
Let

.. math::
  \mathbf{T}=\rho \vec{V}\otimes\vec{V}=\rho \vec{V}\vec{V}=\mathbf{T}^{\text{T}}
  
The Dyad (the tensor product)

.. math::  
  (\mathbf{a}\otimes\mathbf{b})\cdot\mathbf{c}=(\mathbf{a}\mathbf{b})\cdot\mathbf{c}=\mathbf{a}(\mathbf{b}\cdot\mathbf{c})  
  
-
  
.. math:: 
  (\rho \vec{V}\otimes\vec{V})\cdot\vec{n}=(\rho \vec{V}\vec{V})\cdot\vec{n}=(\rho \vec{V})(\vec{V}\cdot\vec{n})  
  
-
  
.. math::  
  \int_{\mathrm{CS}} \boldsymbol\sigma\cdot\vec{n} d A 
  =\int_{\mathrm{CV}}\text{div} \boldsymbol\sigma d V\\  
  
-
  
.. math:: 
  \int_{\text{CV}}\left[\cfrac{\partial \rho}{\partial t}+\text{div}(\rho\vec{V}\otimes\vec{V})-\rho\vec{g}-\text{div}\boldsymbol\sigma\right]dV=0  
  
Hence, we have a general differential equation
for linear momentum, known as Cauchy’s equation,

Cauchy’s equation:

.. math:: 
  \cfrac{\partial \rho}{\partial t}+\text{div}(\rho\vec{V}\otimes\vec{V})=\rho\vec{g}+\text{div}\boldsymbol\sigma

-

.. math:: 
  \sum \vec{F}=\frac{d ({m}  \vec{V})_{\mathrm{sys}}}{d t}=\frac{d}{d t} \int_{\mathrm{CV}} \rho \vec{V} d V+\int_{\mathrm{CS}} (\rho\vec{V}) ((\vec{V}-\vec{V}_{\text{CS}}) \cdot \vec{n}) d A=\int_{\mathrm{CV}} \rho \vec{g}  d V+\int_{\mathrm{CS}} \sigma_{ij}\cdot \vec{n} d A
  