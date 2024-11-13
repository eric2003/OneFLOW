Strain Measures
==================================

Deformation Gradient
------------------------------------
The description of deformation and the measure of strain are essential parts of nonlinear
continuum mechanics. An important variable in the characterization of deformation is the
deformation gradient. The deformation gradient is defined by

.. math::
  F_{ij}=\cfrac{\partial \phi_{i}}{\partial X_{j}}=\cfrac{\partial x_{i}}{\partial X_{j}}

or 

.. math::
  \mathbf{F}=\cfrac{\partial \boldsymbol\phi}{\partial \mathbf{X}}=\cfrac{\partial \mathbf{x}}{\partial \mathbf{X}}
  =(\nabla_{0}\boldsymbol\phi)^{\text{T}}
  
In the terminology of mathematics, the deformation gradient :math:`F` is the Jacobian matrix of the
motion :math:`\boldsymbol\phi(\mathbf{X}, t)`. Note in the above that the first index of :math:`F_{ij}` refers to the motion, the second to the
partial derivative. The operator :math:`\nabla_{0}` is the left gradient with respect to the material coordinates.  

If we consider an infinitesimal line segment :math:`d\mathbf{X}` in the reference configuration, then the corresponding line segment :math:`d\mathbf{x}` in the current configuration is
given by

.. math::
  d\mathbf{x}=\mathbf{F}\cdot d\mathbf{X}
  
or   

.. math::
  d{x}_{i}={F}_{ij}\cdot d{X}_{j}
  
In the above expression, the dot could have been omitted between the :math:`\mathbf{F}` and :math:`d\mathbf{X}`, since the
expression is also valid as a matrix expression. We have retained it to conform to our convention
of always explicitly indicating contractions in tensor expressions.  

In two dimensions, the deformation gradient in a rectangular coordinate system is given by

.. math::
  \mathbf{F}=\left[\begin{array}{ll}
  \cfrac{\partial x_{1}}{\partial X_{1}} & \cfrac{\partial x_{1}}{\partial X_{2}} \\
  \cfrac{\partial x_{2}}{\partial X_{1}} & \cfrac{\partial x_{2}}{\partial X_{2}}
  \end{array}\right]=\left[\begin{array}{ll}
  \cfrac{\partial x}{\partial X} & \cfrac{\partial x}{\partial Y} \\
  \cfrac{\partial y}{\partial X} & \cfrac{\partial y}{\partial Y}
  \end{array}\right]
  
As can be seen in the above, in writing a second-order tensor in matrix form, we use the first
index for the row number, and the second index for the column number. Note that :math:`\mathbf{F}` is the
transpose of the left-gradient.  
  
The determinant of :math:`\mathbf{F}` is denoted by :math:`J` and called the Jacobian determinant or the determinant of the deformation gradient

.. math::
  J=\text{det}(\mathbf{F})
  
The Jacobian determinant can be used to relate integrals in the current and reference configurations by  

.. math::
  \int\limits_{\Omega}f(\mathbf{x},t)\text{d}\Omega=\int\limits_{\Omega}f(\boldsymbol\phi(\mathbf{X},t),t)J\text{d}\Omega_{0}
  
or 

.. math::
  \int\limits_{\Omega}f\text{d}\Omega=\int\limits_{\Omega}fJ\text{d}\Omega_{0}
  
or in two dimensions

.. math::
  \int\limits_{\Omega}f(x,y)\text{d}x\text{d}y=\int\limits_{\Omega}f(X,Y)J\text{d}X\text{d}Y
  
The material derivative of the Jacobian determinant is given by

.. math::
  \cfrac{\text{d}J}{\text{d}t}=\dot{J} =J\text{ div }\mathbf{v}=J\cfrac{\partial{v}_{i}}{\partial{x}_{i}}  
  
  
Strain Measures
----------------------------

In contrast to linear elasticity, many different measures of strain and strain rate are used in
nonlinear continuum mechanics. Only two of these measures are considered here:

1. The Green (Green–Lagrange) strain, E.
2. The rate-of-deformation tensor, D

Green Strain Tensor
----------------------------

The Green strain tensor :math:`E` is defined by

.. math::
  ds^{2}-dS^{2}=2d\mathbf{X}\cdot\mathbf{E}\cdot d\mathbf{X}
  
or 

.. math::
  dx_{i}dx_{i}-dX_{i}dX_{i}=2dX_{i}E_{ij}dX_{j}
  
Rate-of-Deformation  
----------------------------
The second kinematic measure to be considered here is the rate-of-deformation :math:`\mathbf{D}`. It is also
called the velocity strain. In contrast to the Green strain tensor, it is a rate measure of
deformation.

To develop an expression for the rate-of-deformation, we first define the velocity
gradient :math:`\mathbf{L}` by

.. math::
  \mathbf{L}=\cfrac{\partial \mathbf{v}}{\partial \mathbf{x}}=\text{grad }\mathbf{v}=(\nabla \mathbf{v})^{\text{T}}\\
  
or   

.. math::
  L_{ij} = \cfrac{\partial v_{i}}{\partial x_{j}}
  
-
  
.. math::
  \text{d}\mathbf{v}=\mathbf{L}\cdot \text{d}\mathbf{x}  
  
-
  
.. math::
  \text{d}{v}_{i}=L_{ij}\text{d}{x}_{j}
  
We have shown several tensor forms of the definition, but we will primarily use the indicial
form. The symbol :math:`\nabla` preceding the function denotes the left spatial gradient of the function: in a spatial gradient, the derivatives are taken with respect
to the spatial (Eulerian) coordinates. The symbol :math:`\nabla` specifies the spatial gradient; :math:`\nabla_{0}` is the
material gradient.  

The velocity gradient tensor can be decomposed into symmetric and skew-symmetric
parts by

.. math::
  \mathbf{L}=\cfrac{1}{2}(\mathbf{L}+\mathbf{L}^{\text{T}})+\cfrac{1}{2}(\mathbf{L}-\mathbf{L}^{\text{T}})
  
or

.. math::  
  L_{ij}=\cfrac{1}{2}({L}_{ij}+{L}_{ji})+\cfrac{1}{2}({L}_{ij}-{L}_{ji})

This is a standard decomposition of a second-order tensor or square matrix: any secondorder tensor can be expressed as the sum of its symmetric and skew-symmetric parts in this
manner.

The rate-of-deformation :math:`\mathbf{D}` is defined as the symmetric part of :math:`\mathbf{L}` and the spin :math:`\mathbf{W}` is the skew-symmetric part of :math:`\mathbf{L}`. Using these definitions, we can write  

.. math::  
  \mathbf{D}=\cfrac{1}{2}(\mathbf{L}+\mathbf{L}^{\text{T}}) \quad or \quad 
  D_{ij}=\cfrac{1}{2}\left(\cfrac{\partial v_{i}}{\partial x_{j}}+\cfrac{\partial v_{j}}{\partial x_{i}}\right)\\
  
-
  
.. math::  
  \mathbf{W}=\cfrac{1}{2}(\mathbf{L}-\mathbf{L}^{\text{T}}) \quad or \quad 
  W_{ij}=\cfrac{1}{2}\left(\cfrac{\partial v_{i}}{\partial x_{j}}-\cfrac{\partial v_{j}}{\partial x_{i}}\right)

-
  
.. math::  
  \mathbf{L}=(\nabla\mathbf{v})^{\text{T}}=\mathbf{D}+\mathbf{W} \quad or \quad {L}_{ij}=\cfrac{\partial v_{i}}{\partial x_{j}}={D}_{ij}+{W}_{ij}\\  