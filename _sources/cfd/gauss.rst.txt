Gauss’s Theorem
==================================

In the derivation of the conservation equations, Gauss’s theorem is frequently used. This theorem relates integrals over a domain to an integral over the boundary of this domain. It can be
used to relate a volume integral to a surface integral or an area integral to a contour integral.

Gauss’s theorem states that when :math:`f(\mathbf{x})` is piecewise continuously differentiable, that is, a
:math:`C^{0}` function, then

.. math::
  \int\limits_{\Omega}\cfrac{\partial f(\mathbf{x},t)}{\partial x_{i}}\text{d}\Omega=\int\limits_{\Gamma} n_{i}f(\mathbf{x}) d \Gamma 
  
The Divergence Theorem-vector 
------------------------------------
Consider an arbitrary differentiable vector field :math:`v(\mathbf{x},t)` defined in some finite region of
physical space. Let :math:`V` be a volume in this space with a closed surface S bounding the
volume, and let the outward normal to this bounding surface be :math:`\mathbf{n}`. The divergence
theorem of Gauss states that (in symbolic and index notation)

.. math::
  \int\limits_{S}\mathbf{v}\cdot\mathbf{n}\text{d}S=\int\limits_{V}\text{div }\mathbf{v}\text{d}V

-  

.. math::
  \int\limits_{S}{v}_{i}\cdot{n}_{i}\text{d}S=\int\limits_{V}\cfrac{\partial {v}_{i}}{\partial {x}_{i}}\text{d}V\\
  
The divergence theorem is also called Gauss’s theorem, named after the
German mathematician Johann Carl Friedrich Gauss (1777–1855). The
divergence theorem allows us to transform a volume integral of the divergence of a vector into an area integral over the surface that defines the volume.
For any vector :math:`\vec{G}`, the divergence of :math:`\vec{G}` is defined as :math:`\nabla \cdot \vec{G}`, and the divergence theorm is 
written as 

.. math::
  \int_{V}\nabla\cdot\vec{G}=\oint_{A}\vec{G}\cdot\vec{n}dA
  
-
  
.. math::
  \vec{G}=g_{1}(x,y,z)\vec{i}+g_{2}(x,y,z)\vec{j}+g_{3}(x,y,z)\vec{k}
  =\begin{bmatrix}
   g_{1}(x,y,z)\\g_{2}(x,y,z)\\g_{3}(x,y,z)
  \end{bmatrix}\\ 
  
-
  
.. math::
  \nabla \cdot \vec{G}=\cfrac{\partial g_{1}}{\partial x}+\cfrac{\partial g_{2}}{\partial y}+\cfrac{\partial g_{3}}{\partial z} 
  
-
  
.. math::
  \vec{G}\cdot\vec{n}=g_{1}n_{1}+g_{2}n_{2}+g_{3}n_{3}=g_{1}n_{x}+g_{2}n_{y}+g_{3}n_{z}  
  
-
  
.. math::  
  \vec{n}=n_{1}\vec{i}+n_{2}\vec{j}+n_{3}\vec{k}=n_{x}\vec{i}+n_{y}\vec{j}+n_{z}\vec{k}
  =\begin{bmatrix}
   n_{1}\\n_{2}\\n_{3}
  \end{bmatrix}
  =\begin{bmatrix}
   n_{x}\\n_{y}\\n_{z}
  \end{bmatrix}\\  

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

