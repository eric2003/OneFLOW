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

