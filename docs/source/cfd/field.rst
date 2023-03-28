Vector Calculus
==================================

Vectors and Vector Operations
---------------------------------
.. math:: 
  \mathbf{a}=a_{x}\mathbf{i}+a_{y}\mathbf{j}+a_{z}\mathbf{k}

.. math:: 
  \mathbf{a}=\begin{bmatrix}
   a_{x} \\
   a_{y}\\
   a_{z}
  \end{bmatrix}  
  
.. math:: 
  \mathbf{a}^{T}=\begin{bmatrix}
   a_{x} &
   a_{y} &
   a_{z}
  \end{bmatrix}
  
.. math:: 
  \mathbf{v}=u\mathbf{i}+v\mathbf{j}+w\mathbf{k}

.. math:: 
  \mathbf{v}=\begin{bmatrix}
   u \\
   v\\
   w
  \end{bmatrix}  
  
.. math:: 
  \mathbf{v}^{T}=\begin{bmatrix}
   u &
   v &
   w
  \end{bmatrix}  


The Cross Product of Two Vectors
---------------------------------
It is then clear that the cross product of two collinear vectors is zero as they
define no area, and that the cross product of two orthogonal unit vectors is a unit
vector perpendicular to the two unit vectors. Adopting the right hand rule to define
the direction of the resulting vector, the following cross product operations hold:

.. math:: 
  \begin{align}
    \mathbf{i}\times \mathbf{i}&=\mathbf{j}\times \mathbf{j}=\mathbf{k}\times \mathbf{k}=0 \\
    \mathbf{i}\times \mathbf{j}&=\mathbf{k}=-\mathbf{j}\times \mathbf{i}\\
    \mathbf{j}\times \mathbf{k}&=\mathbf{i}=-\mathbf{k}\times \mathbf{j} \\
    \mathbf{k}\times \mathbf{i}&=\mathbf{j}=-\mathbf{i}\times \mathbf{k} \\
  \end{align}

Using the above relations, the cross product of two vectors in terms of their
Cartesian components is given by

.. math::
  \begin{align}
  \mathbf{a}\times \mathbf{b}&=(a_{x}\mathbf{i}+a_{y}\mathbf{j}+a_{z}\mathbf{k}) \times (b_{x}\mathbf{i}+b_{y}\mathbf{j}+b_{z}\mathbf{k})\\
  &=(a_{x}b_{x}\mathbf{i}\times\mathbf{i}+a_{x}b_{y}\mathbf{i}\times\mathbf{j}+a_{x}b_{z}\mathbf{i}\times\mathbf{k})\\
  &+(a_{y}b_{x}\mathbf{j}\times\mathbf{i}+a_{y}b_{y}\mathbf{j}\times\mathbf{j}+a_{y}b_{z}\mathbf{j}\times\mathbf{k})\\
  &+(a_{z}b_{x}\mathbf{k}\times\mathbf{i}+a_{z}b_{y}\mathbf{k}\times\mathbf{j}+a_{z}b_{z}\mathbf{k}\times\mathbf{k})\\
  &=(a_{x}b_{x}0+a_{x}b_{y}\mathbf{k}+a_{x}b_{z}\mathbf{(-j)})\\
  &+(a_{y}b_{x}\mathbf{(-k)}+a_{y}b_{y}0+a_{y}b_{z}\mathbf{i}\\
  &+(a_{z}b_{x}\mathbf{j}+a_{z}b_{y}\mathbf{(-i)}+a_{z}b_{z}0)\\
  &=(a_{y}b_{z}-a_{z}b_{y})\mathbf{i}
  +(a_{z}b_{x}-a_{x}b_{z})\mathbf{j}
  +(a_{x}b_{y}-a_{y}b_{x})\mathbf{k}  
  \end{align}

which can be written using determinant notation as
  
.. math::
  \begin{align}
  \mathbf{a}\times \mathbf{b}&= 
      \begin{vmatrix}
        \mathbf{i}& \mathbf{j} & \mathbf{k}\\
        a_{x} &  a_{y}&a_{z}\\
        b_{x} &  b_{y}&b_{z}
      \end{vmatrix} \\ 
      &=(a_{y}b_{z}-a_{z}b_{y})\mathbf{i}
       +(a_{z}b_{x}-a_{x}b_{z})\mathbf{j}
       +(a_{x}b_{y}-a_{y}b_{x})\mathbf{k}
  \end{align}
  
Gradient
-------------

.. math:: 
  grad\varphi =\frac{\partial \varphi}{\partial x}\mathbf{i}+\frac{\partial \varphi}{\partial y}\mathbf{j}+\frac{\partial \varphi}{\partial z}\mathbf{k}
  
Divergence
-------------

.. math:: 
  div \mathbf{a} =\frac{\partial a_{x} }{\partial x}+\frac{\partial a_{y}}{\partial y}+\frac{\partial a_{z}}{\partial z}
  
Gradient of a Scalar and Directional Derivatives  
--------------------------------------------------

An important vector operator, which arises frequently in fluid dynamics, is the “del”
(or “nabla”) operator defined as

.. math:: 
  \nabla=\frac{\partial}{\partial x} \mathbf{i}+\frac{\partial}{\partial y} \mathbf{j}+\frac{\partial}{\partial z} \mathbf{k}
  
When the “del” operator is applied on a scalar variable s it results in the gradient of :math:`\varphi` given by

.. math:: 
  \nabla\varphi=\frac{\partial\varphi}{\partial x} \mathbf{i}+\frac{\partial\varphi}{\partial y} \mathbf{j}+\frac{\partial\varphi}{\partial z} \mathbf{k}
  
Operations on the Nabla Operator
--------------------------------------------------  

.. math:: 
  \nabla \cdot \mathbf{v}=\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}
  
.. math:: 
  \nabla^{2} =\nabla \cdot \nabla=\frac{\partial^{2} }{\partial x^{2}}+\frac{\partial^{2}}{\partial y^{2}}+\frac{\partial^{2}}{\partial z^{2}}
  
.. math:: 
  \nabla^{2}\varphi  =\nabla \cdot \nabla{\varphi}=\frac{\partial^{2}{\varphi} }{\partial x^{2}}+\frac{\partial^{2}{\varphi}}{\partial y^{2}}+\frac{\partial^{2}{\varphi}}{\partial z^{2}}  
  
.. math::
  \nabla^{2}\mathbf{v}=\nabla^{2}(u\mathbf{i}+v\mathbf{j}+w\mathbf{k}) =(\nabla^{2}u)\mathbf{i}+(\nabla^{2}v)\mathbf{j}+(\nabla^{2}w)\mathbf{k}  
  
.. math::
  \begin{align}
  \nabla\times \mathbf{a}&=(\frac{\partial}{\partial x} \mathbf{i}+\frac{\partial}{\partial y} \mathbf{j}+\frac{\partial}{\partial z} \mathbf{k})
   \times (a_{x} \mathbf{i}+a_{y}\mathbf{j}+a_{z}\mathbf{k}) \\
  &=\begin{vmatrix}
    \mathbf{i}& \mathbf{j} & \mathbf{k}\\
    \frac{\partial}{\partial x}& \frac{\partial}{\partial y} & \frac{\partial}{\partial z}\\
    a_{x} &  a_{y}&a_{z}
  \end{vmatrix}\\  
  &=\left ( \frac{\partial a_{z} }{\partial y} -\frac{\partial a_{y}}{\partial z}\right ) \mathbf{i}+
  \left ( \frac{\partial a_{x} }{\partial y} -\frac{\partial a_{z}}{\partial z}\right ) \mathbf{j}+
  \left ( \frac{\partial a_{y} }{\partial y} -\frac{\partial a_{x}}{\partial z}\right ) \mathbf{k}  
  \end{align}  
  

.. math::  
  \begin{aligned}
  {[(\mathbf{v} \cdot \nabla) \mathbf{v}] } & =(u \mathbf{i}+v \mathbf{j}+w \mathbf{k}) \cdot\left(\frac{\partial}{\partial x} \mathbf{i}+\frac{\partial}{\partial y} \mathbf{j}+\frac{\partial}{\partial z} \mathbf{k}\right)(u \mathbf{i}+v \mathbf{j}+w \mathbf{k}) \\
  & =\left(u \frac{\partial}{\partial x}+v \frac{\partial}{\partial y}+w \frac{\partial}{\partial z}\right)(u \mathbf{i}+v \mathbf{j}+w \mathbf{k}) \\
  & =\left(u \frac{\partial u}{\partial x}+v \frac{\partial u}{\partial y}+w \frac{\partial u}{\partial z}\right) \mathbf{i} \\
  & +\left(u \frac{\partial v}{\partial x}+v \frac{\partial v}{\partial y}+w \frac{\partial v}{\partial z}\right) \mathbf{j} \\
  & +\left(u \frac{\partial w}{\partial x}+v \frac{\partial w}{\partial y}+w \frac{\partial w}{\partial z}\right) \mathbf{k}
  \end{aligned}  
  
Additional Vector Operations  
--------------------------------------------------

.. math::
  \nabla \times ( \mathbf{v_{1}} \times \mathbf{v_{2}}) =(\nabla\cdot\mathbf{v_{2}}+ \mathbf{v_{2}}\cdot\nabla)\mathbf{v_{1}}-(\nabla\cdot\mathbf{v_{1}}+\mathbf{v_{1}}\cdot\nabla)\mathbf{v_{2}}
  
.. math::
  \nabla \times ( \mathbf{v} \times \mathbf{F}) =(\nabla\cdot\mathbf{F}+ \mathbf{F}\cdot\nabla)\mathbf{v}-(\nabla\cdot\mathbf{v}+\mathbf{v}\cdot\nabla)\mathbf{F}

.. math::
  \nabla \cdot(\phi  \mathbf{v})=\phi \nabla \cdot \mathbf{v}+\mathbf{v} \cdot \nabla \phi



  
  