Vector & Tensor
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
  
-  
  
.. math:: 
  \nabla^{2} =\nabla \cdot \nabla=\frac{\partial^{2} }{\partial x^{2}}+\frac{\partial^{2}}{\partial y^{2}}+\frac{\partial^{2}}{\partial z^{2}}

-
  
.. math:: 
  \nabla^{2}\varphi  =\nabla \cdot \nabla{\varphi}=\frac{\partial^{2}{\varphi} }{\partial x^{2}}+\frac{\partial^{2}{\varphi}}{\partial y^{2}}+\frac{\partial^{2}{\varphi}}{\partial z^{2}}  

-
  
.. math::
  \nabla^{2}\mathbf{v}=\nabla^{2}(u\mathbf{i}+v\mathbf{j}+w\mathbf{k}) =(\nabla^{2}u)\mathbf{i}+(\nabla^{2}v)\mathbf{j}+(\nabla^{2}w)\mathbf{k}  

-
  
.. math::
  \begin{align}
  \nabla\times \mathbf{a}&=(\cfrac{\partial}{\partial x} \mathbf{i}+\cfrac{\partial}{\partial y} \mathbf{j}+\cfrac{\partial}{\partial z} \mathbf{k})
   \times (a_{x} \mathbf{i}+a_{y}\mathbf{j}+a_{z}\mathbf{k}) \\
  &=\begin{vmatrix}
    \mathbf{i}& \mathbf{j} & \mathbf{k}\\
    \cfrac{\partial}{\partial x}& \cfrac{\partial}{\partial y} & \cfrac{\partial}{\partial z}\\
    a_{x} &  a_{y}&a_{z}
  \end{vmatrix}\\  
  &=\left ( \cfrac{\partial a_{z} }{\partial y} -\cfrac{\partial a_{y}}{\partial z}\right ) \mathbf{i}+
    \left ( \cfrac{\partial a_{x} }{\partial z} -\cfrac{\partial a_{z}}{\partial x}\right ) \mathbf{j}+
    \left ( \cfrac{\partial a_{y} }{\partial x} -\cfrac{\partial a_{x}}{\partial y}\right ) \mathbf{k}  
  \end{align} 

-
  
.. math::
  \begin{align}
  \nabla\times \mathbf{v}&=
  \begin{vmatrix}
    \mathbf{i}& \mathbf{j} & \mathbf{k}\\
    \cfrac{\partial}{\partial x}& \cfrac{\partial}{\partial y} & \cfrac{\partial}{\partial z}\\
    u &  v&w
  \end{vmatrix}\\
  &=\left ( \cfrac{\partial w }{\partial y} -\cfrac{\partial v}{\partial z}\right ) \mathbf{i}+
    \left ( \cfrac{\partial u }{\partial z} -\cfrac{\partial w}{\partial x}\right ) \mathbf{j}+
    \left ( \cfrac{\partial v }{\partial x} -\cfrac{\partial u}{\partial y}\right ) \mathbf{k}  
  \end{align}      

-
  
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
  
The Dyad (the tensor product)
``````````````````````````````````````
The vector dot product and vector cross product have been considered in previous
sections. A third vector product, the tensor product (or dyadic product), is important in
the analysis of tensors of order 2 or more. The tensor product of two vectors u and v is
written as

.. math::
   \mathbf{u} \otimes \mathbf{v}
   
This tensor product is itself a tensor of order two, and is called dyad:

.. math::
  \begin{align}
   &\mathbf{u} \cdot \mathbf{v}   &\text{ is a scalar}\, \quad  &\text{(a zeroth order tensor)}\\
   &\mathbf{u} \times \mathbf{v}  &\text{ is a vector} \quad  &\text{(a first order tensor)}\quad\\
   &\mathbf{u} \otimes \mathbf{v} &\text{ is a dyad }\,\, \quad  &\text{(a second order tensor)}
  \end{align} 

It is best to define this dyad by what it does: it transforms a vector :math:`\mathbf{w}` into another vector
with the direction of :math:`\mathbf{u}` according to the rule

.. math::
  (\mathbf{u} \otimes \mathbf{v})\mathbf{w}=\mathbf{u}(\mathbf{v}\cdot\mathbf{w}) \quad \text{The Dyad Transformation}
  
The following important relations follow from the above definition

.. math::
  ((\mathbf{u} \otimes \mathbf{v})\mathbf{w})\otimes\mathbf{x}=\mathbf{u}(\mathbf{v}\cdot\mathbf{w})\otimes\mathbf{x}=(\mathbf{v}\cdot\mathbf{w})(\mathbf{u}\otimes\mathbf{x})\\

.. math::
  \begin{align}
  \mathbf{a}=a_{i}\mathbf{e}_{i}=&a_{1}\mathbf{e}_{1}+a_{2}\mathbf{e}_{2}+\cdots +a_{n}\mathbf{e}_{n}\\
  \mathbf{b}=b_{i}\mathbf{e}_{i}=&b_{1}\mathbf{e}_{1}+b_{2}\mathbf{e}_{2}+\cdots +b_{n}\mathbf{e}_{n}\\
  \end{align}
  
-

.. math::
  \mathbf{a}\mathbf{b}=\mathbf{a}\otimes \mathbf{b}=a_{i}b_{j}\mathbf{e}_{i}\otimes\mathbf{e}_{j}

-

.. math::
  \mathbf{a}\mathbf{b}=\mathbf{a}\otimes \mathbf{b}=(\mathbf{a})(\mathbf{b}^{\text{T}})
  =\begin{bmatrix}
   a_{1}\\
   a_{2}\\
   \vdots\\
   a_{n}\\
  \end{bmatrix}
  \begin{bmatrix}
   b_{1}&b_{2}&\cdots &b_{n}\\
  \end{bmatrix}
  =\begin{bmatrix}
   a_{1}b_{1}&a_{1}b_{2}&\cdots &a_{1}b_{n}\\
   a_{2}b_{1}&a_{2}b_{2}&\cdots &a_{2}b_{n}\\
   \vdots &\vdots&\ddots &\vdots\\
   a_{n}b_{1}&a_{n}b_{2}&\cdots &a_{n}b_{n}\\
  \end{bmatrix}
  
The following identities are a direct consequence of the definition of the tensor product:

- Compatible with scalar multiplication:

.. math::
  (\alpha\mathbf{a})\otimes\mathbf{b}=\mathbf{a}\otimes(\alpha\mathbf{b})=\alpha(\mathbf{a}\otimes\mathbf{b})
  
- Distributive over vector addition:

.. math::
  \begin{align}
  \mathbf{a}\otimes(\mathbf{b}+\mathbf{c})=\mathbf{a}\otimes\mathbf{b}+\mathbf{a}\otimes\mathbf{c}\\
  (\mathbf{a}+\mathbf{b})\otimes\mathbf{c}=\mathbf{a}\otimes\mathbf{c}+\mathbf{b}\otimes\mathbf{c}\\
  \end{align}
  
Product of dyadic and vector 
There are four operations defined on a vector and dyadic, constructed from the products defined on vectors. 

Dot product

.. math::
  \begin{align}
  \mathbf{c}\cdot(\mathbf{a}\otimes\mathbf{b})=(\mathbf{c}\cdot\mathbf{a})\mathbf{b}\\
  (\mathbf{a}\otimes\mathbf{b})\cdot\mathbf{c}=\mathbf{a}(\mathbf{b}\cdot\mathbf{c})\\
  \end{align}


Cross product

.. math::
  \begin{align}
  \mathbf{c}\times(\mathbf{a}\otimes\mathbf{b})=(\mathbf{c}\times\mathbf{a})\otimes\mathbf{b}\\
  (\mathbf{a}\otimes\mathbf{b})\times\mathbf{c}=\mathbf{a}\otimes(\mathbf{b}\times\mathbf{c})\\
  \end{align}

Tensors and Tensor Operation 
--------------------------------------------------
Tensors can be thought of as extensions to the ideas already used when defining
quantities like scalars and vectors. A scalar is a tensor of rank zero, and a
vector is a tensor of rank one. Tensors of higher rank (2, 3, etc.) can be developed
and their main use is to manipulate and transform sets of equations. 

Let :math:`x`, :math:`y` and :math:`z` represent the directions in an orthonormal Cartesian coordinate
system, then the stress tensor :math:`\tau` and its transpose designated with superscript :math:`T({\tau}^{T})`
are represented in terms of their components as

.. math::
  \tau=\begin{bmatrix}
  \tau_{xx}& \tau_{xy} & \tau_{xz}\\
  \tau_{yx}& \tau_{yy} & \tau_{yz}\\
  \tau_{zx}& \tau_{zy} & \tau_{zz}\\
  \end{bmatrix}

-
  
.. math::
  \tau^{T}=\begin{bmatrix}
  \tau_{xx}& \tau_{yx} & \tau_{zx}\\
  \tau_{xy}& \tau_{yy} & \tau_{xy}\\
  \tau_{xz}& \tau_{yz} & \tau_{zz}\\
  \end{bmatrix}  
  
Similar to writing a vector in terms of its components, defining the unit vectors :math:`\mathbf{i}`, :math:`\mathbf{j}`,
and :math:`\mathbf{k}` in the :math:`x`, :math:`y` and :math:`z` direction, respectively, the tensor :math:`\tau` given can be
written in terms of its components as 
 
.. math::
  \begin{align}
  \boldsymbol{\tau}
  & = \mathbf{ii}\tau_{xx}+\mathbf{ij}\tau_{xy}+\mathbf{ik}\tau_{xz}\\
     & + \mathbf{ji}\tau_{yx}+\mathbf{jj}\tau_{yy}+\mathbf{jk}\tau_{yz}\\
     & + \mathbf{ki}\tau_{zx}+\mathbf{kj}\tau_{zy}+\mathbf{kk}\tau_{zz}\\
  \end{align}
  
Second Order Tensors
``````````````````````````````````  
A second-order tensor :math:`\mathbf{T}` may be defined as an operator that acts on a vector :math:`\mathbf{u}` generating
another vector :math:`\mathbf{v}`, so that :math:`\mathbf{T}(\mathbf{u})=\mathbf{v}` or

.. math::
  \mathbf{T}(\mathbf{u})  = \mathbf{v} \quad or \quad  \mathbf{T}\cdot {\mathbf{u}} = \mathbf{v} \quad or \quad \mathbf{T}{\mathbf{u}}  = \mathbf{v}
  
  
Second Order Tensor as a Dyadic
``````````````````````````````````
In what follows, it will be shown that a second order tensor can always be written as a
dyadic involving the Cartesian base vectors :math:`\mathbf{e}_{i}`

Consider an arbitrary second-order tensor :math:`\mathbf{T}` which operates on :math:`\mathbf{a}` to produce :math:`\mathbf{b}`, :math:`\mathbf{T}(\mathbf{a})=\mathbf{b}`,
or :math:`\mathbf{T}(a_{i}\mathbf{e}_{i})=\mathbf{b}`. From the linearity of :math:`\mathbf{T}`, 

.. math::
  a_{1}\mathbf{T}(\mathbf{e}_{1})+a_{2}\mathbf{T}(\mathbf{e}_{2})+a_{3}\mathbf{T}(\mathbf{e}_{3})=\mathbf{b}
  
Just as :math:`\mathbf{T}` transforms :math:`\mathbf{a}` into :math:`\mathbf{b}`, it transforms the base vectors :math:`\mathbf{e}_{i}` into some other vectors;
suppose that :math:`\mathbf{T}(\mathbf{e}_{1})=\mathbf{u},\mathbf{T}(\mathbf{e}_{2})=\mathbf{v},\mathbf{T}(\mathbf{e}_{3})=\mathbf{w}`, then  

.. math::
  \begin{aligned}
  \mathbf{b} & =a_{1} \mathbf{u}+a_{2} \mathbf{v}+a_{3} \mathbf{w} \\
  & =\left(\mathbf{a} \cdot \mathbf{e}_{1}\right) \mathbf{u}+\left(\mathbf{a} \cdot \mathbf{e}_{2}\right) \mathbf{v}+\left(\mathbf{a} \cdot \mathbf{e}_{3}\right) \mathbf{w} \\
  & =\left(\mathbf{u} \otimes \mathbf{e}_{1}\right)\cdot  \mathbf{a}+\left(\mathbf{v} \otimes \mathbf{e}_{2}\right)\cdot  \mathbf{a}+\left(\mathbf{w} \otimes \mathbf{e}_{3}\right)\cdot  \mathbf{a} \\
  & =\left[\mathbf{u} \otimes \mathbf{e}_{1}+\mathbf{v} \otimes \mathbf{e}_{2}+\mathbf{w} \otimes \mathbf{e}_{3}\right] \cdot \mathbf{a}
  \end{aligned}
  
and so

.. math::
  \begin{align}
  &\mathbf{T}(\mathbf{a})=\mathbf{T}\cdot\mathbf{a}=\mathbf{b}=(\mathbf{u} \otimes \mathbf{e}_{1}+\mathbf{v} \otimes \mathbf{e}_{2}+\mathbf{w} \otimes \mathbf{e}_{3})\cdot \mathbf{a}\\
  &\Longrightarrow \mathbf{T}=\mathbf{u} \otimes \mathbf{e}_{1}+\mathbf{v} \otimes \mathbf{e}_{2}+\mathbf{w} \otimes \mathbf{e}_{3}
  \end{align}
  
which is indeed a dyadic.

.. math::
  \mathbf{u}=u_{i}\mathbf{e}_{i}\quad \mathbf{v}=v_{i}\mathbf{e}_{i}\quad \mathbf{w}=w_{i}\mathbf{e}_{i}

-
  
.. math::  
  \begin{align}
  \mathbf{T}&=\mathbf{u} \otimes \mathbf{e}_{1}+\mathbf{v} \otimes \mathbf{e}_{2}+\mathbf{w} \otimes \mathbf{e}_{3}\\
  &=(u_{i}\mathbf{e}_{i})\otimes \mathbf{e}_{1}+(v_{i}\mathbf{e}_{i})\otimes \mathbf{e}_{2}+(w_{i}\mathbf{e}_{i})\otimes \mathbf{e}_{3}\\
  &=(T_{i1}\mathbf{e}_{i})\otimes \mathbf{e}_{1}+(T_{i2}\mathbf{e}_{i})\otimes \mathbf{e}_{2}+(T_{i3}\mathbf{e}_{i})\otimes \mathbf{e}_{3}\\
  &=(T_{ij}\mathbf{e}_{i})\otimes \mathbf{e}_{j}\\
  &=T_{ij}(\mathbf{e}_{i}\otimes \mathbf{e}_{j})  
  \end{align}
  
Cartesian components of a Second Order Tensor 

.. math::
  \mathbf{T}\cdot\mathbf{e}_{k} =T_{ij}(\mathbf{e}_{i}\otimes \mathbf{e}_{j})\cdot\mathbf{e}_{k}=
  T_{ij}\mathbf{e}_{i}(\mathbf{e}_{j}\cdot\mathbf{e}_{k})=T_{ij}\mathbf{e}_{i}\delta_{jk} =T_{ik}\mathbf{e}_{i}\\
  
-
  
.. math::
  \mathbf{e}_{m}\cdot(\mathbf{T}\cdot\mathbf{e}_{k}) =\mathbf{e}_{m}\cdot(T_{ik}\mathbf{e}_{i})
  =T_{ik}(\mathbf{e}_{m}\cdot\mathbf{e}_{i})=T_{ik}\delta_{mi}=T_{mk}\\

-
  
.. math::
  \mathbf{e}_{i}\cdot(\mathbf{T}\cdot\mathbf{e}_{j})=\mathbf{e}_{i}\cdot\mathbf{T}\cdot\mathbf{e}_{j}=T_{ij}\\

Matrix multiplication is associative
Even though matrix multiplication is not commutative, it is associative in the following sense. If :math:`\mathbf{A}`
is an :math:`m\times p` matrix, :math:`\mathbf{B}` is a :math:`p\times q` matrix, and :math:`\mathbf{C}` is a :math:`q\times n` matrix, then

.. math::
  \mathbf{A}(\mathbf{B}\mathbf{C})=(\mathbf{A}\mathbf{B})\mathbf{C}

Tensor Product of Matrices
``````````````````````````````````
#. `Tensor Product of Matrices <https://stemandmusic.in/maths/mvt-algebra/matrixTP.php>`_

- Tensor Product is a Special Type of Kronecker Product carried out between a Co-Vector (Row Matrix) and Some Other Tensor whose Rank is Greater than or Equal to 1.
- Tensor Product Always increases the Rank of the Tensor that the Co-Vector is operating upon.
- Tensor Product of 2 Vectors/Column Matrices is also called Outer Product of Vectors.
- Let's consider 2 Column Matrices (or Vectors), a :math:`m \times 1` Matrix :math:`\mathbf{A}` and a :math:`n \times 1` Matrix :math:`\mathbf{B}` given as follows

.. math::
  \begin{align}
  \mathbf{A}=\begin{bmatrix}
   a_{1} \\
   a_{2} \\
   \vdots \\
   a_{m} \\
  \end{bmatrix}
  \quad 
  \mathbf{B}=\begin{bmatrix}
   a_{1} \\
   a_{2} \\
   \vdots \\
   a_{m} \\
  \end{bmatrix}
  \end{align}
  
The Tensor Product of Vector Matrix :math:`\mathbf{A}` with Vector Matrix :math:`\mathbf{B}` is given as

Dyadic Product
`````````````````
#. `Dyadics <https://en.wikipedia.org/wiki/Dyadics>`_
#. `Dyads and Dyadics <https://stemandmusic.in/maths/mvt-algebra/vectorDyad.php>`_

tensor product

.. math::
  \mathbf{u}=\begin{bmatrix}
   u_{1}\\u_{2}\\u_{3}
  \end{bmatrix}=u_{1}\mathbf{i}+u_{2}\mathbf{j}+u_{3}\mathbf{k}
  
- 
 
.. math::  
  \mathbf{v}=\begin{bmatrix}
  v_{1}\\v_{2}\\v_{3}
  \end{bmatrix}=v_{1}\mathbf{i}+v_{2}\mathbf{j}+v_{3}\mathbf{k}  
  
- 
 
.. math::  
  \mathbf{uv}\equiv \mathbf{u\otimes v}=\mathbf{uv^{T} }=\begin{bmatrix}
  u_{1}\\u_{2}\\u_{3}
  \end{bmatrix}\begin{bmatrix}
   v_{1}& v_{2} &v_{3}
  \end{bmatrix}=\begin{bmatrix}
  u_{1}v_{1} &u_{1}v_{2}  &u_{1}v_{3} \\
  u_{2}v_{1} &u_{2}v_{2}  &u_{2}v_{3} \\
  u_{3}v_{1} &u_{3}v_{2}  &u_{3}v_{3} \\
  \end{bmatrix} 
  
- 

.. math::
  \begin{align}
  (\mathbf{a}\mathbf{b})\cdot \mathbf{c} & = \mathbf{a}(\mathbf{b}\cdot \mathbf{c})\\
  \mathbf{a}\cdot (\mathbf{b}\mathbf{c}) & = (\mathbf{a}\cdot \mathbf{b})\mathbf{c}\\
  \end{align}
  
-
  
.. math::  
  \mathbf{i}=\begin{bmatrix}
   1\\0\\0
  \end{bmatrix}\quad
  \mathbf{j}=\begin{bmatrix}
   0\\1\\0
  \end{bmatrix}\quad
  \mathbf{k}=\begin{bmatrix}
   0\\0\\1
  \end{bmatrix}
  
tensor product

.. math:: 
  \begin{matrix}
  \mathbf{ii}\equiv \mathbf{i\otimes i}=\mathbf{ii^{T} }& \mathbf{ij}\equiv \mathbf{i\otimes j}=\mathbf{ij^{T} } & \mathbf{ik}\equiv \mathbf{i\otimes k}=\mathbf{ik^{T} }\\
  \mathbf{ji}\equiv \mathbf{j\otimes i}=\mathbf{ji^{T} }& \mathbf{jj}\equiv \mathbf{j\otimes j}=\mathbf{jj^{T} } & \mathbf{jk}\equiv \mathbf{j\otimes k}=\mathbf{jk^{T} }\\
  \mathbf{ki}\equiv \mathbf{k\otimes i}=\mathbf{ki^{T} }& \mathbf{kj}\equiv \mathbf{k\otimes j}=\mathbf{kj^{T} } & \mathbf{kk}\equiv \mathbf{k\otimes k}=\mathbf{kk^{T} }\\
  \end{matrix}
    
-

.. math:: 
  \mathbf{ii}=\begin{bmatrix}
   1&0&0\\
   0&0&0\\
   0&0&0
  \end{bmatrix} \quad
  \mathbf{ij}=\begin{bmatrix}
   0&1&0\\
   0&0&0\\
   0&0&0
  \end{bmatrix} \quad
  \mathbf{ik}=\begin{bmatrix}
   0&0&1\\
   0&0&0\\
   0&0&0
  \end{bmatrix} 
  
-
  
.. math::  
  \mathbf{ji}=\begin{bmatrix}
   0&0&0\\
   1&0&0\\
   0&0&0
  \end{bmatrix} \quad
  \mathbf{jj}=\begin{bmatrix}
   0&0&0\\
   0&1&0\\
   0&0&0
  \end{bmatrix} \quad
  \mathbf{jk}=\begin{bmatrix}
   0&0&0\\
   0&0&1\\
   0&0&0
  \end{bmatrix}
  
-
  
.. math:: 
  \mathbf{ki}=\begin{bmatrix}
   0&0&0\\
   0&0&0\\
   1&0&0
  \end{bmatrix} \quad
  \mathbf{kj}=\begin{bmatrix}
   0&0&0\\
   0&0&0\\
   0&1&0
  \end{bmatrix} \quad
  \mathbf{kk}=\begin{bmatrix}
   0&0&0\\
   0&0&0\\
   0&0&1
  \end{bmatrix}   

the dyadic product of a vector :math:`\mathbf{v}` by itself, arising in the formulation of the momentum equation of fluid
flow, gives

.. math::
  \begin{align}
  \mathbf{vv} & = (u\mathbf{i}+v\mathbf{j}+w\mathbf{k})\\
  &=uu\mathbf{ii}+uv\mathbf{ij}+uw\mathbf{ik}\\
  &+vu\mathbf{ji}+vv\mathbf{jj}+vw\mathbf{jk}\\
  &+wu\mathbf{ki}+wv\mathbf{kj}+ww\mathbf{kk}\\
  \end{align}
  
-

.. math::
  \mathbf{vv}=\begin{bmatrix}
  uu& uv & uw\\
  vu& vv & vw\\
  wu& wv & ww\\
  \end{bmatrix}  
  
-

.. math::  
  \mathbf{vv}=\begin{bmatrix}
  v_{x}v_{x}& v_{x}v_{y} & v_{x}v_{z}\\
  v_{y}v_{x}& v_{y}v_{y} & v_{y}v_{z}\\
  v_{z}v_{x}& v_{z}v_{y} & v_{z}v_{z}\\
  \end{bmatrix}  
  
-

.. math::  
  \mathbf{vv}=\begin{bmatrix}
    v_{1}v_{1}& v_{1}v_{2} & v_{1}v_{3}\\
    v_{2}v_{1}& v_{2}v_{2} & v_{2}v_{3}\\
    v_{3}v_{1}& v_{3}v_{2} & v_{3}v_{3}\\
  \end{bmatrix} 
  
The Gradient of a Vector Functions
`````````````````````````````````````` 

The gradient of a vector field is defined to be the second-order tensor

.. math::  
  \text {grad }\mathbf{a}=\cfrac{\partial \mathbf{a}}{\partial x_{j}}\otimes \mathbf{e}_{j}
  =\cfrac{\partial {a}_{i}}{\partial x_{j}}\mathbf{e}_{i}\otimes \mathbf{e}_{j} 
  
In matrix notation,

.. math:: 
  \text {grad }\mathbf{a}=\left[\begin{array}{lll}
  \cfrac{\partial a_{1}}{\partial x_{1}} & \cfrac{\partial a_{1}}{\partial x_{2}} & \cfrac{\partial a_{1}}{\partial x_{3}} \\
  \cfrac{\partial a_{2}}{\partial x_{1}} & \cfrac{\partial a_{2}}{\partial x_{2}} & \cfrac{\partial a_{2}}{\partial x_{3}} \\
  \cfrac{\partial a_{3}}{\partial x_{1}} & \cfrac{\partial a_{3}}{\partial x_{2}} & \cfrac{\partial a_{3}}{\partial x_{3}}
  \end{array}\right]
  
Although for a scalar field :math:`\text{grad}\phi` is equivalent to :math:`\nabla \phi`, note that the gradient  of a vector field 
is not the same as :math:`\nabla \otimes \mathbf{a}`. In fact 

.. math::
  (\nabla \otimes \mathbf{a})^{\text{T}}=\text{grad }\mathbf{a} 
  
since

.. math::
  \nabla \otimes \mathbf{a}=\mathbf{e}_{i} \frac{\partial}{\partial x_{i}} \otimes a_{j} \mathbf{e}_{j}=\frac{\partial a_{j}}{\partial x_{i}} \mathbf{e}_{i} \otimes \mathbf{e}_{j}
  
These two different definitions of the gradient of a vector :math:`({\partial a_{i}}/{\partial x_{j}}) \mathbf{e}_{i} \otimes \mathbf{e}_{j}` 
and :math:`({\partial a_{j}}/{\partial x_{i}}) \mathbf{e}_{i} \otimes \mathbf{e}_{j}` , are both commonly used. In what follows, they will be distinguished by
labeling the former as :math:`\text{grad }\mathbf{a}` (which will be called the gradient of :math:`\mathbf{a}` ) and the latter as :math:`\nabla \otimes \mathbf{a}`.

Note the following:

- in much of the literature, :math:`\nabla \otimes \mathbf{a}` is written in the contracted form :math:`\nabla \mathbf{a}`, but the more explicit version is used here.
- some authors define the operation of through :math:`\nabla \otimes(\cdot)\equiv ({\partial (\cdot)}/{\partial x_{i}})\otimes \mathbf{e}_{i}` so that :math:`\nabla \otimes \mathbf{a}=\text{grad}\mathbf{a}=({\partial a_{i}}/{\partial x_{j}}) \mathbf{e}_{i} \otimes \mathbf{e}_{j}`

.. math::
  \begin{align}
  \nabla \otimes(\mathbf{a})&=({\partial (\mathbf{a})}/{\partial x_{i}})\otimes \mathbf{e}_{i}\\
  &=({\partial (a_{j}\mathbf{e}_{j})}/{\partial x_{i}})\otimes \mathbf{e}_{i}\\
  &=({\partial (a_{j})}/{\partial x_{i}})\mathbf{e}_{j}\otimes \mathbf{e}_{i}\\
  &=({\partial a_{i}}/{\partial x_{j}})\mathbf{e}_{i}\otimes \mathbf{e}_{j}\\
  \end{align}
  
The :math:`\nabla {\mathbf{v}}` is a tensor given by

.. math:: 
  \begin{align}
  \nabla {\mathbf{v}}&=\nabla \otimes{\mathbf{v}}=(\nabla) ({\mathbf{v}})^{\text T}\\
  &=\begin{bmatrix}
  \cfrac{\partial }{\partial x}\\
  \cfrac{\partial }{\partial y}\\
  \cfrac{\partial }{\partial z}\\
  \end{bmatrix}
  \begin{bmatrix}
  u&v&w
  \end{bmatrix}=
  \begin{bmatrix}
  \cfrac{\partial u}{\partial x}& \cfrac{\partial v}{\partial x} & \cfrac{\partial w}{\partial x}\\
  \cfrac{\partial u}{\partial y}& \cfrac{\partial v}{\partial y} & \cfrac{\partial w}{\partial y}\\
  \cfrac{\partial u}{\partial z}& \cfrac{\partial v}{\partial z} & \cfrac{\partial w}{\partial z}\\
  \end{bmatrix}
  \end{align}
  
-

.. math::  
  \begin{align}
  \nabla \mathbf{v} & = (\cfrac {\partial }{\partial x}\mathbf{i}+\cfrac {\partial }{\partial y}\mathbf{j}+\cfrac {\partial }{\partial z}\mathbf{k})(u\mathbf{i}+v\mathbf{j}+w\mathbf{k})\\
  &=\cfrac{\partial u}{\partial x}\mathbf{ii}+\cfrac{\partial v}{\partial x}\mathbf{ij}+\cfrac{\partial w}{\partial x}\mathbf{ik}\\
  &+\cfrac{\partial u}{\partial y}\mathbf{ji}+\cfrac{\partial v}{\partial y}\mathbf{jj}+\cfrac{\partial w}{\partial y}\mathbf{jk}\\
  &+\cfrac{\partial u}{\partial z}\mathbf{ki}+\cfrac{\partial v}{\partial z}\mathbf{kj}+\cfrac{\partial w}{\partial z}\mathbf{kk}\\
  \end{align}  
  
-

.. math::   
  \nabla \mathbf{v}=\begin{bmatrix}
  \cfrac{\partial u}{\partial x}& \cfrac{\partial v}{\partial x} & \cfrac{\partial w}{\partial x}\\
  \cfrac{\partial u}{\partial y}& \cfrac{\partial v}{\partial y} & \cfrac{\partial w}{\partial y}\\
  \cfrac{\partial u}{\partial z}& \cfrac{\partial v}{\partial z} & \cfrac{\partial w}{\partial z}\\
  \end{bmatrix} 
  
-

.. math:: 
  \text{grad}\mathbf{v}=(\nabla \mathbf{v})^{\text T}=\begin{bmatrix}
  \cfrac{\partial u}{\partial x}& \cfrac{\partial v}{\partial x} & \cfrac{\partial w}{\partial x}\\
  \cfrac{\partial u}{\partial y}& \cfrac{\partial v}{\partial y} & \cfrac{\partial w}{\partial y}\\
  \cfrac{\partial u}{\partial z}& \cfrac{\partial v}{\partial z} & \cfrac{\partial w}{\partial z}\\
  \end{bmatrix}^{\text T} = \begin{bmatrix}
  \cfrac{\partial u}{\partial x}& \cfrac{\partial u}{\partial y} & \cfrac{\partial u}{\partial z}\\
  \cfrac{\partial v}{\partial x}& \cfrac{\partial v}{\partial y} & \cfrac{\partial v}{\partial z}\\
  \cfrac{\partial w}{\partial x}& \cfrac{\partial w}{\partial y} & \cfrac{\partial w}{\partial z}\\
  \end{bmatrix}
  
-

.. math:: 
  \nabla \mathbf{v}=\begin{bmatrix}
  \cfrac{\partial v_{1}}{\partial x_{1}}& \cfrac{\partial v_{2}}{\partial x_{1}} & \cfrac{\partial v_{3}}{\partial x_{1}}\\
  \cfrac{\partial v_{1}}{\partial x_{2}}& \cfrac{\partial v_{2}}{\partial x_{2}} & \cfrac{\partial v_{3}}{\partial x_{2}}\\
  \cfrac{\partial v_{1}}{\partial x_{3}}& \cfrac{\partial v_{2}}{\partial x_{3}} & \cfrac{\partial v_{3}}{\partial x_{3}}\\
  \end{bmatrix} \quad or\quad \left.(\nabla \mathbf{v})\right|_{i,j}=\cfrac{\partial v_{j}}{\partial x_{i}}  
  
Divergence of Tensor Functions
`````````````````````````````````````` 
the divergence of a second order tensor :math:`\mathbf{T}` is defined to be the vector

.. math::
  \begin{align}
  \text{div}\mathbf{T} & = \text{grad}\mathbf{T}:\mathbf{I}  = \cfrac{\partial \mathbf{T}}{\partial x_{i}}\mathbf{e}_{i}  = \cfrac{\partial ({T_{jk}\mathbf{e}_{j}\otimes \mathbf{e}_{k})}}{\partial x_{i}}\mathbf{e}_{i}  = \cfrac{\partial ({T_{jk})}}{\partial x_{i}}(\mathbf{e}_{j}\otimes \mathbf{e}_{k})\mathbf{e}_{i} \\
  &= \cfrac{\partial ({T_{jk})}}{\partial x_{i}}\mathbf{e}_{j}(\mathbf{e}_{k}\cdot\mathbf{e}_{i})\\
  &= \cfrac{\partial ({T_{ji})}}{\partial x_{i}}\mathbf{e}_{j}\\
  &= \cfrac{\partial ({T_{ij})}}{\partial x_{j}}\mathbf{e}_{i}\\
  \end{align}
  
One also has

.. math::
  \begin{align}
  \nabla \cdot \mathbf{T} & = \left ( \mathbf{e}_{i}\cfrac{\partial}{\partial x_{i}} \right )\cdot  \left (  {T_{jk}\mathbf{e}_{j}\otimes \mathbf{e}_{k}}\right ) \\
  &= \left ( \cfrac{\partial T_{jk}}{\partial x_{i}} \right )\mathbf{e}_{i}\cdot  \left (  {\mathbf{e}_{j}\otimes \mathbf{e}_{k}}\right ) \\
  &= \left ( \cfrac{\partial T_{jk}}{\partial x_{i}} \right )(\mathbf{e}_{i}\cdot \mathbf{e}_{j}) \mathbf{e}_{k} \\
  &= \left ( \cfrac{\partial T_{ik}}{\partial x_{i}} \right ) \mathbf{e}_{k} \\
  &= \cfrac{\partial\left ( { T_{ji}} \right )}{{\partial x_{j}}} \mathbf{e}_{i} \\
  \end{align}

so that

.. math::
  \text{div}\mathbf{T}  =\nabla \cdot \mathbf{T}^{\text{T}}
  
-
  
.. math::  
  \begin{align}
  \text{div}\mathbf{T} & =
  \cfrac{\partial ({T_{ij})}}{\partial x_{j}}\mathbf{e}_{i}=
  \left \{ \cfrac{\partial ({T_{1j})}}{\partial x_{j}}\right \}\mathbf{e}_{1}+
  \left \{ \cfrac{\partial ({T_{2j})}}{\partial x_{j}}\right \}\mathbf{e}_{2}+
  \left \{ \cfrac{\partial ({T_{3j})}}{\partial x_{j}}\right \}\mathbf{e}_{3} \\
  &=\left \{ \cfrac{\partial ({T_{11})}}{\partial x_{1}}
            +\cfrac{\partial ({T_{12})}}{\partial x_{2}}
            +\cfrac{\partial ({T_{13})}}{\partial x_{3}}\right \}\mathbf{e}_{1}\\
  &+\left \{ \cfrac{\partial ({T_{21})}}{\partial x_{1}}
            +\cfrac{\partial ({T_{22})}}{\partial x_{2}}
            +\cfrac{\partial ({T_{23})}}{\partial x_{3}}\right \}\mathbf{e}_{2}\\
  &+\left \{ \cfrac{\partial ({T_{31})}}{\partial x_{1}}
            +\cfrac{\partial ({T_{32})}}{\partial x_{2}}
            +\cfrac{\partial ({T_{33})}}{\partial x_{3}}\right \}\mathbf{e}_{3}\\
  \end{align}

-
  
.. math::  
  \text{div}\mathbf{T} = \begin{bmatrix}
   \cfrac{\partial {T_{11}}}{\partial x_{1}}
  +\cfrac{\partial {T_{12}}}{\partial x_{2}}
  +\cfrac{\partial {T_{13}}}{\partial x_{3}}\\
   \cfrac{\partial {T_{21}}}{\partial x_{1}}
  +\cfrac{\partial {T_{22}}}{\partial x_{2}}
  +\cfrac{\partial {T_{23}}}{\partial x_{3}}\\
   \cfrac{\partial {T_{31}}}{\partial x_{1}}
  +\cfrac{\partial {T_{32}}}{\partial x_{2}}
  +\cfrac{\partial {T_{33}}}{\partial x_{3}}\\
  \end{bmatrix}
  
-
  
.. math::  
  \begin{align}
  \nabla \cdot \mathbf{T} & =
  \cfrac{\partial\left ( { T_{ji}} \right )}{{\partial x_{j}}} \mathbf{e}_{i}=
  \left \{ \cfrac{\partial ({T_{j1})}}{\partial x_{j}}\right \}\mathbf{e}_{1}+
  \left \{ \cfrac{\partial ({T_{j2})}}{\partial x_{j}}\right \}\mathbf{e}_{2}+
  \left \{ \cfrac{\partial ({T_{j3})}}{\partial x_{j}}\right \}\mathbf{e}_{3} \\
  &=\left \{ \cfrac{\partial ({T_{11})}}{\partial x_{1}}
            +\cfrac{\partial ({T_{21})}}{\partial x_{2}}
            +\cfrac{\partial ({T_{31})}}{\partial x_{3}}\right \}\mathbf{e}_{1}\\
  &+\left \{ \cfrac{\partial ({T_{12})}}{\partial x_{1}}
            +\cfrac{\partial ({T_{22})}}{\partial x_{2}}
            +\cfrac{\partial ({T_{32})}}{\partial x_{3}}\right \}\mathbf{e}_{2}\\
  &+\left \{ \cfrac{\partial ({T_{13})}}{\partial x_{1}}
            +\cfrac{\partial ({T_{23})}}{\partial x_{2}}
            +\cfrac{\partial ({T_{33})}}{\partial x_{3}}\right \}\mathbf{e}_{3}\\
  \end{align}
  
-
  
.. math::    
  \nabla \cdot \mathbf{T} = \begin{bmatrix}
   \cfrac{\partial {T_{11}}}{\partial x_{1}}
  +\cfrac{\partial {T_{21}}}{\partial x_{2}}
  +\cfrac{\partial {T_{31}}}{\partial x_{3}}\\
   \cfrac{\partial {T_{12}}}{\partial x_{1}}
  +\cfrac{\partial {T_{22}}}{\partial x_{2}}
  +\cfrac{\partial {T_{32}}}{\partial x_{3}}\\
   \cfrac{\partial {T_{13}}}{\partial x_{1}}
  +\cfrac{\partial {T_{23}}}{\partial x_{2}}
  +\cfrac{\partial {T_{33}}}{\partial x_{3}}\\
  \end{bmatrix}
  
-
  
.. math::    
  \nabla \cdot \mathbf{T} = \begin{bmatrix}
   \cfrac{\partial {T_{xx}}}{\partial x}
  +\cfrac{\partial {T_{yx}}}{\partial y}
  +\cfrac{\partial {T_{zx}}}{\partial z}\\
   \cfrac{\partial {T_{xy}}}{\partial x}
  +\cfrac{\partial {T_{yy}}}{\partial y}
  +\cfrac{\partial {T_{zy}}}{\partial z}\\
   \cfrac{\partial {T_{xz}}}{\partial x}
  +\cfrac{\partial {T_{yz}}}{\partial y}
  +\cfrac{\partial {T_{zz}}}{\partial z}\\
  \end{bmatrix}    
  
-
  
.. math::    
  \text{div}\mathbf{T} = \nabla \cdot (\mathbf{T}^{\text T})=
  \begin{bmatrix}
   \cfrac{\partial {T_{11}}}{\partial x_{1}}
  +\cfrac{\partial {T_{12}}}{\partial x_{2}}
  +\cfrac{\partial {T_{13}}}{\partial x_{3}}\\
   \cfrac{\partial {T_{21}}}{\partial x_{1}}
  +\cfrac{\partial {T_{22}}}{\partial x_{2}}
  +\cfrac{\partial {T_{23}}}{\partial x_{3}}\\
   \cfrac{\partial {T_{31}}}{\partial x_{1}}
  +\cfrac{\partial {T_{32}}}{\partial x_{2}}
  +\cfrac{\partial {T_{33}}}{\partial x_{3}}\\
  \end{bmatrix}
  
 
-
  
.. math::   
  \text{div}\mathbf{T} = \nabla \cdot (\mathbf{T}^{\text T})=
  \begin{bmatrix}
   \cfrac{\partial {T_{xx}}}{\partial x}
  +\cfrac{\partial {T_{xy}}}{\partial y}
  +\cfrac{\partial {T_{xz}}}{\partial z}\\
   \cfrac{\partial {T_{yx}}}{\partial x}
  +\cfrac{\partial {T_{yy}}}{\partial y}
  +\cfrac{\partial {T_{yz}}}{\partial z}\\
   \cfrac{\partial {T_{zx}}}{\partial x}
  +\cfrac{\partial {T_{zy}}}{\partial y}
  +\cfrac{\partial {T_{zz}}}{\partial z}\\
  \end{bmatrix}  
  

 
The dot product of a tensor :math:`\boldsymbol{\tau}` by a vector :math:`\mathbf{v}` results in the following vector:  

.. math:: 
  [\boldsymbol\tau \cdot \mathbf{v}]& = \begin{pmatrix}
     \mathbf{ii}\tau_{xx}+\mathbf{ij}\tau_{xy}+\mathbf{ik}\tau_{xz}\\
   + \mathbf{ji}\tau_{yx}+\mathbf{jj}\tau_{yy}+\mathbf{jk}\tau_{yz}\\
   + \mathbf{ki}\tau_{zx}+\mathbf{kj}\tau_{zy}+\mathbf{kk}\tau_{zz}\\
  \end{pmatrix} 
  \cdot (u\mathbf{i}+v\mathbf{j}+w\mathbf{k})
  
which upon expanding becomes

.. math:: 
  \begin{align}
  [\boldsymbol\tau \cdot \mathbf{v}]& = \begin{pmatrix}
     \mathbf{ii}\tau_{xx}+\mathbf{ij}\tau_{xy}+\mathbf{ik}\tau_{xz}\\
   + \mathbf{ji}\tau_{yx}+\mathbf{jj}\tau_{yy}+\mathbf{jk}\tau_{yz}\\
   + \mathbf{ki}\tau_{zx}+\mathbf{kj}\tau_{zy}+\mathbf{kk}\tau_{zz}\\
  \end{pmatrix} 
  \cdot (u\mathbf{i}+v\mathbf{j}+w\mathbf{k})\\
  &=\begin{pmatrix}
     \mathbf{ii}\tau_{xx}+\mathbf{ij}\tau_{xy}+\mathbf{ik}\tau_{xz}\\
   + \mathbf{ji}\tau_{yx}+\mathbf{jj}\tau_{yy}+\mathbf{jk}\tau_{yz}\\
   + \mathbf{ki}\tau_{zx}+\mathbf{kj}\tau_{zy}+\mathbf{kk}\tau_{zz}\\
  \end{pmatrix} \cdot u\mathbf{i}\\
  &+\begin{pmatrix}
     \mathbf{ii}\tau_{xx}+\mathbf{ij}\tau_{xy}+\mathbf{ik}\tau_{xz}\\
   + \mathbf{ji}\tau_{yx}+\mathbf{jj}\tau_{yy}+\mathbf{jk}\tau_{yz}\\
   + \mathbf{ki}\tau_{zx}+\mathbf{kj}\tau_{zy}+\mathbf{kk}\tau_{zz}\\
  \end{pmatrix} \cdot v\mathbf{j}\\
  &+\begin{pmatrix}
     \mathbf{ii}\tau_{xx}+\mathbf{ij}\tau_{xy}+\mathbf{ik}\tau_{xz}\\
   + \mathbf{ji}\tau_{yx}+\mathbf{jj}\tau_{yy}+\mathbf{jk}\tau_{yz}\\
   + \mathbf{ki}\tau_{zx}+\mathbf{kj}\tau_{zy}+\mathbf{kk}\tau_{zz}\\
  \end{pmatrix} \cdot w\mathbf{k}\\
  \end{align}
  
- 
 
.. math:: 
  \begin{bmatrix}
  \mathbf{ii}\cdot \mathbf{i}=\mathbf{i}& \mathbf{ij}\cdot \mathbf{i}=\mathbf{0} & \mathbf{ik}\cdot \mathbf{i}=\mathbf{0}\\
  \mathbf{ji}\cdot \mathbf{i}=\mathbf{j}&  \mathbf{jj}\cdot \mathbf{i}=\mathbf{0}& \mathbf{jk}\cdot \mathbf{i}=\mathbf{0}\\
  \mathbf{ki}\cdot \mathbf{i}=\mathbf{k}&  \mathbf{kj}\cdot \mathbf{i}=\mathbf{0}&\mathbf{kk}\cdot \mathbf{i}=\mathbf{0}\\
  \mathbf{ii}\cdot \mathbf{j}=\mathbf{0}&\mathbf{ij}\cdot \mathbf{j}=\mathbf{i}&\mathbf{ik}\cdot \mathbf{j}=\mathbf{0}\\
  \mathbf{ji}\cdot \mathbf{j}=\mathbf{0}&\mathbf{jj}\cdot \mathbf{j}=\mathbf{j}&\mathbf{jk}\cdot \mathbf{j}=\mathbf{0}\\
  \mathbf{ki}\cdot \mathbf{j}=\mathbf{0}&\mathbf{kj}\cdot \mathbf{j}=\mathbf{k}&\mathbf{kk}\cdot \mathbf{j}=\mathbf{0}\\
  \mathbf{ii}\cdot \mathbf{k}=\mathbf{0}&\mathbf{ij}\cdot \mathbf{k}=\mathbf{0}&\mathbf{ik}\cdot \mathbf{k}=\mathbf{i}\\
  \mathbf{ji}\cdot \mathbf{k}=\mathbf{0}&\mathbf{jj}\cdot \mathbf{k}=\mathbf{0}&\mathbf{jk}\cdot \mathbf{k}=\mathbf{j}\\
  \mathbf{ki}\cdot \mathbf{k}=\mathbf{0}&\mathbf{kj}\cdot \mathbf{k}=\mathbf{0}&\mathbf{kk}\cdot \mathbf{k}=\mathbf{k}\\
  \end{bmatrix} 
  
- 
 
.. math:: 
  \begin{pmatrix}
     \mathbf{ii}\tau_{xx}+\mathbf{ij}\tau_{xy}+\mathbf{ik}\tau_{xz}\\
   + \mathbf{ji}\tau_{yx}+\mathbf{jj}\tau_{yy}+\mathbf{jk}\tau_{yz}\\
   + \mathbf{ki}\tau_{zx}+\mathbf{kj}\tau_{zy}+\mathbf{kk}\tau_{zz}\\
  \end{pmatrix} \cdot u\mathbf{i}=u(\tau_{xx}\mathbf{i}+\tau_{yx}\mathbf{j}+\tau_{zx}\mathbf{k}) 
  
- 
 
.. math::
  \begin{pmatrix}
   \mathbf{ii}\tau_{xx}+\mathbf{ij}\tau_{xy}+\mathbf{ik}\tau_{xz}\\
  + \mathbf{ji}\tau_{yx}+\mathbf{jj}\tau_{yy}+\mathbf{jk}\tau_{yz}\\
  + \mathbf{ki}\tau_{zx}+\mathbf{kj}\tau_{zy}+\mathbf{kk}\tau_{zz}\\
  \end{pmatrix} \cdot v\mathbf{j}=v(\tau_{xy}\mathbf{i}+\tau_{yy}\mathbf{j}+\tau_{zy}\mathbf{k})\\
  
- 
 
.. math::
  \begin{pmatrix}
   \mathbf{ii}\tau_{xx}+\mathbf{ij}\tau_{xy}+\mathbf{ik}\tau_{xz}\\
  + \mathbf{ji}\tau_{yx}+\mathbf{jj}\tau_{yy}+\mathbf{jk}\tau_{yz}\\
  + \mathbf{ki}\tau_{zx}+\mathbf{kj}\tau_{zy}+\mathbf{kk}\tau_{zz}\\
  \end{pmatrix} \cdot w\mathbf{k}=w(\tau_{xz}\mathbf{i}+\tau_{yz}\mathbf{j}+\tau_{zz}\mathbf{k})\\  
  
- 
 
.. math::
  \begin{align}
  [\boldsymbol\tau \cdot \mathbf{v}] & = u(\tau_{xx}\mathbf{i}+\tau_{yx}\mathbf{j}+\tau_{zx}\mathbf{k})\\
  &+v(\tau_{xy}\mathbf{i}+\tau_{yy}\mathbf{j}+\tau_{zy}\mathbf{k})\\
  &+w(\tau_{xz}\mathbf{i}+\tau_{yz}\mathbf{j}+\tau_{zz}\mathbf{k})\\
  &=(\tau_{xx}u+\tau_{xy}v+\tau_{xz}w)\mathbf{i}\\
  &+(\tau_{yx}u+\tau_{yy}v+\tau_{yz}w)\mathbf{j}\\
  &+(\tau_{zx}u+\tau_{zy}v+\tau_{zz}w)\mathbf{k}
  \end{align}
  
The above equation can be derived using matrix multiplication as 
 
.. math::
  [\boldsymbol{\tau}\cdot \mathbf{v}]=
  \begin{bmatrix}
    {\tau}_{xx}& {\tau}_{xy} & {\tau}_{xz}\\
    {\tau}_{yx}& {\tau}_{yy} & {\tau}_{yz}\\
    {\tau}_{zx}& {\tau}_{zy} & {\tau}_{zz}\\
  \end{bmatrix}
  \begin{bmatrix}
    u\\ v\\  w\\
  \end{bmatrix}
  =\begin{bmatrix}
    {\tau}_{xx}u+{\tau}_{xy}v+{\tau}_{xz}w\\
    {\tau}_{yx}u+{\tau}_{yy}v+{\tau}_{yz}w\\
    {\tau}_{zx}u+{\tau}_{zy}v+{\tau}_{zz}w\\
  \end{bmatrix}
  
In a similar way the divergence of a tensor :math:`\boldsymbol{\tau}` is found to be a vector given by  

.. math::
  \text{div }\boldsymbol{\tau} = \nabla \cdot (\boldsymbol{\tau}^{\text T})=
  \begin{bmatrix}
   \cfrac{\partial {{\tau}_{xx}}}{\partial x}
  +\cfrac{\partial {{\tau}_{xy}}}{\partial y}
  +\cfrac{\partial {{\tau}_{xz}}}{\partial z}\\
   \cfrac{\partial {{\tau}_{yx}}}{\partial x}
  +\cfrac{\partial {{\tau}_{yy}}}{\partial y}
  +\cfrac{\partial {{\tau}_{yz}}}{\partial z}\\
   \cfrac{\partial {{\tau}_{zx}}}{\partial x}
  +\cfrac{\partial {{\tau}_{zy}}}{\partial y}
  +\cfrac{\partial {{\tau}_{zz}}}{\partial z}\\
  \end{bmatrix}
  
-

.. math::
  [\nabla \cdot \boldsymbol{\tau}^{\text{T}}]=
  \begin{bmatrix}
    \cfrac{\partial }{\partial x}&\cfrac{\partial }{\partial y}  &\cfrac{\partial }{\partial z}
  \end{bmatrix}
  \begin{bmatrix}
    {\tau}_{xx}& {\tau}_{yx} & {\tau}_{zx}\\
    {\tau}_{xy}& {\tau}_{yy} & {\tau}_{zy}\\
    {\tau}_{xz}& {\tau}_{yz} & {\tau}_{zz}\\
  \end{bmatrix}
  =\begin{bmatrix}
    \cfrac{\partial{\tau}_{xx} }{\partial x}+
    \cfrac{\partial{\tau}_{xy} }{\partial y}+
    \cfrac{\partial{\tau}_{xz} }{\partial z}\\
    \cfrac{\partial{\tau}_{yx} }{\partial x}+
    \cfrac{\partial{\tau}_{yy} }{\partial y}+
    \cfrac{\partial{\tau}_{yz} }{\partial z}\\
    \cfrac{\partial{\tau}_{zx} }{\partial x}+
    \cfrac{\partial{\tau}_{zy} }{\partial y}+
    \cfrac{\partial{\tau}_{zz} }{\partial z}\\
  \end{bmatrix}  
  
-

.. math::
  \begin{align}
  [\nabla \cdot \boldsymbol{\tau}] & = (\cfrac{\partial{\tau}_{xx} }{\partial x}+
   \cfrac{\partial{\tau}_{yx} }{\partial y}+
   \cfrac{\partial{\tau}_{zx} }{\partial z})\mathbf{i}\\
  &+(\cfrac{\partial{\tau}_{xy} }{\partial x}+
   \cfrac{\partial{\tau}_{yy} }{\partial y}+
   \cfrac{\partial{\tau}_{zy} }{\partial z})\mathbf{j}\\
  &+(\cfrac{\partial{\tau}_{xz} }{\partial x}+
   \cfrac{\partial{\tau}_{yz} }{\partial y}+
   \cfrac{\partial{\tau}_{zz} }{\partial z})\mathbf{k}
  \end{align}

-

.. math::
  [\nabla \cdot \boldsymbol{\tau}]=
  \begin{bmatrix}
    \cfrac{\partial }{\partial x}&\cfrac{\partial }{\partial y}  &\cfrac{\partial }{\partial z}
  \end{bmatrix}
  \begin{bmatrix}
    {\tau}_{xx}& {\tau}_{xy} & {\tau}_{xz}\\
    {\tau}_{yx}& {\tau}_{yy} & {\tau}_{yz}\\
    {\tau}_{zx}& {\tau}_{zy} & {\tau}_{zz}\\
  \end{bmatrix}
  =\begin{bmatrix}
    \cfrac{\partial{\tau}_{xx} }{\partial x}+
    \cfrac{\partial{\tau}_{yx} }{\partial y}+
    \cfrac{\partial{\tau}_{zx} }{\partial z}\\
    \cfrac{\partial{\tau}_{xy} }{\partial x}+
    \cfrac{\partial{\tau}_{yy} }{\partial y}+
    \cfrac{\partial{\tau}_{zy} }{\partial z}\\
    \cfrac{\partial{\tau}_{xz} }{\partial x}+
    \cfrac{\partial{\tau}_{yz} }{\partial y}+
    \cfrac{\partial{\tau}_{zz} }{\partial z}\\
  \end{bmatrix}\begin{aligned}
  \end{aligned}
  
Double Contraction  
`````````````````````````````````
Double contraction, as the name implies, contracts the tensors twice as much a simple
contraction. Thus, where the sum of the orders of two tensors is reduced by two in the
simple contraction, the sum of the orders is reduced by four in double contraction. The
double contraction is denoted by a colon (:), e.g. :math:`\mathbf{A}:\mathbf{B}`.

.. math::
  \mathbf{A}:\mathbf{B}={A}_{ij}{B}_{ij}

Hadamard product 
`````````````````````````````````

#. `Hadamard product (matrices) <https://en.wikipedia.org/wiki/Hadamard_product_(matrices)>`_  


For two matrices A and B of the same dimension m × n, the Hadamard product 
is a matrix of the same dimension as the operands, with elements given by

.. math::
  (A\circ B)_{ij}=(A\odot B)_{ij}=(A)_{ij}(B)_{ij}
  
Hadamard Product or Element Wise Matrix Multiplication between any number of Matrices can be done if all the Matrices have Same Number of Rows and Same Number of Columns. The Resultant Matrix also has the Same Number of Rows and Columns as the input Matrices.  

.. math::
  A=\begin{bmatrix}
  a_{11}&a_{12}  & \cdots  &a_{1n} \\
  a_{21}&a_{22}  &  \cdots&a_{2n} \\
  \vdots & \vdots  & \ddots  & \vdots\\
  a_{m1}&a_{m2}  & \cdots &a_{mn}
  \end{bmatrix}
  
-
  
.. math::
  B=\begin{bmatrix}
  b_{11}&b_{12}  & \cdots  &b_{1n} \\
  b_{21}&b_{22}  &  \cdots &b_{2n} \\
  \vdots & \vdots  & \ddots  & \vdots\\
  b_{m1}&b_{m2}  & \cdots &b_{mn}
  \end{bmatrix}  
  
-
  
.. math::  
  \begin{align}
  A\circ B  = B\circ A & = \begin{bmatrix}
  a_{11}&a_{12}  & \cdots  &a_{1n} \\
  a_{21}&a_{22}  &  \cdots&a_{2n} \\
  \vdots & \vdots  & \ddots  & \vdots\\
  a_{m1}&a_{m2}  & \cdots &a_{mn}
  \end{bmatrix}\circ
 \begin{bmatrix}
  b_{11}&b_{12}  & \cdots  &b_{1n} \\
  b_{21}&b_{22}  &  \cdots&b_{2n} \\
  \vdots & \vdots  & \ddots  & \vdots\\
  b_{m1}&b_{m2}  & \cdots &b_{mn}
  \end{bmatrix} \\
  & = \begin{bmatrix}
  a_{11}\circ b_{11}&a_{12}\circ b_{12}  & \cdots  &a_{1n}\circ b_{1n} \\
  a_{21}\circ b_{21}&a_{22}\circ b_{22}  &  \cdots&a_{2n}\circ b_{2n} \\
  \vdots & \vdots  & \ddots  & \vdots\\
  a_{m1}\circ b_{m1}&a_{m2}\circ b_{m2}  & \cdots &a_{mn}\circ b_{mn}
  \end{bmatrix}
  \end{align}
  
Double-Dot Product of 2 Matrices
`````````````````````````````````
#. `Double-Dot Product of 2 Matrices <https://stemandmusic.in/maths/mvt-algebra/matrixDDP.php>`_

- Double-Dot Product between any 2 Matrices can be done if Both the Matrices have Same Number of Rows and Same Number of Columns. The Double-Dot Product of 2 Matrices is a Scalar Value.
- The Double-Dot Product of 2 Matrices is calculated by Calculating their Hadamard Product and Adding up all the Elements of the Resulting Matrix. 
- Given 2 :math:`m\times n` Matrices, Matrix :math:`A` having elements :math:`a_{ij}` Matrix :math:`B` having elements :math:`b_{ij}` as following

.. math::  
  A=\begin{bmatrix}
  a_{11}&a_{12}  & \cdots  &a_{1n} \\
  a_{21}&a_{22}  &  \cdots&a_{2n} \\
  \vdots & \vdots  & \ddots  & \vdots\\
  a_{m1}&a_{m2}  & \cdots &a_{mn}
  \end{bmatrix}
  \quad 
  B=\begin{bmatrix}
  b_{11}&b_{12}  & \cdots  &b_{1n} \\
  b_{21}&b_{22}  &  \cdots &b_{2n} \\
  \vdots & \vdots  & \ddots  & \vdots\\
  b_{m1}&b_{m2}  & \cdots &b_{mn}
  \end{bmatrix}  
  
The Double-Dot Product of Matrices :math:`A` and Matrices :math:`B` is calculated as

.. math::  
  \begin{align}
  \mathbf{A}: \mathbf{B}  = \mathbf{B}: \mathbf{A} =a_{ij}\circ b_{ij}& = a_{11} \circ b_{11}+a_{12} \circ b_{12}+\cdots+a_{1 n} \circ b_{1 n}\\
  &+a_{21} \circ b_{21}+a_{22} \circ b_{22}+\cdots+a_{2 n} \circ b_{2 n}\\
  &+\cdots\\
  &+a_{m 1} \circ b_{m 1}+a_{m 2} \circ b_{m 2}+\cdots+a_{m n} \circ b_{m n}
  \end{align}
  
-
  
.. math::  
  \begin{align}
  \mathbf{A}: \mathbf{B}  = \mathbf{B}: \mathbf{A} =a_{ij}\cdot b_{ij}& = a_{11}\cdot   b_{11}+a_{12} \cdot b_{12}+\cdots+a_{1 n} \cdot  b_{1 n}\\
  &+a_{21} \cdot  b_{21}+a_{22} \cdot  b_{22}+\cdots+a_{2 n} \cdot b_{2 n}\\
  &+\cdots\\
  &+a_{m 1} \cdot  b_{m 1}+a_{m 2} \cdot  b_{m 2}+\cdots+a_{m n} \cdot  b_{m n}
  \end{align} 

-
  
.. math::  
  \begin{align}
  \mathbf{A}: \mathbf{B}  = \mathbf{B}: \mathbf{A} =a_{ij} b_{ij}& = a_{11} b_{11}+a_{12}  b_{12}+\cdots+a_{1 n}  b_{1 n}\\
  &+a_{21} b_{21}+a_{22} b_{22}+\cdots+a_{2 n} \cdot b_{2 n}\\
  &+\cdots\\
  &+a_{m 1}  b_{m 1}+a_{m 2}   b_{m 2}+\cdots+a_{m n}  b_{m n}
  \end{align}
  
- Alternative Way to Calculate the Double-Dot Product of 2 Matrices is to find the Trace of following Inner/Dot Product of the 2 Matrices  

.. math:: 
  \begin{align}
  \mathbf{A}: \mathbf{B}  = \mathbf{B}: \mathbf{A}
   =Trace(\mathbf{A}\mathbf{B}^T)=Trace(\mathbf{B}^T\mathbf{A})
   =Trace(\mathbf{B}\mathbf{A}^T)=Trace(\mathbf{A}^T\mathbf{B})
  \end{align}
  
- 

.. math::  
  \begin{align}
  (\mathbf{a}\mathbf{b}):(\mathbf{c}\mathbf{d})=(\mathbf{a}\cdot \mathbf{c})(\mathbf{b}\cdot\mathbf{d})
  \end{align}  
  
Trace (linear algebra)
`````````````````````````````````
#. `Trace (linear algebra) <https://en.wikipedia.org/wiki/Trace_(linear_algebra)>`_

The trace of an :math:`n\times n` square matrix :math:`\mathbf{A}` is defined as

.. math::
  \mathrm{tr}(\mathbf{A})=a_{ii}=a_{11}+a_{22}+\cdots +a_{nn}

The double dot product of two tensors :math:`\boldsymbol{\tau}`  and :math:`\nabla \mathbf{v}` is a scalar computed as

The Divergence Theorem
`````````````````````````````````
Consider an arbitrary differentiable vector field :math:`\mathbf{v}(\mathbf{x},t)` defined in some finite region of
physical space. Let :math:`V` be a volume in this space with a closed surface :math:`S` bounding the
volume, and let the outward normal to this bounding surface be :math:`\mathbf{n}`. The **divergence
theorem of Gauss** states that (in symbolic and index notation)

.. math::
  \int\limits_{S} \mathbf{v}\cdot \mathbf{n}\text{d}S=\int\limits_{V} \text{div }\mathbf{v}\text{d}V
  
-

.. math::
  \int\limits_{S} v_{i}n_{i}\text{d}S = \int\limits_{V} \cfrac{\partial v_{i}}{\partial x_{i}}\text{d}V
  
and one has the following useful identities 

.. math::
  \begin{align}
  \int\limits_{S} \phi\mathbf{v}\cdot \mathbf{n}\text{d}S & = \int\limits_{V} \text{div}(\phi\mathbf{v})\text{d}V\\
  \int\limits_{S} \phi \mathbf{n}\text{d}S & = \int\limits_{V} \text{grad}(\phi)\text{d}V\\
  \int\limits_{S} \mathbf{n}\times \mathbf{v}\text{d}S & = \int\limits_{V} \text{curl}(\mathbf{v})\text{d}V\\
  \end{align}
  
The divergence theorem can be extended to the case of higher-order tensors.

Consider an arbitrary differentiable tensor field :math:`T_{ij\cdots k}(\mathbf{x},t)` defined in some finite region of
physical space. Let :math:`S` be a closed surface bounding a volume :math:`V` in this space, and let the
outward normal to :math:`S` be :math:`\mathbf{n}`. The divergence theorem of Gauss then states that  

.. math::
  \int\limits_{S} T_{i j \cdots k} n_{k} \text{d} S=\int\limits_{V} \frac{\partial T_{i j \cdots k}}{\partial x_{k}} \text{d} V
  
For a second order tensor,

.. math::
  \int\limits_{S} \mathbf{T}\cdot \mathbf{n} \text{d} S  = \int\limits_{V} \text{div } \mathbf{T} \text{d} V
  
-

.. math::
  \int\limits_{S} {T}_{ij}n_{ij} \text{d} S  = \int\limits_{V} \frac{\partial T_{ij}}{\partial x_{j}} \text{d} V
  
One then has the important identities

.. math::
  \begin{align}
  \int\limits_{S} (\phi\mathbf{T})\cdot \mathbf{n} \text{d} S  &= \int\limits_{V} \text{div }(\phi\mathbf{T}) \text{d} V\\
  \int\limits_{S} (\mathbf{u}\otimes \mathbf{n})\text{d} S  &= \int\limits_{V} \text{grad }\mathbf{u} \text{d} V\\
  \int\limits_{S} \mathbf{u}\cdot\mathbf{T}\cdot \mathbf{n} \text{d} S  &= \int\limits_{V} \text{div }(\mathbf{T}^{\text{T}}\cdot \mathbf{u}) \text{d} V\\
  \end{align}

Some formulation:

.. math::
  \text{div}(\mathbf{T}\cdot \mathbf{a})=\mathbf{a}\cdot\text{div}(\mathbf{T}^{\text{T}})+\text{tr}(\mathbf{T}\text{ grad }\mathbf{a})  
  
-
  
.. math::
  \mathbf{T}=\begin{bmatrix}
  T_{11}& T_{12} & T_{13}\\
  T_{21}& T_{22} & T_{23}\\
  T_{31}& T_{32} & T_{33}\\
  \end{bmatrix} \quad
  \mathbf{T}^{\text{T}}=\begin{bmatrix}
  T_{11}& T_{21} & T_{31}\\
  T_{12}& T_{22} & T_{32}\\
  T_{13}& T_{23} & T_{33}\\
  \end{bmatrix} \quad
  \mathbf{a}=\begin{bmatrix}
  a_{1}\\  a_{2}\\  a_{3}\\
  \end{bmatrix}
  
-
  
.. math::
  \text{grad } \mathbf{a}=\cfrac{\partial a_{i}}{\partial x_{j}}=\begin{bmatrix}
  \cfrac{\partial a_{1}}{\partial x_{1}}&
  \cfrac{\partial a_{1}}{\partial x_{2}}&
  \cfrac{\partial a_{1}}{\partial x_{3}}\\
  \cfrac{\partial a_{2}}{\partial x_{1}}&
  \cfrac{\partial a_{2}}{\partial x_{2}}&
  \cfrac{\partial a_{2}}{\partial x_{3}}\\
  \cfrac{\partial a_{3}}{\partial x_{1}}&
  \cfrac{\partial a_{3}}{\partial x_{2}}&
  \cfrac{\partial a_{3}}{\partial x_{3}}\\
  \end{bmatrix}  
  
-
  
.. math::
  \mathbf{T}\text{ grad } \mathbf{a}=T_{ik}\cfrac{\partial a_{k}}{\partial x_{j}}  
  
-
  
.. math::  
  \text{tr}(\mathbf{T}\text{ grad } \mathbf{a})=T_{ik}\cfrac{\partial a_{k}}{\partial x_{i}}=
  \begin{bmatrix}
  \quad(T_{11}\cfrac{\partial a_{1}}{\partial x_{1}}+
  T_{12}\cfrac{\partial a_{2}}{\partial x_{1}}+
  T_{13}\cfrac{\partial a_{3}}{\partial x_{1}})\\
  +(T_{21}\cfrac{\partial a_{1}}{\partial x_{2}}+
  T_{22}\cfrac{\partial a_{2}}{\partial x_{2}}+
  T_{23}\cfrac{\partial a_{3}}{\partial x_{2}})\\
  +(T_{31}\cfrac{\partial a_{1}}{\partial x_{3}}+
  T_{32}\cfrac{\partial a_{2}}{\partial x_{3}}+
  T_{33}\cfrac{\partial a_{3}}{\partial x_{3}})\\
  \end{bmatrix}

-
  
.. math::    
  \text{div}\mathbf{T}=\cfrac{\partial T_{ij}}{\partial x_{j}}=\begin{bmatrix}
  \cfrac{\partial T_{11}}{\partial x_{1}}
  +\cfrac{\partial T_{12}}{\partial x_{2}}
  +\cfrac{\partial T_{13}}{\partial x_{3}}\\
  \cfrac{\partial T_{21}}{\partial x_{1}}
  +\cfrac{\partial T_{22}}{\partial x_{2}}
  +\cfrac{\partial T_{23}}{\partial x_{3}}\\
  \cfrac{\partial T_{31}}{\partial x_{1}}
  +\cfrac{\partial T_{32}}{\partial x_{2}}
  +\cfrac{\partial T_{33}}{\partial x_{3}}\\
  \end{bmatrix}  
  
  
-
  
.. math::     
  \mathbf{T}^{\text{T}}=\begin{bmatrix}
  T_{11}& T_{21} & T_{31}\\
  T_{12}& T_{22} & T_{32}\\
  T_{13}& T_{23} & T_{33}\\
  \end{bmatrix} 
  
-
  
.. math::  
  \text{div}\mathbf{T}^{\text{T}}=\cfrac{\partial T_{ji}}{\partial x_{j}}=\begin{bmatrix}
   \cfrac{\partial T_{11}}{\partial x_{1}}
  +\cfrac{\partial T_{21}}{\partial x_{2}}
  +\cfrac{\partial T_{31}}{\partial x_{3}}\\
   \cfrac{\partial T_{12}}{\partial x_{1}}
  +\cfrac{\partial T_{22}}{\partial x_{2}}
  +\cfrac{\partial T_{32}}{\partial x_{3}}\\
   \cfrac{\partial T_{13}}{\partial x_{1}}
  +\cfrac{\partial T_{23}}{\partial x_{2}}
  +\cfrac{\partial T_{33}}{\partial x_{3}}\\
  \end{bmatrix}
  
-
  
.. math::   
  \mathbf{a}\cdot \text{div}\mathbf{T}^{\text{T}}=a_{i}\cfrac{\partial T_{ji}}{\partial x_{j}}=\begin{bmatrix}
  \quad a_{1}(\cfrac{\partial T_{11}}{\partial x_{1}}
  +\cfrac{\partial T_{21}}{\partial x_{2}}
  +\cfrac{\partial T_{31}}{\partial x_{3}})\\
  +a{2}(\cfrac{\partial T_{12}}{\partial x_{1}}
  +\cfrac{\partial T_{22}}{\partial x_{2}}
  +\cfrac{\partial T_{32}}{\partial x_{3}})\\
  +a{3}(\cfrac{\partial T_{13}}{\partial x_{1}}
  +\cfrac{\partial T_{23}}{\partial x_{2}}
  +\cfrac{\partial T_{33}}{\partial x_{3}})\\
  \end{bmatrix}
  
-
  
.. math:: 
  \begin{array}{c}
  \mathbf{T}\cdot\mathbf{a}=T_{ij}a_{j}\\
  \text{div }(\mathbf{T}\cdot\mathbf{a})=\cfrac{\partial (T_{ij}a_{j})}{\partial x_{i}}
  =T_{ij}\cfrac{\partial (a_{j})}{\partial x_{i}}+a_{j}\cfrac{\partial (T_{ij})}{\partial x_{i}}
  \end{array} 
  
-
  
.. math:: 
  \begin{align}
  \text{div }(\mathbf{T}\cdot\mathbf{a})
  & = \cfrac{\partial (T_{ij}a_{j})}{\partial x_{i}}\\
  & = T_{ij}\cfrac{\partial (a_{j})}{\partial x_{i}}+a_{j}\cfrac{\partial (T_{ij})}{\partial x_{i}}\\
  & = \text{tr}(\mathbf{T}\text{ grad } \mathbf{a})+\mathbf{a}\cdot \text{div}(\mathbf{T}^{\text{T}})
  \end{align}  