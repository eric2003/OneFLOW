Matrix
==============================================


Classical Adjoint (Adjugate) Matrix
----------------------------------------

Cofactor Formula
----------------------------------------
A Cofactor, in mathematics, is used to find the inverse of the matrix, adjoined. The Cofactor is the number you get when you remove the column and row of a designated element in a matrix, which is just a numerical grid in the form of rectangle or a square.
The cofactor is always preceded by a positive (+) or negative (-) sign. Let :math:`\mathbf{A}` be an :math:`n\times n` matrix and let
:math:`M_{ij}` be the :math:`(n-1)\times (n-1)` matrix obtained by deleting the :math:`i^{th}` row and :math:`j^{th}` column. Then,
:math:`\text{det}M_{ij}` is called the minor of :math:`a_{ij}`. The cofactor :math:`A_{ij}`of :math:`a_{ij}` is defined by:

.. math::
  A_{ij}=(-1)^{i+j}\text{det}M_{ij}
  
Minors and Cofactors  
-----------------------------
- `Minors and Cofactors <https://byjus.com/jee/minors-and-cofactors/>`_
- `Adjoint Of a Matrix <https://byjus.com/maths/adjoint-of-a-matrix/>`_


What Are Minors?
`````````````````````````
The minor of an element in a matrix is defined as the determinant obtained by deleting the row and column in which that element lies. For example, in the determinant

.. math::
  D=\begin{vmatrix}
  a_{11}& a_{12} & a_{13}\\
  a_{21}& a_{22} & a_{23}\\
  a_{31}& a_{32} & a_{33}\\
  \end{vmatrix}
  
minor of :math:`a_{12}` is denoted as :math:`M_{12}`. Here, 

.. math::
  M_{12}=\begin{vmatrix}
  a_{21}& a_{23}\\
  a_{31}& a_{33}\\
  \end{vmatrix}
  
What Are Cofactors?
`````````````````````````
Cofactor of an element aij is related to its minor as  

.. math::
  C_{ij}=(-1)^{i+j}M_{ij}
  
where :math:`i` denotes the :math:`i^{th}` row and :math:`j` denotes the :math:`i^{jth}` column to which the element :math:`a_{ij}` belongs.
Now, we define the value of the determinant of order three in terms of ‘Minor’ and ‘Cofactor’ as

.. math::
  D=a_{11}M_{11}-a_{12}M_{12}+a_{13}M_{13}
  
- 
 
.. math::
  D=a_{11}C_{11}+a_{12}C_{12}+a_{13}C_{13}  
  
The classical adjoint matrix should not be confused with the adjoint matrix. The adjoint is the conjugate transpose of a matrix while the classical adjoint is another name for the adjugate matrix or cofactor transpose of a matrix.

.. math::
  \mathbf{A}^{*}=\begin{bmatrix}
  A_{11}&A_{21}  &\cdots   & A_{n1}\\
  A_{12}&A_{22}  &\cdots   & A_{n2}\\
  \vdots& \vdots &  &\vdots \\
  A_{1n}&A_{2n}  &\cdots   & A_{nn}\\
  \end{bmatrix}  
  
Transpose
----------------------
`Transpose <https://en.wikipedia.org/wiki/Transpose/>`_

Formally, the :math:`i`-th row, :math:`j`-th column element of :math:`\mathbf{A}^{\text{T}}` is the :math:`j`-th row, :math:`i`-th column element of :math:`\mathbf{A}`:  

.. math::
  [\mathbf{A}^{\text{T}}]_{ij}=[\mathbf{A}]_{ji}
  
Properties
````````````````````

Let :math:`\mathbf{A}` and :math:`\mathbf{B}` be matrices and :math:`c` be a scalar.  

1. 

.. math::
  {\displaystyle \left(\mathbf {A} ^{\operatorname {T} }\right)^{\operatorname {T} }=\mathbf {A} .}
  
2. 

.. math::
  {\displaystyle \left(\mathbf {A} +\mathbf {B} \right)^{\operatorname {T} }=\mathbf {A} ^{\operatorname {T} }+\mathbf {B} ^{\operatorname {T} }.}
  
3. 

.. math::  
  {\displaystyle \left(\mathbf {AB} \right)^{\operatorname {T} }=\mathbf {B} ^{\operatorname {T} }\mathbf {A} ^{\operatorname {T} }.}
  
  
What are Eigenvalues?
---------------------------------------
The eigenvalue is explained to be a scalar associated with a linear set of equations which, when multiplied by a nonzero vector, equals to the vector obtained by transformation operating on the vector.

Let us consider :math:`k \times k` square matrix :math:`A` and :math:`\mathbf{v}` be a vector, then :math:`\lambda` is a scalar quantity represented in the following way:

.. math::
  A\mathbf{v} = \lambda\mathbf{v}

Here, :math:`\lambda` is considered to be the eigenvalue of matrix :math:`A`.

The above equation can also be written as:

.. math::
  (A – \lambda I) = 0

Where “:math:`I`” is the identity matrix of the same order as :math:`A`.

This equation can be represented in the determinant of matrix form.

.. math::
  |A – \lambda I| = 0
 
The above relation enables us to calculate eigenvalues :math:`\lambda` easily.

