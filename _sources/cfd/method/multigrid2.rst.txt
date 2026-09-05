Multigrid Poisson Solver
==================================

#. `Multigrid Solver for 1D Poisson Problem <https://people.math.sc.edu/Burkardt/c_src/multigrid_poisson_1d/multigrid_poisson_1d.html>`_
#. `Final Project: A Parallel Multigrid Poisson Solver <https://web.stanford.edu/class/cs315b/projects/multigrid_poisson/multigrid_poisson_slides.pdf>`_


Problem Description
--------------------------
We are interested in solving the Poisson equation, which can be expressed
in a domain :math:`\Omega\subset \mathbb{R}^{n}` with a boundary :math:`\partial\Omega` as

.. math::
  \left\{\begin{array}{l}
  -a\nabla^{2}\phi=f& \text{in }\Omega\\
  \phi=g& \text{on }\partial\Omega
  \end{array}\right.
  
where :math:`a` is a known constant, :math:`\nabla^{2}=\nabla\cdot\nabla(\cdot)` is the Laplacian operator, :math:`\phi=\phi(\vec{x})`
is our scalar variable that is a function of space, and :math:`f=f(\vec{x})` is a forcing function. The second line represents a Dirichlet
condition on the boundary :math:`\partial\Omega`, where :math:`g=g(\vec{x})` is a prescribed function along the boundary.

Discretization
--------------------------
Let :math:`\Omega` represent a 2D square domain given by :math:`[0, 1] \times [0, 1]`. We choose to discretize the above equation using a 2nd-order
finite difference approximation on a cartesian mesh composed of a number :math:`N` nodes in the x- and y-directions
with uniform spacing :math:`h`. The value of :math:`\psi` on the 2D cartesian mesh can then be approximated
for each node :math:`\{i, j\}` in the interior of the computational domain as

.. math::
  -\cfrac{a}{h^{2}} (\phi_{i,j-1}+\phi_{i-1,j}-4\phi_{i,j}+\phi_{i+1,j}+\phi_{i,j+1})=f_{i,j},\quad i,j=2,3,\cdots,N-1
  
where the subscripts :math:`i` and :math:`j` represent the indices of the current node in the computational domain for the
:math:`x`- and :math:`y` -directions, respectively. Assuming a Dirichlet boundary condition specifies the values for :math:`\phi` on the
boundary of the domain, we can consider the values for :math:`\phi_{1,j}`, :math:`\phi_{N,j}` , :math:`\phi_{i,1}`, and :math:`\phi_{i,N}`
known. The result is a system of :math:`(N − 2) \times (N − 2)` linear equations for the unknown
values of :math:`\phi_{i,j}` in the interior of the domain.

Forming the Linear System
--------------------------

.. figure:: ../images/multigrid6.png
   :width: 600
   :align: center
   
   2D uniform mesh featuring a discretization with a 5-point stencil.
   
We have freedom in deciding how to organize and ultimately solve the system of linear equations. Given a natural ordering of the unknowns and the 5-point stencil, we can express the full
linear system:  
 
.. math::
  A\Phi=\mathbf{f}
  
where the matrix :math:`A` is banded as a result of the structured mesh discretization and is given by  

.. figure:: ../images/multigrid7.png
   :width: 800
   :align: center
   
   Mesh
   
.. figure:: ../images/multigrid8.png
   :width: 800
   :align: center
   
   Mesh   
   
   
.. math::
  \cfrac{\phi_{i-1,j}-2\phi_{i,j}+\phi_{i+1,j}}{\Delta x^{2}}+\cfrac{\phi_{i,j-1}-2\phi_{i,j}+\phi_{i,j+1}}{\Delta y^{2}}=f_{i,j}

-
  
.. math::  
  \begin{array}{l}
  \Delta x=\cfrac{L_{x}}{M}=\cfrac{x_{\max}-x_{\min}}{M}=h_{x}\\
  \Delta y=\cfrac{L_{y}}{N}=\cfrac{y_{\max}-y_{\min}}{N}=h_{y}\\
  x_{i,j}=x_{\min}+i*h_{x},\quad i=0,1,\cdots,M\\
  y_{i,j}=y_{\min}+j*h_{y},\quad j=0,1,\cdots,N\\
  \end{array}  

Let 

.. math::  
  \begin{array}{l}
  L_{x}=L_{y}=1 \\
  M=N=4 \\
  \end{array}  
  
then  1d case

.. math::
  \cfrac{\phi_{i-1}-2\phi_{i}+\phi_{i+1}}{\Delta x^{2}}=f_{i}

-

.. math::
  \begin{array}{l}
  \cfrac{\phi_{0}-2\phi_{1}+\phi_{2}}{\Delta x^{2}}=f_{1}\\
  \cfrac{\phi_{1}-2\phi_{2}+\phi_{3}}{\Delta x^{2}}=f_{2}\\
  \cfrac{\phi_{2}-2\phi_{3}+\phi_{4}}{\Delta x^{2}}=f_{3}\\
  \end{array}  

-

.. math::
  \Delta x=\cfrac{1}{N}  
  
-

.. math::
  A=\begin{bmatrix}
  \cfrac{-2}{\Delta x^{2}}&\cfrac{1}{\Delta x^{2}}&0\\
  \cfrac{1}{\Delta x^{2}}&\cfrac{-2}{\Delta x^{2}}&\cfrac{1}{\Delta x^{2}}\\
  0&\cfrac{1}{\Delta x^{2}}&\cfrac{-2}{\Delta x^{2}}
  \end{bmatrix}=\cfrac{1}{\Delta x^{2}}\begin{bmatrix}
  {-2}&{1}&0\\
  {1}&{-2}&{1}\\
  0&{1}&{-2}
  \end{bmatrix}
  
-

.. math::
  \boldsymbol{\phi}=\begin{bmatrix}
  \phi_{1}\\
  \phi_{2}\\
  \phi_{3}\\
  \end{bmatrix}\quad 
  \mathbf{f}=\begin{bmatrix}
  f_{1}\\
  f_{2}\\
  f_{3}\\
  \end{bmatrix}  
  
2d case  

.. figure:: ../images/multigrid7.png
   :width: 800
   :align: center
   
   Mesh

.. math::
  \begin{array}{l}
  \cfrac{\phi_{i-1,j}-2\phi_{i,j}+\phi_{i+1,j}}{\Delta x^{2}}+\cfrac{\phi_{i,j-1}-2\phi_{i,j}+\phi_{i,j+1}}{\Delta y^{2}}=f_{i,j}\\
  \cfrac{\phi_{0,j}-2\phi_{1,j}+\phi_{2,j}}{\Delta x^{2}}+\cfrac{\phi_{1,j-1}-2\phi_{1,j}+\phi_{1,j+1}}{\Delta y^{2}}=f_{1,j}\\
  \cfrac{\phi_{1,j}-2\phi_{2,j}+\phi_{3,j}}{\Delta x^{2}}+\cfrac{\phi_{2,j-1}-2\phi_{2,j}+\phi_{2,j+1}}{\Delta y^{2}}=f_{2,j}\\
  \cfrac{\phi_{2,j}-2\phi_{3,j}+\phi_{4,j}}{\Delta x^{2}}+\cfrac{\phi_{3,j-1}-2\phi_{3,j}+\phi_{3,j+1}}{\Delta y^{2}}=f_{3,j}\\
  \end{array}
  
:math:`j=1`:

.. math::
  \begin{array}{l}
  \cfrac{\phi_{0,1}-2\phi_{1,1}+\phi_{2,1}}{\Delta x^{2}}+\cfrac{\phi_{1,0}-2\phi_{1,1}+\phi_{1,2}}{\Delta y^{2}}=f_{1,1}\\
  \cfrac{\phi_{1,1}-2\phi_{2,1}+\phi_{3,1}}{\Delta x^{2}}+\cfrac{\phi_{2,0}-2\phi_{2,1}+\phi_{2,2}}{\Delta y^{2}}=f_{2,1}\\
  \cfrac{\phi_{2,1}-2\phi_{3,1}+\phi_{4,1}}{\Delta x^{2}}+\cfrac{\phi_{3,0}-2\phi_{3,1}+\phi_{3,2}}{\Delta y^{2}}=f_{3,1}\\
  \end{array}
  
:math:`j=2`:  

.. math::
  \begin{array}{l}
  \cfrac{\phi_{0,2}-2\phi_{1,2}+\phi_{2,2}}{\Delta x^{2}}+\cfrac{\phi_{1,1}-2\phi_{1,2}+\phi_{1,3}}{\Delta y^{2}}=f_{1,2}\\
  \cfrac{\phi_{1,2}-2\phi_{2,2}+\phi_{3,2}}{\Delta x^{2}}+\cfrac{\phi_{2,1}-2\phi_{2,2}+\phi_{2,3}}{\Delta y^{2}}=f_{2,2}\\
  \cfrac{\phi_{2,2}-2\phi_{3,2}+\phi_{4,2}}{\Delta x^{2}}+\cfrac{\phi_{3,1}-2\phi_{3,2}+\phi_{3,3}}{\Delta y^{2}}=f_{3,2}\\
  \end{array}
  
:math:`j=3`: 

.. math::
  \begin{array}{l}
  \cfrac{\phi_{0,3}-2\phi_{1,3}+\phi_{2,3}}{\Delta x^{2}}+\cfrac{\phi_{1,2}-2\phi_{1,3}+\phi_{1,4}}{\Delta y^{2}}=f_{1,3}\\
  \cfrac{\phi_{1,3}-2\phi_{2,3}+\phi_{3,3}}{\Delta x^{2}}+\cfrac{\phi_{2,2}-2\phi_{2,3}+\phi_{2,4}}{\Delta y^{2}}=f_{2,3}\\
  \cfrac{\phi_{2,3}-2\phi_{3,3}+\phi_{4,3}}{\Delta x^{2}}+\cfrac{\phi_{3,2}-2\phi_{3,3}+\phi_{3,4}}{\Delta y^{2}}=f_{3,3}\\
  \end{array}  
  
-
  
.. math::
  \Phi=\begin{bmatrix}
  {\phi}_{1,1}\\{\phi}_{2,1}\\{\phi}_{3,1}\\
  {\phi}_{1,2}\\{\phi}_{2,2}\\{\phi}_{3,2}\\
  {\phi}_{1,3}\\{\phi}_{2,3}\\{\phi}_{3,3}\\
  \end{bmatrix}=
  \begin{bmatrix}
   \hat{\phi}_{1}\\\hat{\phi}_{2} \\\hat{\phi}_{3}\\
   \hat{\phi}_{4}\\\hat{\phi}_{5} \\\hat{\phi}_{6}\\
   \hat{\phi}_{7}\\\hat{\phi}_{8} \\\hat{\phi}_{9}\\
  \end{bmatrix}
  \quad \mathbf{f}=\begin{bmatrix}
  {f}_{1,1}\\{f}_{2,1}\\{f}_{3,1}\\
  {f}_{1,2}\\{f}_{2,2}\\{f}_{3,2}\\
  {f}_{1,3}\\{f}_{2,3}\\{f}_{3,3}\\
  \end{bmatrix}=
  \begin{bmatrix}
   \hat{f}_{1}\\\hat{f}_{2} \\\hat{f}_{3}\\
   \hat{f}_{4}\\\hat{f}_{5} \\\hat{f}_{6}\\
   \hat{f}_{7}\\\hat{f}_{8} \\\hat{f}_{9}\\
  \end{bmatrix}
  \quad A\Phi=\mathbf{f}
  
-
  
.. math::
  A=\begin{bmatrix}
  a_{11}&a_{12}&a_{13}&a_{14}&a_{15}&a_{16}&a_{17}&a_{18}&a_{19} \\
  a_{21}&a_{22}&a_{23}&a_{24}&a_{25}&a_{26}&a_{27}&a_{28}&a_{29} \\
  a_{31}&a_{32}&a_{33}&a_{34}&a_{35}&a_{36}&a_{37}&a_{38}&a_{39} \\
  a_{41}&a_{42}&a_{43}&a_{44}&a_{45}&a_{46}&a_{47}&a_{48}&a_{49} \\
  a_{51}&a_{52}&a_{53}&a_{54}&a_{55}&a_{56}&a_{57}&a_{58}&a_{59} \\
  a_{61}&a_{62}&a_{63}&a_{64}&a_{65}&a_{66}&a_{67}&a_{68}&a_{69} \\
  a_{71}&a_{72}&a_{73}&a_{74}&a_{75}&a_{76}&a_{77}&a_{78}&a_{79} \\
  a_{81}&a_{82}&a_{83}&a_{84}&a_{85}&a_{86}&a_{87}&a_{88}&a_{89} \\
  a_{91}&a_{92}&a_{93}&a_{94}&a_{95}&a_{96}&a_{97}&a_{98}&a_{99} \\
  \end{bmatrix}  
  
or
  
.. math::
  A=\begin{bmatrix}
  a_{1,(1,1)}&a_{1,(2,1)}&a_{1,(3,1)}&a_{1,(1,2)}&a_{1,(2,2)}&a_{1,(3,2)}&a_{1,(1,3)}&a_{1,(2,3)}&a_{1,(3,3)} \\
  a_{2,(1,1)}&a_{2,(2,1)}&a_{2,(3,1)}&a_{2,(1,2)}&a_{2,(2,2)}&a_{2,(3,2)}&a_{2,(1,3)}&a_{2,(2,3)}&a_{2,(3,3)} \\
  a_{3,(1,1)}&a_{3,(2,1)}&a_{3,(3,1)}&a_{3,(1,2)}&a_{3,(2,2)}&a_{3,(3,2)}&a_{3,(1,3)}&a_{3,(2,3)}&a_{3,(3,3)} \\
  a_{4,(1,1)}&a_{4,(2,1)}&a_{4,(3,1)}&a_{4,(1,2)}&a_{4,(2,2)}&a_{4,(3,2)}&a_{4,(1,3)}&a_{4,(2,3)}&a_{4,(3,3)} \\
  a_{5,(1,1)}&a_{5,(2,1)}&a_{5,(3,1)}&a_{5,(1,2)}&a_{5,(2,2)}&a_{5,(3,2)}&a_{5,(1,3)}&a_{5,(2,3)}&a_{5,(3,3)} \\
  a_{6,(1,1)}&a_{6,(2,1)}&a_{6,(3,1)}&a_{6,(1,2)}&a_{6,(2,2)}&a_{6,(3,2)}&a_{6,(1,3)}&a_{6,(2,3)}&a_{6,(3,3)} \\
  a_{7,(1,1)}&a_{7,(2,1)}&a_{7,(3,1)}&a_{7,(1,2)}&a_{7,(2,2)}&a_{7,(3,2)}&a_{7,(1,3)}&a_{7,(2,3)}&a_{7,(3,3)} \\
  a_{8,(1,1)}&a_{8,(2,1)}&a_{8,(3,1)}&a_{8,(1,2)}&a_{8,(2,2)}&a_{8,(3,2)}&a_{8,(1,3)}&a_{8,(2,3)}&a_{8,(3,3)} \\
  a_{9,(1,1)}&a_{9,(2,1)}&a_{9,(3,1)}&a_{9,(1,2)}&a_{9,(2,2)}&a_{9,(3,2)}&a_{9,(1,3)}&a_{9,(2,3)}&a_{9,(3,3)} \\
  \end{bmatrix}  
  
-
  
.. math::
  \begin{array}{l}
  \cfrac{\phi_{i-1,j}-2\phi_{i,j}+\phi_{i+1,j}}{\Delta x^{2}}+\cfrac{\phi_{i,j-1}-2\phi_{i,j}+\phi_{i,j+1}}{\Delta y^{2}}=f_{i,j}\\
  \cfrac{1}{\Delta y^{2}}\phi_{i,j-1}+\cfrac{1}{\Delta x^{2}} \phi_{i-1,j}-\bigg(\cfrac{2}{\Delta x^{2}}+\cfrac{2}{\Delta y^{2}} \bigg)\phi_{i,j}
  +\cfrac{1}{\Delta x^{2}} \phi_{i+1,j}+\cfrac{1}{\Delta y^{2}} \phi_{i,j+1}=f_{i,j}\\
  a\phi_{i,j-1}+b \phi_{i-1,j}+c\phi_{i,j}
  +d \phi_{i+1,j}+e\phi_{i,j+1}=f_{i,j}\\
  a=\cfrac{1}{\Delta y^{2}},b=\cfrac{1}{\Delta x^{2}},c=-\bigg(\cfrac{2}{\Delta x^{2}}+\cfrac{2}{\Delta y^{2}} \bigg),  d=\cfrac{1}{\Delta x^{2}},e=\cfrac{1}{\Delta y^{2}}\\
  \end{array}
  
-
  
.. math::
  \begin{array}{l}
  a\phi_{i,j-1}+b \phi_{i-1,j}+c\phi_{i,j}+d \phi_{i+1,j}+e\phi_{i,j+1}=f_{i,j}\\
  a\phi_{1,0}+b \phi_{0,1}+c\phi_{1,1}+d \phi_{2,1}+e\phi_{1,2}=f_{1,1}\\
  a\phi_{2,0}+b \phi_{1,1}+c\phi_{2,1}+d \phi_{3,1}+e\phi_{2,2}=f_{2,1}\\
  a\phi_{3,0}+b \phi_{2,1}+c\phi_{3,1}+d \phi_{4,1}+e\phi_{3,2}=f_{3,1}\\
  a\phi_{1,1}+b \phi_{0,2}+c\phi_{1,2}+d \phi_{2,2}+e\phi_{1,3}=f_{1,2}\\
  a\phi_{2,1}+b \phi_{1,2}+c\phi_{2,2}+d \phi_{3,2}+e\phi_{2,3}=f_{2,2}\\
  a\phi_{3,1}+b \phi_{2,2}+c\phi_{3,2}+d \phi_{4,2}+e\phi_{3,3}=f_{3,2}\\
  a\phi_{1,2}+b \phi_{0,3}+c\phi_{1,3}+d \phi_{2,3}+e\phi_{1,4}=f_{1,3}\\
  a\phi_{2,2}+b \phi_{1,3}+c\phi_{2,3}+d \phi_{3,3}+e\phi_{2,4}=f_{2,3}\\
  a\phi_{3,2}+b \phi_{2,3}+c\phi_{3,3}+d \phi_{4,3}+e\phi_{3,4}=f_{3,3}\\
  \end{array}  
  
-
  
.. math::
  A=\begin{bmatrix}
  c&d&0&e&0&0&0&0&0 \\
  b&c&d&0&e&0&0&0&0 \\
  0&b&c&0&0&e&0&0&0 \\
  a&0&0&c&d&0&e&0&0 \\
  0&a&0&b&c&d&0&e&0 \\
  0&0&a&0&b&c&0&0&e \\
  0&0&0&a&0&0&c&d&0 \\
  0&0&0&0&a&0&b&c&d \\
  0&0&0&0&0&a&0&b&c \\
  \end{bmatrix}    
  
Classical Iterative Methods
--------------------------------
At this point, a number of methods can be employed for solving the linear system in Eqn. 3, but we will
focus on the classical iterative methods in this project. The Gauss-Seidel will serve as an example in the
discussion below. Starting with a linear system

Geometric Multigrid
For typical iterative numerical solution methods, high-frequency (local) errors in the solution are well-damped, while lower frequency (global) errors are poorly damped. Therefore, the low-frequency errors are
difficult to eliminate, which leads to slower solver convergence, especially on fine meshes. The key idea behind
multigrid is that effective rates of convergence at all scales can be maintained in a solver by leveraging a
sequence of grids at various resolutions. With geometric multigrid, multiple levels of physical grids with
varying resolution are used to provide better approximations of the solution with each step of an iterative
solution method (i.e., a multigrid cycle).

To illustrate the basic components of linear multigrid for elliptic problems, define the error in the
solution to be the difference between the solution :math:`\Phi` and the approximation to the solution :math:`\tilde{\Phi}`, or

.. math::
  \mathbf{e}={\Phi} -\tilde{\Phi} 

where :math:`e` is the error vector (one value per node in the computational mesh). We can also define a residual
vector :math:`r`, which is a measure of how well the discretized governing equations are being satisfied by our
numerical solution procedure, as

.. math::
  \mathbf{r}=\mathbf{f}-A\tilde{\Phi} 
  
-
  
.. math::
  \begin{array}{l}
  \mathbf{e}={\Phi} -\tilde{\Phi} \\
  A(\tilde{\Phi}+\mathbf{e})=\mathbf{f}\\
  A\mathbf{e}=\mathbf{r}\\
  \end{array}  
  
which relates the error in the solution to the residual. Eqn. :math:`A\mathbf{e}=\mathbf{r}` allows us to compute a measure of the error
on coarser mesh levels after transferring the values of the residual from the fine mesh level onto the coarse
level (restriction). After calculating e on a coarse level, we can form a correction to the solution on the fine
mesh as  

.. math::
  {\Phi}=\tilde{\Phi}+\mathbf{e}

upon transferring the error up to a fine mesh from the coarse mesh level below (prolongation). Furthermore,
we can apply these ideas recursively over an entire set of grids of various resolutions to complete a full
multigrid cycle, such as the V-cycle detailed in Alg. 1.

During a multigrid V-cycle, the solution is first approximated using several smoothing iterations with a
method like Gauss-Seidel on the finest mesh (pre-smoothing), and then the residual is transferred to the first
coarse level, where additional smoothing iterations occur. This restriction followed by smoothing continues
recursively until the coarsest mesh level is reached (the downstroke of the cycle). After performing some
smoothing iterations on the coarsest level, a correction for the solution values is transferred to the finer mesh
level above. This upward stroke of the cycle with prolongation and smoothing continues until a correction
is finally applied to the solution on the finest mesh. Typically, several final smoothing iterations (post-smoothing) are performed on the finest mesh before moving on to the next multigrid cycle. The downstroke
and upstroke of the cycle form a V-shape when viewed graphically, as in Fig. 3. Other cycles are possible,
and W-cycles are common, for instance.

Algorithm 1 Multigrid V-Cycle
-----------------------------------
1. :math:`\text{procedure MULTIGRID_CYCLE}(\Phi,\mathbf{f},l), \quad \text{ l is the current mesh level}`
2. :math:`\quad\quad\Phi^{l}\leftarrow \text{GAUSS_SEIDEL}(\Phi^{l},\mathbf{f}^{l}) \quad \text{ Pre-smoothing of the solution }`
3. :math:`\quad\quad\text{if l< n_levels then}`
4. :math:`\quad\quad\quad\Phi^{l+1}\leftarrow 0`
5. :math:`\quad\quad\quad\mathbf{f}^{l+1}\leftarrow \text{RESTRICT}(\Phi^{l},\mathbf{f}^{l})`
6. :math:`\quad\quad\quad\text{MULTIGRID_CYCLE}(\Phi,\mathbf{f},l+1)\quad\quad\quad \text{Recursive call}`
7. :math:`\quad\quad\quad\Phi^{l}\leftarrow\text{PROLONGATE}(\Phi^{l+1},\mathbf{f}^{l+1})`
8. :math:`\quad\quad\Phi^{l}\leftarrow \text{GAUSS_SEIDEL}(\Phi^{l},\mathbf{f}^{l}) \quad \text{ Post-smoothing of the solution }`


Restriction Operator
-----------------------------------
For our discretized Poisson problem, we can express the value of the residual :math:`\mathbf{r}` from Eqn. :math:`\mathbf{r}=\mathbf{f}-A\tilde{\Phi}` at each node as

.. math::
  \begin{array}{l}
  -\cfrac{a}{h^{2}} (\phi_{i,j-1}+\phi_{i-1,j}-4\phi_{i,j}+\phi_{i+1,j}+\phi_{i,j+1})=f_{i,j}\\
  -(\phi_{i,j-1}+\phi_{i-1,j}-4\phi_{i,j}+\phi_{i+1,j}+\phi_{i,j+1})=\cfrac{h^{2}}{a}f_{i,j}=\tilde{f}_{i,j}\\
  A\Phi=\mathbf{f}\\
  \mathbf{r}=\mathbf{f}-A\Phi\\
  r_{i,j}=\cfrac{h^{2}}{a}f_{i,j}+(\phi_{i,j-1}+\phi_{i-1,j}-4\phi_{i,j}+\phi_{i+1,j}+\phi_{i,j+1})
  \end{array}
  
After computing :math:`\mathbf{r}`, we restrict the values down to the next coarse level to form the right-hand side of
Eqn. :math:`A\mathbf{e}=\mathbf{r}`. For a weighted restriction, we will include information from all of the fine nodes that surround a
particular coarse mesh node. To accomplish this, we will set the residual at a coarse mesh node to be the
sum of a contribution from the coincident node (1/4 of the value), the nodes that are part of the stencil in
the north, south, east, and west directions (1/8 of the value), and the diagonal neighbors, i.e., north-east,
north-west, south-east, and south-west (1/16 of the value).  

Prolongation Operator
-----------------------------
A prolongation operation is one that transfers the correction from a coarse mesh to a fine mesh. Similar to
the weighted method for restriction, we will perform a weighted prolongation by setting the correction at a
fine mesh node to be the value of the correction at a coincident coarse node, if applicable, or as the sum of
a contribution from the coarse nodes that are nearest in the north, south, east, and west directions (1/2 of
the value) and the nearest diagonal neighbors, i.e., north-east, north-west, south-east, and south-west (1/4
of the value).

Restriction Operator Continue 
------------------------------------
The second class of intergrid transfer operations involves moving vectors from
a fine grid to a coarse grid. They are generally called restriction operators and are
denoted by :math:`I_{h}^{2h}`. The most obvious restriction operator is injection. It is defined by
:math:`I_{h}^{2h}\mathbf{v}^{h}=\mathbf{v}^{2h}`, where

.. math::
  v_{j}^{2h}=v_{2j}^{h}
  
.. figure:: ../images/multigrid9.png
   :width: 800
   :align: center
   
   1d case
  
.. math::
  \begin{array}{l}
  v_{0}^{2h}=v_{0}^{h}\\
  v_{1}^{2h}=v_{2}^{h}\\
  v_{2}^{2h}=v_{4}^{h}\\
  \end{array}  
  
In other words, with injection, the coarse-grid vector simply takes its value directly
from the corresponding fine-grid point.
An alternate restriction operator, called full weighting, is defined by  :math:`I_{h}^{2h}\mathbf{v}^{h}=\mathbf{v}^{2h}`, where

.. math::
  v_{j}^{2h}=\cfrac{1}{4}(v_{2j-1}^{h}+2v_{2j}^{h}+v_{2j+1}^{h}),\quad 1\le j \le \cfrac{n}{2}-1.
  
-

.. math::
  \begin{array}{l}
  v_{j}^{2h}=\cfrac{1}{4}(v_{2j-1}^{h}+2v_{2j}^{h}+v_{2j+1}^{h}),\quad 1\le j \le \cfrac{n}{2}-1.\\
  v_{0}^{2h}=\cfrac{1}{4}(v_{-1}^{h}+2v_{0}^{h}+v_{1}^{h})\\
  v_{1}^{2h}=\cfrac{1}{4}(v_{1}^{h}+2v_{2}^{h}+v_{3}^{h})\\
  v_{2}^{2h}=\cfrac{1}{4}(v_{3}^{h}+2v_{4}^{h}+v_{5}^{h})\\
  \end{array}  
  
The fact that the interpolation operator and the full weighting operator are transposes of each other up to a constant is called a variational property and will soon
be of importance.

For the sake of completeness, we give the full weighting operator in two dimensions. It is just an averaging of the fine-grid nearest neighbors. Letting :math:`I_{h}^{2h}\mathbf{v}^{h}=\mathbf{v}^{2h}`,
we have that  

.. math::
  \begin{array}{l}
  v_{i,j}^{2h}=\cfrac{1}{16}[v_{2i-1,2j-1}^{h}+v_{2i-1,2j+1}^{h}+v_{2i+1,2j-1}^{h}+v_{2i+1,2j+1}^{h}\\
  +2(v_{2i,2j-1}^{h}+v_{2i,2j+1}^{h}+v_{2i+1,2j}^{h}+v_{2i-1,2j}^{h})\\
  +4v_{2i,2j}^{h}],\quad 1\le i,j \le \cfrac{n}{2}-1.\\
  \end{array}
  
-

.. math::
  \begin{bmatrix}
  \cfrac{1}{16}v_{2i-1,2j-1}^{h}& \cfrac{1}{8}v_{2i,2j-1}^{h} &\cfrac{1}{16}v_{2i+1,2j-1}^{h} \\
  \cfrac{1}{8}v_{2i-1,2j}^{h} & \cfrac{1}{4}v_{2i,2j}^{h}& \cfrac{1}{8}v_{2i+1,2j}^{h}\\
  \cfrac{1}{16}v_{2i-1,2j+1}^{h}& \cfrac{1}{8}v_{2i,2j+1}^{h} &\cfrac{1}{16}v_{2i+1,2j+1}^{h}
  \end{bmatrix}  
  
-

.. math::
  \begin{bmatrix}
  {2i-1,2j+1}& {2i,2j+1} &{2i+1,2j+1} \\
  {2i-1,2j} & {2i,2j}& {2i+1,2j}\\
  {2i-1,2j-1}& {2i,2j-1} &{2i+1,2j-1} \\
  \end{bmatrix}
  
index

.. math::
  \begin{bmatrix}
  {2i-1,2j-1}& {2i,2j-1} &{2i+1,2j-1} \\
  {2i-1,2j} & {2i,2j}& {2i+1,2j}\\
  {2i-1,2j+1}& {2i,2j+1} &{2i+1,2j+1}
  \end{bmatrix}  
  
.. figure:: ../images/multigrid10.png
   :width: 800
   :align: center
   
   2d case  
  
(1,1)

.. math::
  \begin{bmatrix}
  {1,3}& {2,3} &{3,3}\\
  {1,2} & {2,2}& {3,2}\\
  {1,1}& {2,1} &{3,1} \\  
  \end{bmatrix}
  
left bottom corner point (0,0)

.. math::
  \begin{bmatrix}
  {-1,1}& {0,1} &{1,1}\\  
  {-1,0} & {0,0}& {1,0}\\
  {-1,-1}& {0,-1} &{1,-1} \\  
  \end{bmatrix}  
  
bottom boundary (1,0)

.. math::
  \begin{bmatrix}
  {1,1}& {2,1} &{3,1}\\  
  {1,0} & {2,0}& {3,0}\\
  {1,-1}& {2,-1} &{3,-1} \\
  \end{bmatrix}  
  
right bottom corner point (2,0)

.. math::
  \begin{bmatrix}
  {3,1}& {4,1} &{5,1}\\  
  {3,0}& {4,0}& {5,0}\\
  {3,-1}& {4,-1} &{5,-1} \\
  \end{bmatrix}  
  
left top corner point (0,2)

.. math::
  \begin{bmatrix}
  {-1,5}& {0,5} &{1,5}\\  
  {-1,4} & {0,4}& {1,4}\\
  {-1,3}& {0,3} &{1,3} \\
  \end{bmatrix}  
  
  
right top corner point (2,2)  

.. math::
  \begin{bmatrix}
  {3,5}& {4,5}&{5,5}\\
  {3,4}& {4,4}&{5,4}\\
  {3,3}& {4,3}&{5,3}\\
  \end{bmatrix}
  
top boundary (1,2)  

.. math::
  \begin{bmatrix}
  {1,5}& {2,5} &{3,5} \\
  {1,4} & {2,4}& {3,4}\\
  {1,3}& {2,3} &{3,3} \\
  \end{bmatrix}
  
left boundary (0,1)  

.. math::
  \begin{bmatrix}
  {-1,3}& {0,3} &{1,3} \\
  {-1,2} & {0,2}& {1,2}\\
  {-1,1}& {0,1} &{1,1} \\
  \end{bmatrix}
  
  
right boundary (2,1)  

.. math::
  \begin{bmatrix}
  {3,3}& {4,3} &{5,3} \\
  {3,2} & {4,2}& {5,2}\\
  {3,1}& {4,1} &{5,1} \\
  \end{bmatrix}
  
Interpolation or Prolongation
--------------------------------------
In our discussion of intergrid transfers, we consider only the case in which the
coarse grid has twice the grid spacing of the next finest grid. This is a nearly
universal practice, because there is usually no advantage in using grid spacings
with ratios other than 2. Think for a moment about the step in the correction
scheme that requires transferring the error approximation :math:`\mathbf{e}^{2h}` from the coarse grid
:math:`\Omega^{2h}` to the fine grid :math:`\Omega^{h}`. This is a common procedure in numerical analysis and is
generally called interpolation or prolongation. Many interpolation methods could
be used. Fortunately, for most multigrid purposes, the simplest of these is quite
effective. For this reason, we consider only linear interpolation.

The linear interpolation operator will be denoted :math:`I_{2h}^{h}`. It takes coarse-grid vectors and produces fine-grid vectors according to the rule :math:`I_{2h}^{h}\mathbf{v}^{2h}=\mathbf{v}^{h}`, where

.. math::
  \begin{array}{l}
  v_{2j}^{h}=v_{j}^{2h}\\
  v_{2j+1}^{h}=\cfrac{1}{2}( v_{j}^{2h}+v_{j+1}^{2h}),\quad 0\le j\le \cfrac{n}{2}-1.\\
  \end{array}

.. figure:: ../images/multigrid9.png
   :width: 800
   :align: center
   
   1d case
  
.. math::
  \begin{array}{l}
  v_{0}^{h}=v_{0}^{2h}\\
  v_{1}^{h}=\cfrac{1}{2}( v_{0}^{2h}+v_{1}^{2h})\\
  v_{2}^{h}=v_{1}^{2h}\\
  v_{3}^{h}=\cfrac{1}{2}( v_{1}^{2h}+v_{2}^{2h})\\
  v_{4}^{h}=v_{2}^{2h}\\
  \end{array} 
  
For two-dimensional problems, the interpolation operator may be defined in a
similar way. If we let  :math:`I_{2h}^{h}\mathbf{v}^{2h}=\mathbf{v}^{h}`, then the components of :math:`\mathbf{v}^{h}` are given by  

.. math::
  \begin{array}{l}
  v_{2i,2j}^{h}=v_{i,j}^{2h}\\
  v_{2i+1,2j}^{h}=\cfrac{1}{2}( v_{i,j}^{2h}+v_{i+1,j}^{2h})\\
  v_{2i,2j+1}^{h}=\cfrac{1}{2}( v_{i,j}^{2h}+v_{i,j+1}^{2h})\\
  v_{2i+1,2j+1}^{h}=\cfrac{1}{4}( v_{i,j}^{2h}+v_{i+1,j}^{2h}+v_{i,j+1}^{2h}+v_{i+1,j+1}^{2h})\quad 0\le j\le \cfrac{n}{2}-1.\\
  \end{array}
  
.. figure:: ../images/multigrid10.png
   :width: 800
   :align: center
   
   2d case  
   
.. math::
  \begin{array}{l}
  v_{2i,2j}^{h}=v_{i,j}^{2h}\\
  v_{2i+1,2j}^{h}=\cfrac{1}{2}( v_{i,j}^{2h}+v_{i+1,j}^{2h})\\
  v_{2i,2j+1}^{h}=\cfrac{1}{2}( v_{i,j}^{2h}+v_{i,j+1}^{2h})\\
  v_{2i+1,2j+1}^{h}=\cfrac{1}{4}( v_{i,j}^{2h}+v_{i+1,j}^{2h}+v_{i,j+1}^{2h}+v_{i+1,j+1}^{2h})\quad 0\le j\le \cfrac{n}{2}-1.\\
  \end{array}      

-
   
.. math::
  \begin{array}{l}
  v_{0,0}^{h}=v_{0,0}^{2h}\\
  v_{1,0}^{h}=\cfrac{1}{2}( v_{0,0}^{2h}+v_{1,0}^{2h})\\
  v_{0,1}^{h}=\cfrac{1}{2}( v_{0,0}^{2h}+v_{0,1}^{2h})\\
  v_{1,1}^{h}=\cfrac{1}{4}( v_{0,0}^{2h}+v_{1,0}^{2h}+v_{0,1}^{2h}+v_{1,1}^{2h})\\
  \end{array}   

-
    
.. math::
  \begin{array}{l}
  v_{2,0}^{h}=v_{1,0}^{2h}\\
  v_{3,0}^{h}=\cfrac{1}{2}( v_{1,0}^{2h}+v_{2,0}^{2h})\\
  v_{2,1}^{h}=\cfrac{1}{2}( v_{1,0}^{2h}+v_{1,1}^{2h})\\
  v_{3,1}^{h}=\cfrac{1}{4}( v_{1,0}^{2h}+v_{2,0}^{2h}+v_{1,1}^{2h}+v_{2,1}^{2h})\\
  \end{array}
  
-
    
.. math::
  \begin{array}{l}
  v_{4,0}^{h}=v_{2,0}^{2h}\\
  v_{5,0}^{h}=\cfrac{1}{2}( v_{2,0}^{2h}+v_{3,0}^{2h})\\
  v_{4,1}^{h}=\cfrac{1}{2}( v_{2,0}^{2h}+v_{2,1}^{2h})\\
  v_{5,1}^{h}=\cfrac{1}{4}( v_{2,0}^{2h}+v_{3,0}^{2h}+v_{2,1}^{2h}+v_{3,1}^{2h})\\
  \end{array}  
  
-
    
.. math::
  \begin{array}{l}
  v_{0,2}^{h}=v_{0,j}^{2h}\\
  v_{1,2}^{h}=\cfrac{1}{2}( v_{0,1}^{2h}+v_{1,1}^{2h})\\
  v_{0,3}^{h}=\cfrac{1}{2}( v_{0,1}^{2h}+v_{0,2}^{2h})\\
  v_{1,3}^{h}=\cfrac{1}{4}( v_{0,1}^{2h}+v_{1,1}^{2h}+v_{0,2}^{2h}+v_{1,2}^{2h})\\
  \end{array}  
  
-
    
.. math::
  \begin{array}{l}
  v_{2,2}^{h}=v_{1,1}^{2h}\\
  v_{3,2}^{h}=\cfrac{1}{2}( v_{1,1}^{2h}+v_{2,1}^{2h})\\
  v_{2,3}^{h}=\cfrac{1}{2}( v_{1,1}^{2h}+v_{1,2}^{2h})\\
  v_{3,3}^{h}=\cfrac{1}{4}( v_{1,1}^{2h}+v_{2,1}^{2h}+v_{1,2}^{2h}+v_{2,2}^{2h})\\
  \end{array}  
  
-
    
.. math::
  \begin{array}{l}
  v_{4,2}^{h}=v_{2,1}^{2h}\\
  v_{5,2}^{h}=\cfrac{1}{2}( v_{2,1}^{2h}+v_{3,1}^{2h})\\
  v_{4,3}^{h}=\cfrac{1}{2}( v_{2,1}^{2h}+v_{2,2}^{2h})\\
  v_{5,3}^{h}=\cfrac{1}{4}( v_{2,1}^{2h}+v_{3,1}^{2h}+v_{2,2}^{2h}+v_{3,2}^{2h})\\
  \end{array}  
  
-
    
.. math::
  \begin{array}{l}
  v_{0,4}^{h}=v_{0,2}^{2h}\\
  v_{1,4}^{h}=\cfrac{1}{2}( v_{0,2}^{2h}+v_{1,2}^{2h})\\
  v_{0,5}^{h}=\cfrac{1}{2}( v_{0,2}^{2h}+v_{0,3}^{2h})\\
  v_{1,5}^{h}=\cfrac{1}{4}( v_{0,2}^{2h}+v_{1,2}^{2h}+v_{0,3}^{2h}+v_{1,3}^{2h})\\
  \end{array}  
  
-
    
.. math::
  \begin{array}{l}
  v_{2,4}^{h}=v_{1,2}^{2h}\\
  v_{3,4}^{h}=\cfrac{1}{2}( v_{1,2}^{2h}+v_{2,2}^{2h})\\
  v_{2,5}^{h}=\cfrac{1}{2}( v_{1,2}^{2h}+v_{1,3}^{2h})\\
  v_{3,5}^{h}=\cfrac{1}{4}( v_{1,2}^{2h}+v_{2,2}^{2h}+v_{1,3}^{2h}+v_{2,3}^{2h})\\
  \end{array}  
  
-
    
.. math::
  \begin{array}{l}
  v_{4,4}^{h}=v_{2,2}^{2h}\\
  v_{5,4}^{h}=\cfrac{1}{2}( v_{2,2}^{2h}+v_{3,2}^{2h})\\
  v_{4,5}^{h}=\cfrac{1}{2}( v_{2,2}^{2h}+v_{2,3}^{2h})\\
  v_{5,5}^{h}=\cfrac{1}{4}( v_{2,2}^{2h}+v_{3,2}^{2h}+v_{2,3}^{2h}+v_{3,3}^{2h})\\
  \end{array}  
  
Left, Right Boundary

.. math::
  \begin{array}{l}
  v_{2i,2j}^{h}=v_{i,j}^{2h}\\
  v_{2i,2j+1}^{h}=\cfrac{1}{2}( v_{i,j}^{2h}+v_{i,j+1}^{2h})\\
  \end{array}
  
Bottom, Top Boundary  

.. math::
  \begin{array}{l}
  v_{2i,2j}^{h}=v_{i,j}^{2h}\\
  v_{2i+1,2j}^{h}=\cfrac{1}{2}( v_{i,j}^{2h}+v_{i+1,j}^{2h})\\
  \end{array}
  
Left Boundary  

.. math::
  \begin{array}{l}
  v_{0,2j}^{h}=v_{0,j}^{2h}\\
  v_{0,2j+1}^{h}=\cfrac{1}{2}( v_{0,j}^{2h}+v_{0,j+1}^{2h})\\
  \end{array}
  
Right Boundary

.. math::
  \begin{array}{l}
  v_{2*nxc,2j}^{h}=v_{nxc,j}^{2h}\\
  v_{2*nxc,2j+1}^{h}=\cfrac{1}{2}( v_{nxc,j}^{2h}+v_{nxc,j+1}^{2h})\\
  \end{array}  

-

.. math::
  \begin{array}{l}
  v_{nxf,2j}^{h}=v_{nxc,j}^{2h}\\
  v_{nxf,2j+1}^{h}=\cfrac{1}{2}( v_{nxc,j}^{2h}+v_{nxc,j+1}^{2h})\\
  \end{array}    
  
Bottom Boundary  

.. math::
  \begin{array}{l}
  v_{2i,0}^{h}=v_{i,0}^{2h}\\
  v_{2i+1,0}^{h}=\cfrac{1}{2}( v_{i,0}^{2h}+v_{i+1,0}^{2h})\\
  \end{array}  
  
Top Boundary  

.. math::
  \begin{array}{l}
  v_{2i,2*nyc}^{h}=v_{i,nyc}^{2h}\\
  v_{2i+1,2*nyc}^{h}=\cfrac{1}{2}( v_{i,nyc}^{2h}+v_{i+1,nyc}^{2h})\\
  \end{array}  
  
-
  
.. math::
  \begin{array}{l}
  v_{2i,nyf}^{h}=v_{i,nyc}^{2h}\\
  v_{2i+1,nyf}^{h}=\cfrac{1}{2}( v_{i,nyc}^{2h}+v_{i+1,nyc}^{2h})\\
  \end{array}    
  
Left Boundary  

.. math::
  \begin{array}{l}
  v_{0,2j}^{h}=v_{0,j}^{2h},\quad j=0,1,\cdots,nyc;\quad 2j=0,2,\cdots,nyf;\\
  v_{0,2j+1}^{h}=\cfrac{1}{2}( v_{0,j}^{2h}+v_{0,j+1}^{2h}),\\
  j=0,1,\cdots,nyc-1;\\
  j+1=1,2,\cdots,nyc;\\
  \{j,j+1\}=\{0,1\},\{1,2\},\cdots,\{nyc-1,nyc\};\\
  2j+1=1,3,\cdots,2nyc-1=1,3,\cdots,nyf-1;
  \end{array}
  
Right Boundary  

.. math::
  \begin{array}{l}
  v_{nxf,2j}^{h}=v_{nxc,j}^{2h},\quad j=0,1,\cdots,nyc;\quad 2j=0,2,\cdots,nyf;\\
  v_{nxf,2j+1}^{h}=\cfrac{1}{2}( v_{nxc,j}^{2h}+v_{nxc,j+1}^{2h}),\\
  j=0,1,\cdots,nyc-1;\\
  j+1=1,2,\cdots,nyc;\\
  \{j,j+1\}=\{0,1\},\{1,2\},\cdots,\{nyc-1,nyc\};\\
  2j+1=1,3,\cdots,2nyc-1=1,3,\cdots,nyf-1;
  \end{array}
  
Bottom Boundary  

.. math::
  \begin{array}{l}
  v_{2i,0}^{h}=v_{i,0}^{2h},\quad i=0,1,\cdots,nxc;\quad 2i=0,2,\cdots,nxf;\\
  v_{2i+1,0}^{h}=\cfrac{1}{2}( v_{i,0}^{2h}+v_{i+1,0}^{2h})\\
  i=0,1,\cdots,nxc-1;\\
  i+1=1,2,\cdots,nxc;\\
  \{i,i+1\}=\{0,1\},\{1,2\},\cdots,\{nxc-1,nxc\};\\
  2i+1=1,3,\cdots,2nxc-1=1,3,\cdots,nxf-1;
  \end{array}
  
Top Boundary  

.. math::
  \begin{array}{l}
  v_{2i,nyf}^{h}=v_{i,nyc}^{2h},\quad i=0,1,\cdots,nxc;\quad 2i=0,2,\cdots,nxf;\\
  v_{2i+1,nyf}^{h}=\cfrac{1}{2}( v_{i,nyc}^{2h}+v_{i+1,nyc}^{2h})\\
  i=0,1,\cdots,nxc-1;\\
  i+1=1,2,\cdots,nxc;\\
  \{i,i+1\}=\{0,1\},\{1,2\},\cdots,\{nxc-1,nxc\};\\
  2i+1=1,3,\cdots,2nxc-1=1,3,\cdots,nxf-1;
  \end{array}