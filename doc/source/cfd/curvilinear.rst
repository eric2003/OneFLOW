Generalized Curvilinear Coordinate System
===========================================

- The Geometric Conservation Law - A link between finite-difference and finite-volume methods of flow computation on moving grids
- Development of CFD Algorithms for Transient and Steady Aerodynamics- Fergal Boyle [Thesis] 

two dimension static grid
----------------------------

.. math::
  \begin{array}{c}
  x=x(\xi,\eta)=\phi_{1}(\xi,\eta)\\
  y=y(\xi,\eta)=\phi_{2}(\xi,\eta)\\
  \end{array}

-  
  
.. math:: 
  \begin{array}{c}
  \xi=\xi(x,y)=\psi_{1}(x,y)\\
  \eta=\eta(x,y)=\psi_{2}(x,y)\\
  \end{array}  
  
-  
  
.. math::  
  \begin{array}{c}
  \cfrac{\partial f}{\partial x}=\cfrac{\partial f}{\partial \xi}\cfrac{\partial \xi}{\partial x}+
  \cfrac{\partial f}{\partial \eta}\cfrac{\partial \eta}{\partial x}\\
  \cfrac{\partial f}{\partial y}=\cfrac{\partial f}{\partial \xi}\cfrac{\partial \xi}{\partial y}+
  \cfrac{\partial f}{\partial \eta}\cfrac{\partial \eta}{\partial y}\\
  \end{array}  

-  
  
.. math::  
  \begin{array}{c}
  dx=x_{\xi}d\xi+x_{\eta}d\eta=({\phi_{1}})_{\xi}d\xi+({\phi_{1}})_{\eta}d\eta\\
  dy=y_{\xi}d\xi+y_{\eta}d\eta=({\phi_{2}})_{\xi}d\xi+({\phi_{2}})_{\eta}d\eta\\
  \end{array}  
  
-  
  
.. math::  
  \begin{array}{c}
  \begin{bmatrix}
   dx\\dy
  \end{bmatrix}
  =\begin{bmatrix}
    x_{\xi}& x_{\eta}\\
    y_{\xi}& y_{\eta}\\
  \end{bmatrix}\begin{bmatrix}
   d\xi\\d\eta
  \end{bmatrix}
  \end{array} 
  
-  
  
.. math::  
  \begin{array}{c}
  d{\xi}={\xi}_{x}dx+{\xi}_{y}dy=(\psi_{1})_{x}dx+(\psi_{1})_{y}dy\\
  d{\eta}={\eta}_{x}dx+{\eta}_{y}dy=(\psi_{2})_{x}dx+(\psi_{2})_{y}dy\\
  \end{array}  
  
-  
  
.. math:: 
  \begin{bmatrix}
   d\xi\\d\eta
  \end{bmatrix}
  =\begin{bmatrix}
   {\xi}_{x}& {\xi}_{y}\\
   {\eta}_{x}& {\eta}_{y}\\
  \end{bmatrix}\begin{bmatrix}
   dx\\dy
  \end{bmatrix}
  
-  
  
.. math:: 
  \begin{bmatrix}
   d\xi\\d\eta
  \end{bmatrix}
  =\begin{bmatrix}
   {\xi}_{x}& {\xi}_{y}\\
   {\eta}_{x}& {\eta}_{y}\\
  \end{bmatrix}\begin{bmatrix}
   dx\\dy
  \end{bmatrix}
  = \begin{bmatrix}
   {\xi}_{x}& {\xi}_{y}\\
   {\eta}_{x}& {\eta}_{y}\\
  \end{bmatrix}
  \begin{bmatrix}
    x_{\xi}& x_{\eta}\\
    y_{\xi}& y_{\eta}\\
  \end{bmatrix}\begin{bmatrix}
   d\xi\\d\eta
  \end{bmatrix} 
  
-  
  
.. math:: 
  \begin{bmatrix}
   {\xi}_{x}& {\xi}_{y}\\
   {\eta}_{x}& {\eta}_{y}\\
  \end{bmatrix}
  \begin{bmatrix}
    x_{\xi}& x_{\eta}\\
    y_{\xi}& y_{\eta}\\
  \end{bmatrix}=I  
  
-  
  
.. math:: 
  \begin{bmatrix}
   {\xi}_{x}& {\xi}_{y}\\
   {\eta}_{x}& {\eta}_{y}\\
  \end{bmatrix}=
  \begin{bmatrix}
    x_{\xi}& x_{\eta}\\
    y_{\xi}& y_{\eta}\\
  \end{bmatrix}^{-1}
  
-  
  
.. math:: 
  \begin{bmatrix}
    x_{\xi}& x_{\eta}\\
    y_{\xi}& y_{\eta}\\
  \end{bmatrix}=\begin{bmatrix}
  {\xi}_{x}& {\xi}_{y}\\
  {\eta}_{x}& {\eta}_{y}\\
  \end{bmatrix}^{-1} 
  
inverse of matrix

.. math:: 
  \mathbf{A}^{-1}=\cfrac{1}{\text{det}(\mathbf{A})}\text{adj}(\mathbf{A})=\cfrac{1}{|\mathbf{A}|}\mathbf{A}^{*}  
  
For example,find the inverse matrix of the second-order matrix :math:`\mathbf{A}=\begin{pmatrix}a& b\\c&d\end{pmatrix}`

.. math:: 
   |\mathbf{A}|=ad-bc
   
-
   
.. math:: 
   \mathbf{A}^{*}=\begin{pmatrix}d& -b\\-c&a\end{pmatrix}\\
   
-
   
.. math:: 
  \mathbf{A}^{-1}=\cfrac{1}{|\mathbf{A}|}\mathbf{A}^{*} =\cfrac{1}{ad-bc}\begin{pmatrix}d& -b\\-c&a\end{pmatrix}   
  
Let

.. math:: 
  \mathbf{A}=\begin{bmatrix}x_{\xi}& x_{\eta}\\y_{\xi}& y_{\eta}\\\end{bmatrix}  
  
then

.. math::
  |\mathbf{A}|=x_{\xi}y_{\eta}-x_{\eta}y_{\xi}\\

-

.. math:: 
  \mathbf{A}^{*}=\begin{pmatrix}y_{\eta}& -x_{\eta}\\-y_{\xi}&x_{\xi}\end{pmatrix}\\
  
-

.. math::   
  \mathbf{A}^{-1}=\cfrac{1}{|\mathbf{A}|}\mathbf{A}^{*} =\cfrac{1}{x_{\xi}y_{\eta}-x_{\eta}y_{\xi}}\begin{pmatrix}y_{\eta}& -x_{\eta}\\-y_{\xi}&x_{\xi}\end{pmatrix}\\  
  
Let  

.. math:: 
  J=\begin{vmatrix}x_{\xi}& x_{\eta}\\y_{\xi}& y_{\eta}\\\end{vmatrix} 
  
then  

.. math:: 
  \begin{bmatrix}{\xi}_{x}& {\xi}_{y}\\{\eta}_{x}& {\eta}_{y}\\\end{bmatrix} =\cfrac{1}{J}\begin{bmatrix}y_{\eta}& -x_{\eta}\\-y_{\xi}&x_{\xi}\end{bmatrix}\\

-

.. math::     
  \begin{array}{c}
  {\xi}_{x}=\cfrac{1}{J}(y_{\eta})\quad{\xi}_{y}=\cfrac{1}{J}(-x_{\eta})\\
  {\eta}_{x}=\cfrac{1}{J}(-y_{\xi})\quad{\eta}_{y}=\cfrac{1}{J}(x_{\xi})\\
  \end{array}

Equations in Cartesian Coordinates

.. math::       
  \cfrac{\partial u}{\partial t}+\cfrac{\partial E}{\partial x}+\cfrac{\partial F}{\partial y}=0  

-

.. math::
  \begin{array}{c}
  \cfrac{\partial E}{\partial x}=\cfrac{\partial E}{\partial \xi}\cfrac{\partial \xi}{\partial x}+\cfrac{\partial E}{\partial \eta}\cfrac{\partial \eta}{\partial x}\\
  \cfrac{\partial F}{\partial y}=\cfrac{\partial F}{\partial \xi}\cfrac{\partial \xi}{\partial y}+\cfrac{\partial F}{\partial \eta}\cfrac{\partial \eta}{\partial y}\\
  \end{array}
  
-

.. math::
  \cfrac{\partial u}{\partial t}+\cfrac{\partial E}{\partial \xi}\cfrac{\partial \xi}{\partial x}+\cfrac{\partial E}{\partial \eta}\cfrac{\partial \eta}{\partial x}
  +\cfrac{\partial F}{\partial \xi}\cfrac{\partial \xi}{\partial y}+\cfrac{\partial F}{\partial \eta}\cfrac{\partial \eta}{\partial y}=0
  
-

.. math::
  \cfrac{\partial u}{\partial t}+\cfrac{\partial E}{\partial \xi}\cfrac{\partial \xi}{\partial x}+\cfrac{\partial F}{\partial \xi}\cfrac{\partial \xi}{\partial y}+\cfrac{\partial E}{\partial \eta}\cfrac{\partial \eta}{\partial x}
  +\cfrac{\partial F}{\partial \eta}\cfrac{\partial \eta}{\partial y}=0 
  
-

.. math::
  J\cfrac{\partial u}{\partial t}+J\cfrac{\partial E}{\partial \xi}\cfrac{\partial \xi}{\partial x}+J\cfrac{\partial F}{\partial \xi}\cfrac{\partial \xi}{\partial y}+J\cfrac{\partial E}{\partial \eta}\cfrac{\partial \eta}{\partial x}
  +J\cfrac{\partial F}{\partial \eta}\cfrac{\partial \eta}{\partial y}=0  
  
-

.. math::
  \cfrac{\partial (EJ\cfrac{\partial \xi}{\partial x})}{\partial \xi}=\cfrac{\partial E}{\partial \xi}(J\cfrac{\partial \xi}{\partial x})+E\cfrac{\partial (J\cfrac{\partial \xi}{\partial x})}{\partial \xi}  
  
-

.. math::
  \cfrac{\partial (EJ\cfrac{\partial \eta}{\partial x})}{\partial \eta}=\cfrac{\partial E}{\partial \eta}(J\cfrac{\partial \eta}{\partial x})+E\cfrac{\partial (J\cfrac{\partial \eta}{\partial x})}{\partial \eta}\\
  
-

.. math::  
  \cfrac{\partial (FJ\cfrac{\partial \xi}{\partial y})}{\partial \xi}=\cfrac{\partial F}{\partial \xi}(J\cfrac{\partial \xi}{\partial y})+F\cfrac{\partial (J\cfrac{\partial \xi}{\partial y})}{\partial \xi}\\  
  
-

.. math::  
  \cfrac{\partial (FJ\cfrac{\partial \eta}{\partial y})}{\partial \eta}=\cfrac{\partial F}{\partial \eta}(J\cfrac{\partial \eta}{\partial y})+F\cfrac{\partial (J\cfrac{\partial \eta}{\partial y})}{\partial \eta}\\  
  
-

.. math:: 
  E\cfrac{\partial (J\cfrac{\partial \xi}{\partial x})}{\partial \xi}
  +E\cfrac{\partial (J\cfrac{\partial \eta}{\partial x})}{\partial \eta}
  =E(\cfrac{\partial (y_{\eta})}{\partial \xi}+\cfrac{\partial (-y_{\xi})}{\partial \eta})=0
  
-

.. math:: 
  F\cfrac{\partial (J\cfrac{\partial \xi}{\partial y})}{\partial \xi}
  +F\cfrac{\partial (J\cfrac{\partial \eta}{\partial y})}{\partial \eta}
  =F(\cfrac{\partial (-x_{\eta})}{\partial \xi}+\cfrac{\partial (x_{\xi})}{\partial \eta})=0

-

.. math:: 
  J\cfrac{\partial u}{\partial t}
  +\cfrac{\partial (EJ\cfrac{\partial \xi}{\partial x})}{\partial \xi}
  +\cfrac{\partial (FJ\cfrac{\partial \xi}{\partial y})}{\partial \xi}
  +\cfrac{\partial (EJ\cfrac{\partial \eta}{\partial x})}{\partial \eta}
  +\cfrac{\partial (FJ\cfrac{\partial \eta}{\partial y})}{\partial \eta}=0  
  
-

.. math:: 
  J\cfrac{\partial u}{\partial t}
  +\cfrac{\partial (EJ\cfrac{\partial \xi}{\partial x})}{\partial \xi}
  +\cfrac{\partial (FJ\cfrac{\partial \xi}{\partial y})}{\partial \xi}
  +\cfrac{\partial (EJ\cfrac{\partial \eta}{\partial x})}{\partial \eta}
  +\cfrac{\partial (FJ\cfrac{\partial \eta}{\partial y})}{\partial \eta}=0

-

.. math:: 
  J\cfrac{\partial u}{\partial t}
  +\cfrac{\partial (J(\xi_{x}E+\xi_{y}F))}{\partial \xi}
  +\cfrac{\partial (J(\eta_{x}E+\eta_{y}F))}{\partial \eta}
  =0
  
two dimension dynamic grid
---------------------------- 

.. math::
  \begin{align}
  x & = x(\xi,\eta,\tau)  = \phi_{1}(\xi,\eta,\tau)\\
  y & = y(\xi,\eta,\tau)  = \phi_{2}(\xi,\eta,\tau)\\
  t & = \tau
  \end{align}
  
-

.. math::
  \begin{align}
  \xi & = \xi(x,y,t)  = \psi_{1}(x,y,t)\\
  \eta & = \eta(x,y,t)  = \psi_{2}(x,y,t)\\
  \tau&=t
  \end{align} 

-

.. math::
  J=\cfrac{\partial (x,y)}{\partial (\xi,\eta)}=\begin{vmatrix}x_{\xi}& x_{\eta}\\y_{\xi}& y_{\eta}\\\end{vmatrix}   
  
Equations in general curvilinear coordinate system 

.. math::
  J(y_{1},y_{2},y_{3})=\cfrac{\partial (y_{1},y_{2},y_{3})}{\partial (\xi_{1},\xi_{2},\xi_{3})}
  =\begin{vmatrix}
  \cfrac{\partial y_{1}}{\partial \xi_{1}}&
  \cfrac{\partial y_{1}}{\partial \xi_{2}}&
  \cfrac{\partial y_{1}}{\partial \xi_{3}} \\
  \cfrac{\partial y_{2}}{\partial \xi_{1}}&
  \cfrac{\partial y_{2}}{\partial \xi_{2}}&
  \cfrac{\partial y_{2}}{\partial \xi_{3}} \\
  \cfrac{\partial y_{3}}{\partial \xi_{1}}&
  \cfrac{\partial y_{3}}{\partial \xi_{2}}&
  \cfrac{\partial y_{3}}{\partial \xi_{3}} \\
  \end{vmatrix}
 
-

.. math::
  \cfrac{\partial (y_{1},y_{2},\cdots,y_{i},y_{i+1},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{n})}
  =-\cfrac{\partial (y_{1},y_{2},\cdots,y_{i+1},y_{i},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{n})} 
  
-

.. math::
  \cfrac{\partial (y_{1},y_{2},\cdots,y_{i},y_{i+1},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{n})}
  =\begin{vmatrix}
  \cfrac{\partial y_{1}}{\partial \xi_{1}}& \cfrac{\partial y_{1}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{1}}{\partial \xi_{n}}\\
  \cfrac{\partial y_{2}}{\partial \xi_{1}}& \cfrac{\partial y_{2}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{2}}{\partial \xi_{n}}\\
  \vdots & \vdots&\ddots   & \vdots\\
  \cfrac{\partial y_{i}}{\partial \xi_{1}}& \cfrac{\partial y_{i}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{i}}{\partial \xi_{n}}\\
  \cfrac{\partial y_{i+1}}{\partial \xi_{1}}& \cfrac{\partial y_{i+1}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{i+1}}{\partial \xi_{n}}\\
  \vdots & \vdots&\ddots   & \vdots\\
  \cfrac{\partial y_{n}}{\partial \xi_{1}}& \cfrac{\partial y_{n}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{n}}{\partial \xi_{n}}\\
  \end{vmatrix}  

-

.. math::  
  \cfrac{\partial (y_{1},y_{2},\cdots,y_{i+1},y_{i},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{n})}
  =\begin{vmatrix}
  \cfrac{\partial y_{1}}{\partial \xi_{1}}& \cfrac{\partial y_{1}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{1}}{\partial \xi_{n}}\\
  \cfrac{\partial y_{2}}{\partial \xi_{1}}& \cfrac{\partial y_{2}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{2}}{\partial \xi_{n}}\\
  \vdots & \vdots&\ddots   & \vdots\\
  \cfrac{\partial y_{i+1}}{\partial \xi_{1}}& \cfrac{\partial y_{i+1}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{i+1}}{\partial \xi_{n}}\\
  \cfrac{\partial y_{i}}{\partial \xi_{1}}& \cfrac{\partial y_{i}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{i}}{\partial \xi_{n}}\\
  \vdots & \vdots&\ddots   & \vdots\\
  \cfrac{\partial y_{n}}{\partial \xi_{1}}& \cfrac{\partial y_{n}}{\partial \xi_{2}}&\cdots  & \cfrac{\partial y_{n}}{\partial \xi_{n}}\\
  \end{vmatrix}  
  
-

.. math::  
  \cfrac{\partial (y_{1},y_{2},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{i},\xi_{i+1},\cdots,\xi_{n})}
  =-\cfrac{\partial (y_{1},y_{2},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{i+1},\xi_{i},\cdots,\xi_{n})}  
  
-

.. math::  
  \cfrac{\partial (y_{1},y_{2},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{i},\xi_{i+1},\cdots,\xi_{n})}
  =\begin{vmatrix}
  \cfrac{\partial y_{1}}{\partial \xi_{1}}&
  \cfrac{\partial y_{1}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{1}}{\partial \xi_{i}}&
  \cfrac{\partial y_{1}}{\partial \xi_{i+1}}&
  \cdots&
  \cfrac{\partial y_{1}}{\partial \xi_{n}}\\
  \cfrac{\partial y_{2}}{\partial \xi_{1}}&
  \cfrac{\partial y_{2}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{2}}{\partial \xi_{i}}&
  \cfrac{\partial y_{2}}{\partial \xi_{i+1}}&
  \cdots&
  \cfrac{\partial y_{2}}{\partial \xi_{n}}\\ 
  \vdots & \vdots&\ddots   & \vdots& \vdots&\ddots& \vdots\\
  \cfrac{\partial y_{n}}{\partial \xi_{1}}&
  \cfrac{\partial y_{n}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{n}}{\partial \xi_{i}}&
  \cfrac{\partial y_{n}}{\partial \xi_{i+1}}&
  \cdots&
  \cfrac{\partial y_{n}}{\partial \xi_{n}}\\
  \end{vmatrix} 

-

.. math::  
  \cfrac{\partial (y_{1},y_{2},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{i+1},\xi_{i},\cdots,\xi_{n})}
  =\begin{vmatrix}
  \cfrac{\partial y_{1}}{\partial \xi_{1}}&
  \cfrac{\partial y_{1}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{1}}{\partial \xi_{i+1}}&
  \cfrac{\partial y_{1}}{\partial \xi_{i}}&
  \cdots&
  \cfrac{\partial y_{1}}{\partial \xi_{n}}\\
  \cfrac{\partial y_{2}}{\partial \xi_{1}}&
  \cfrac{\partial y_{2}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{2}}{\partial \xi_{i+1}}&
  \cfrac{\partial y_{2}}{\partial \xi_{i}}&
  \cdots&
  \cfrac{\partial y_{2}}{\partial \xi_{n}}\\ 
  \vdots & \vdots&\ddots   & \vdots& \vdots&\ddots& \vdots\\
  \cfrac{\partial y_{n}}{\partial \xi_{1}}&
  \cfrac{\partial y_{n}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{n}}{\partial \xi_{i+1}}&
  \cfrac{\partial y_{n}}{\partial \xi_{i}}&
  \cdots&
  \cfrac{\partial y_{n}}{\partial \xi_{n}}\\
  \end{vmatrix} 
  
-

.. math:: 
  \cfrac{\partial (y_{1},y_{2},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{n})}
  =\begin{vmatrix}
  \cfrac{\partial y_{1}}{\partial \xi_{1}}&
  \cfrac{\partial y_{1}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{1}}{\partial \xi_{n}}\\
  \cfrac{\partial y_{2}}{\partial \xi_{1}}&
  \cfrac{\partial y_{2}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{2}}{\partial \xi_{n}}\\
  \vdots & \vdots&\ddots   & \vdots\\
  \cfrac{\partial y_{n}}{\partial \xi_{1}}&
  \cfrac{\partial y_{n}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{n}}{\partial \xi_{n}}\\
  \end{vmatrix}  
  
When there are common variables between :math:`y_{i}` and :math:`\xi_{i}`, then determinant reduction occurs, such as  

-

.. math::  
  \cfrac{\partial (\xi_{1},y_{2},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{n})}
  =\begin{vmatrix}
  \cfrac{\partial \xi_{1}}{\partial \xi_{1}}&
  \cfrac{\partial \xi_{1}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial \xi_{1}}{\partial \xi_{n}}\\
  \cfrac{\partial y_{2}}{\partial \xi_{1}}&
  \cfrac{\partial y_{2}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{2}}{\partial \xi_{n}}\\
  \vdots & \vdots&\ddots   & \vdots\\
  \cfrac{\partial y_{n}}{\partial \xi_{1}}&
  \cfrac{\partial y_{n}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{n}}{\partial \xi_{n}}\\
  \end{vmatrix}=
  \begin{vmatrix}
  1&
  0&
  \cdots&
  0\\
  \cfrac{\partial y_{2}}{\partial \xi_{1}}&
  \cfrac{\partial y_{2}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{2}}{\partial \xi_{n}}\\
  \vdots & \vdots&\ddots   & \vdots\\
  \cfrac{\partial y_{n}}{\partial \xi_{1}}&
  \cfrac{\partial y_{n}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{n}}{\partial \xi_{n}}\\
  \end{vmatrix}
  
-

.. math:: 
  \cfrac{\partial (\xi_{1},y_{2},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{n})}
  =\cfrac{\partial (y_{2},\cdots,y_{n})}{\partial (\xi_{2},\cdots,\xi_{n})}=
  \begin{vmatrix}
  \cfrac{\partial y_{2}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{2}}{\partial \xi_{n}}\\
  \vdots & \ddots   & \vdots\\
  \cfrac{\partial y_{n}}{\partial \xi_{2}}&
  \cdots&
  \cfrac{\partial y_{n}}{\partial \xi_{n}}\\
  \end{vmatrix}
  
-

.. math:: 
  \cfrac{\partial (y_{1},y_{2},\cdots,y_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{n})}
  =\cfrac{\partial (y_{1},y_{2},\cdots,y_{n})/\partial (\zeta_{1},\zeta_{2},\cdots,\zeta_{n})}{\partial (\xi_{1},\xi_{2},\cdots,\xi_{n})/\partial (\zeta_{1},\zeta_{2},\cdots,\zeta_{n})}=\cfrac{J(y_{1},y_{2},\cdots,y_{n})}{J(\xi_{1},\xi_{2},\cdots,\xi_{n})}  
  
-

.. math:: 
  \cfrac{\partial U}{\partial t}+\cfrac{\partial E}{\partial x}+\cfrac{\partial F}{\partial y}=0 
  
-

.. math:: 
  \begin{array}{c}
  q_{0}=U,q_{1}=E,q_{2}=F\\
  x_{0}=t,x_{1}=x,x_{2}=y\\
  \xi_{0}=t,\xi_{1}=\xi,\xi_{2}=\eta\\
  \end{array}  
  
-

.. math:: 
  \cfrac{\partial q_{i}}{\partial x_{i}}=0 
  
-

.. math:: 
  Q_{k}=Jq_{i}\cfrac{\partial \xi_{k}}{\partial x_{i}}
  
-

.. math:: 
  \begin{align}
  Q_{k} & = Jq_{0}\cfrac{\partial \xi_{k}}{\partial x_{0}}
     +Jq_{1}\cfrac{\partial \xi_{k}}{\partial x_{1}}
     +Jq_{2}\cfrac{\partial \xi_{k}}{\partial x_{2}}\\
  & = Jq_{0}\cfrac{\partial \xi_{k}}{\partial t}
     +Jq_{1}\cfrac{\partial \xi_{k}}{\partial x}
     +Jq_{2}\cfrac{\partial \xi_{k}}{\partial y}
  \end{align}
  
-

.. math:: 
  \begin{align}
  Q_{0} 
  & = Jq_{0}\cfrac{\partial t}{\partial t}
     +Jq_{1}\cfrac{\partial t}{\partial x}
     +Jq_{2}\cfrac{\partial t}{\partial y}=Jq_{0}
  \end{align}
  
-

.. math:: 
  \begin{align}
  Q_{1}
  & = Jq_{0}\cfrac{\partial \xi}{\partial t}
     +Jq_{1}\cfrac{\partial \xi}{\partial x}
     +Jq_{2}\cfrac{\partial \xi}{\partial y}
  \end{align} 
  
-

.. math:: 
  \begin{align}
  Q_{2}
  & = Jq_{0}\cfrac{\partial \eta}{\partial t}
     +Jq_{1}\cfrac{\partial \eta}{\partial x}
     +Jq_{2}\cfrac{\partial \eta}{\partial y}
  \end{align} 
  
-

.. math:: 
  \cfrac{\partial U}{\partial t}
  =\cfrac{\partial U}{\partial \tau}\cfrac{\partial \tau}{\partial t}
  +\cfrac{\partial U}{\partial \xi}\cfrac{\partial \xi}{\partial t}
  +\cfrac{\partial U}{\partial \eta}\cfrac{\partial \eta}{\partial t}
  
-

.. math:: 
  \cfrac{\partial U(x,y,t)}{\partial t}
  =\cfrac{\partial \hat{U}(\xi,\eta,\tau)}{\partial \tau}\cfrac{\partial \tau}{\partial t}
  +\cfrac{\partial \hat{U}(\xi,\eta,\tau)}{\partial \xi}\cfrac{\partial \xi}{\partial t}
  +\cfrac{\partial \hat{U}(\xi,\eta,\tau)}{\partial \eta}\cfrac{\partial \eta}{\partial t}  
  
-

.. math:: 
  \begin{align}
  Q_{0} 
  & = JU\cfrac{\partial t}{\partial t}
     +JE\cfrac{\partial t}{\partial x}
     +JF\cfrac{\partial t}{\partial y}=JU
  \end{align} 
  
-

.. math:: 
  \begin{align}
  Q_{1}
  & = UJ\cfrac{\partial \xi}{\partial t}
     +EJ\cfrac{\partial \xi}{\partial x}
     +FJ\cfrac{\partial \xi}{\partial y}
  \end{align}

-

.. math:: 
  \begin{align}
  Q_{2}
  & = UJ\cfrac{\partial \eta}{\partial t}
     +EJ\cfrac{\partial \eta}{\partial x}
     +FJ\cfrac{\partial \eta}{\partial y}
  \end{align} 

-

.. math::   
  \cfrac{\partial (UJ\cfrac{\partial \xi}{\partial t})}{\partial \xi}=
  \cfrac{\partial (U)}{\partial \xi}(J\cfrac{\partial \xi}{\partial t})+
  (U)\cfrac{\partial }{\partial \xi}(J\cfrac{\partial \xi}{\partial t})  
  
-

.. math::  
  \begin{align}
  x & = x(\xi,\eta,\tau)  = \phi_{1}(\xi,\eta,\tau)\\
  y & = y(\xi,\eta,\tau)  = \phi_{2}(\xi,\eta,\tau)\\
  t & = t(\xi,\eta,\tau)=\phi_{0}(\xi,\eta,\tau)=\tau\\
  \end{align} 
  
-

.. math::
  \begin{align}
  \xi & = \xi(x,y,t)  = \psi_{1}(x,y,t)\\
  \eta & = \eta(x,y,t)  = \psi_{2}(x,y,t)\\
  \tau&=\tau(x,y,t)= \psi_{0}(x,y,t)=t
  \end{align}   