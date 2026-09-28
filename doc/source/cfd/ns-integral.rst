ALE Form of Conservation Equations
========================================

Differential forms
----------------------------

.. math::
  \begin{align}
  &Mass:     &\frac{\mathrm{d} \rho}{\mathrm{d} t}  &= \left.\frac{\partial \rho}{\partial t}\right|_{\boldsymbol\chi}+\mathbf{c} \cdot \nabla \rho  = -\rho {\nabla} \cdot \mathbf{v} \\
  &Momentum:  &\rho \frac{\mathrm{d} \mathbf{v}}{\mathrm{d} t}  &= \rho\left(\left.\frac{\partial \mathbf{v}}{\partial t}\right|_{\boldsymbol\chi}+(\mathbf{c} \cdot {\nabla}) \mathbf{v}\right)  = {\nabla} \cdot \boldsymbol{\sigma}+\rho \mathbf{b} \\
  &Energy:   &\rho \frac{\mathrm{d} E}{\mathrm{~d} t}  &= \rho\left(\left.\frac{\partial E}{\partial t}\right|_{\boldsymbol\chi}+\mathbf{c} \cdot \nabla E\right)  = \nabla \cdot(\boldsymbol{\sigma} \cdot \mathbf{v})-\nabla\cdot\mathbf{q}+\mathbf{v} \cdot \rho \mathbf{b} \\
  \end{align}
  
-  

.. math:: 
  \begin{align}
  Mass\quad\quad\quad:     &\frac{\partial \rho}{\partial t}\Bigg|_{\boldsymbol\chi}+{\nabla}\cdot(\rho\mathbf{v})-(\mathbf{\hat{v}}\cdot \nabla)( \rho) = 0\\
  Momentum:  & \cfrac{\partial (\rho\mathbf{v})}{\partial t}\Bigg|_{\boldsymbol{\chi}}  +\text{div }(\rho\mathbf{v}\otimes\mathbf{v}) -(\mathbf{\hat{v}} \cdot \nabla)(\rho\mathbf{v})= {\nabla} \cdot \boldsymbol{\sigma}+\rho \mathbf{b} \\
  Energy\quad\quad:   &\frac{\partial (\rho E)}{\partial t}\Bigg|_{\boldsymbol\chi}+{\nabla} \cdot (\rho E \mathbf{v})+(\mathbf{\hat{v}} \cdot \nabla) (\rho E) = \nabla \cdot(\boldsymbol{\sigma} \cdot \mathbf{v})-\nabla\cdot\mathbf{q}+\mathbf{v} \cdot \rho \mathbf{b} \\
  \end{align}
  
-  

.. math::
  \begin{align}
  Mass\quad\quad\quad&:     \frac{\partial \rho}{\partial t}\Bigg|_{\boldsymbol\chi}+{\nabla}\cdot(\rho\mathbf{c})+(\rho) ({\nabla} \cdot \mathbf{\hat{v}}) = 0  \\
  Momentum&:   \cfrac{\partial (\rho\mathbf{v})}{\partial t}\Bigg|_{\boldsymbol{\chi}}+\text{div }(\rho\mathbf{v}\otimes\mathbf{c}) +(\rho\mathbf{v})( {\nabla} \cdot \mathbf{\hat{v}})= {\nabla} \cdot \boldsymbol{\sigma}+\rho \mathbf{b} \\
  Energy\quad\quad&:    \frac{\partial (\rho E)}{\partial t}\Bigg|_{\boldsymbol\chi}+{\nabla} \cdot (\rho E \mathbf{c})+(\rho E )({\nabla} \cdot \mathbf{\hat{v}}) = \nabla \cdot(\boldsymbol{\sigma} \cdot \mathbf{v})-\nabla\cdot\mathbf{q}+\mathbf{v} \cdot \rho \mathbf{b} \\
  \end{align}
  
Integral forms
----------------------------

.. math:: 
  \begin{array}{l}
  \displaystyle \cfrac{\text{d}(m)_{sys}}{\text{d}t}=\cfrac{\text{d}}{\text{d}t}\int_{\text{V}_t}\rho \text{d}V
  +\int_{\text{S}_t}\rho\mathbf{c}\cdot\mathbf{n}\text{d}S=0\\
  \displaystyle \cfrac{\text{d}(mV)_{sys}}{\text{d}t}=\cfrac{\text{d}}{\text{d}t}\int_{\text{V}_t}(\rho\mathbf{v}) \text{d}V
  +\int_{\text{S}_t}(\rho\mathbf{v})\mathbf{c}\cdot\mathbf{n}\text{d}S
  =\int_{\text{S}_t}(\boldsymbol\sigma\cdot\mathbf{n})\text{d}S\\
  \displaystyle\cfrac{\text{d}(mE)_{sys}}{\text{d}t}=\cfrac{\text{d}}{\text{d}t}\int_{\text{V}_t}(\rho E) \text{d}V
  +\int_{\text{S}_t}(\rho E)\mathbf{c}\cdot\mathbf{n}\text{d}S
  =\int_{\text{S}_t}(\boldsymbol\sigma\cdot\mathbf{v}\cdot\mathbf{n}-\mathbf{q}\cdot\mathbf{n})\text{d}S
  \end{array}

where
  
.. math:: 
  \boldsymbol\sigma =\begin{bmatrix}
  \sigma_{11}& \sigma_{12} & \sigma_{13}\\
  \sigma_{21}& \sigma_{22} & \sigma_{23}\\
  \sigma_{31}& \sigma_{32} & \sigma_{33}\\
  \end{bmatrix}
  =\begin{bmatrix}
  \sigma_{xx}& \sigma_{xy} & \sigma_{xz}\\
  \sigma_{yx}& \sigma_{yy} & \sigma_{yz}\\
  \sigma_{zx}& \sigma_{zy} & \sigma_{zz}\\
  \end{bmatrix}  

-
  
.. math::   
  [\boldsymbol{\tau}\cdot \mathbf{n}]=
  \begin{bmatrix}
    {\tau}_{xx}& {\tau}_{xy} & {\tau}_{xz}\\
    {\tau}_{yx}& {\tau}_{yy} & {\tau}_{yz}\\
    {\tau}_{zx}& {\tau}_{zy} & {\tau}_{zz}\\
  \end{bmatrix}
  \begin{bmatrix}
    n_{x}\\ n_{y}\\  n_{z}\\
  \end{bmatrix}
  =\begin{bmatrix}
    {\tau}_{xx}n_{x}+{\tau}_{xy}n_{y}+{\tau}_{xz}n_{z}\\
    {\tau}_{yx}n_{x}+{\tau}_{yy}n_{y}+{\tau}_{yz}n_{z}\\
    {\tau}_{zx}n_{x}+{\tau}_{zy}n_{y}+{\tau}_{zz}n_{z}\\
  \end{bmatrix}  
  
-
  
.. math::    
  [\boldsymbol{\sigma}\cdot \mathbf{n}]=
  \begin{bmatrix}
    {\sigma}_{xx}& {\sigma}_{xy} & {\sigma}_{xz}\\
    {\sigma}_{yx}& {\sigma}_{yy} & {\sigma}_{yz}\\
    {\sigma}_{zx}& {\sigma}_{zy} & {\sigma}_{zz}\\
  \end{bmatrix}
  \begin{bmatrix}
    n_{x}\\ n_{y}\\  n_{z}\\
  \end{bmatrix}
  =\begin{bmatrix}
    {\sigma}_{xx}n_{x}+{\sigma}_{xy}n_{y}+{\sigma}_{xz}n_{z}\\
    {\sigma}_{yx}n_{x}+{\sigma}_{yy}n_{y}+{\sigma}_{yz}n_{z}\\
    {\sigma}_{zx}n_{x}+{\sigma}_{zy}n_{y}+{\sigma}_{zz}n_{z}\\
  \end{bmatrix}  
  
-
  
.. math:: 
  \begin{align}
  \begin{bmatrix}
    {\sigma}_{xx}n_{x}+{\sigma}_{xy}n_{y}+{\sigma}_{xz}n_{z}\\
    {\sigma}_{yx}n_{x}+{\sigma}_{yy}n_{y}+{\sigma}_{yz}n_{z}\\
    {\sigma}_{zx}n_{x}+{\sigma}_{zy}n_{y}+{\sigma}_{zz}n_{z}\\
  \end{bmatrix} & = \begin{bmatrix}
    ({\tau}_{xx}-p)n_{x}+{\tau}_{xy}n_{y}+{\tau}_{xz}n_{z}\\
    {\tau}_{yx}n_{x}+({\tau}_{yy}-p)n_{y}+{\tau}_{yz}n_{z}\\
    {\tau}_{zx}n_{x}+{\tau}_{zy}n_{y}+({\tau}_{zz}n_{z}-p)\\
  \end{bmatrix} \\& = \begin{bmatrix}
    {\tau}_{xx}n_{x}+{\tau}_{xy}n_{y}+{\tau}_{xz}n_{z}\\
    {\tau}_{yx}n_{x}+{\tau}_{yy}n_{y}+{\tau}_{yz}n_{z}\\
    {\tau}_{zx}n_{x}+{\tau}_{zy}n_{y}+{\tau}_{zz}n_{z}\\
  \end{bmatrix}
  -\begin{bmatrix}
    n_{x}\\ n_{y}\\n_{z}\\
  \end{bmatrix}p
  \end{align}
  
-
  
.. math:: 
  \begin{array}{c}
  {\tau}_{nx}={\tau}_{xx}n_{x}+{\tau}_{xy}n_{y}+{\tau}_{xz}n_{z}\\
  {\tau}_{ny}={\tau}_{yx}n_{x}+{\tau}_{yy}n_{y}+{\tau}_{yz}n_{z}\\
  {\tau}_{nz}={\tau}_{zx}n_{x}+{\tau}_{zy}n_{y}+{\tau}_{zz}n_{z}\\
  \end{array} 
  
-
  
.. math:: 
  [\boldsymbol{\sigma}\cdot \mathbf{n}]
  =\begin{bmatrix}
    {\tau}_{nx}-n_{x}p\\
    {\tau}_{ny}-n_{y}p\\
    {\tau}_{nz}-n_{z}p\\
  \end{bmatrix} 
  
-
  
.. math::  
  [\mathbf{q}\cdot \mathbf{n}]
  =q_{n}={q}_{x}n_{x}+{q}_{y}n_{y}+{q}_{z}n_{z}
  
-
  
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
  
-
  
.. math:: 
  \begin{array}{c}
  [\mathbf{q}\cdot \mathbf{n}]
  =q_{n}={q}_{x}n_{x}+{q}_{y}n_{y}+{q}_{z}n_{z}\\
  {\tau}_{vx}={\tau}_{xx}u+{\tau}_{xy}v+{\tau}_{xz}w\\
  {\tau}_{vy}={\tau}_{yx}u+{\tau}_{yy}v+{\tau}_{yz}w\\
  {\tau}_{vz}={\tau}_{zx}u+{\tau}_{zy}v+{\tau}_{zz}w\\
  \end{array}  
  
-
  
.. math::  
  [\boldsymbol{\tau}\cdot \mathbf{v}]\cdot\mathbf{n}
  =\begin{bmatrix}
    {\tau}_{xx}u+{\tau}_{xy}v+{\tau}_{xz}w\\
    {\tau}_{yx}u+{\tau}_{yy}v+{\tau}_{yz}w\\
    {\tau}_{zx}u+{\tau}_{zy}v+{\tau}_{zz}w\\
  \end{bmatrix}\cdot 
  \begin{bmatrix}
    {n}_{x}\\{n}_{y}\\{n}_{z}
  \end{bmatrix}
  =\begin{bmatrix}
    \quad {n}_{x}({\tau}_{xx}u+{\tau}_{xy}v+{\tau}_{xz}w)\\
    +{n}_{y}({\tau}_{yx}u+{\tau}_{yy}v+{\tau}_{yz}w)\\
    +{n}_{z}({\tau}_{zx}u+{\tau}_{zy}v+{\tau}_{zz}w)\\
  \end{bmatrix}  
  
-
  
.. math::   
  [\boldsymbol{\sigma}\cdot \mathbf{v}]\cdot\mathbf{n}
  =\begin{bmatrix}
    \quad {n}_{x}(({\tau}_{xx}-p)u+{\tau}_{xy}v+{\tau}_{xz}w)\\
    +{n}_{y}({\tau}_{yx}u+({\tau}_{yy}-p)v+{\tau}_{yz}w)\\
    +{n}_{z}({\tau}_{zx}u+{\tau}_{zy}v+({\tau}_{zz}-p)w)\\
  \end{bmatrix}
  =\begin{bmatrix}
    \quad {n}_{x}({\tau}_{xx}u+{\tau}_{xy}v+{\tau}_{xz}w)- {n}_{x}up\\
    +{n}_{y}({\tau}_{yx}u+{\tau}_{yy}v+{\tau}_{yz}w)- {n}_{y}vp\\
    +{n}_{z}({\tau}_{zx}u+{\tau}_{zy}v+{\tau}_{zz}w)- {n}_{z}wp\\
  \end{bmatrix} 
  
-
  
.. math:: 
  [\boldsymbol{\sigma}\cdot \mathbf{v}]\cdot\mathbf{n}
  ={n}_{x}({\tau}_{vx})+{n}_{y}({\tau}_{vy})+{n}_{z}({\tau}_{vz})- v_{n}p 
  
-
  
.. math::
  \begin{align}
  \rho E (v_{n}-{v}_{gn})+ v_{n}p & = \rho E (v_{n}-{v}_{gn})+ (v_{n}-{v}_{gn})p+{v}_{gn}p\\
  &= (\rho E +p)(v_{n}-{v}_{gn})+{v}_{gn}p\\
  &= (\rho H)(v_{n}-{v}_{gn})+{v}_{gn}p\\
  \end{align} 

-

.. math::  
  (\rho\mathbf{v})\mathbf{c}\cdot\mathbf{n}+p\mathbf{n}
  =\begin{bmatrix}
  \rho u\widetilde{v}_{n}+n_{x}p\\
  \rho v\widetilde{v}_{n}+n_{y}p\\
  \rho w\widetilde{v}_{n}+n_{z}p\\
  \end{bmatrix} 
  
-

.. math:: 
  \begin{align}
  v_{n}\ \ &=n_{x}u+n_{y}v+n_{z}w\\
  v_{gn}\ &=n_{x}u_{g}+n_{y}v_{g}+n_{z}w_{g}\\
  \widetilde{v}_{n}\ \ &=v_{n}-v_{gn}
  \end{align}