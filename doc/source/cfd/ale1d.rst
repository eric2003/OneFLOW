ALE 1D
==================================


Substitution for Integrals
---------------------------------------------------------------------------
- `Substitution for Definite Integrals <https://web.ma.utexas.edu/users/m408n/m408c/CurrentWeb/LM5-5-7.php/>`_
- `Substitution for Definite Integrals-math.libretexts.org <https://math.libretexts.org/Courses/Mount_Royal_University/MATH_1200%3A_Calculus_for_Scientists_I/4%3A_Integral_Calculus/4.7%3A_Definite_integrals_by_substitution.>`_
- `Lesson Explainer: Integration by Substitution: Definite Integrals <https://www.nagwa.com/en/explainers/816185750769/>`_
- `Substitution Rule for Indefinite Integrals <https://www.sfu.ca/math-coursenotes/Math%20158%20Course%20Notes/sec_SubRule.html>`_

In the ALE description of motion, neither the material configuration :math:`R_{\boldsymbol\xi}` nor the spatial
configuration :math:`R_{\mathbf{x}}` is taken as the reference. Thus, a third domain is needed: the referential
configuration :math:`R_{\boldsymbol\eta}` where reference coordinates :math:`{\boldsymbol\eta}` are introduced to identify the grid points.
The following figure shows these domains and the one-to-one transformations relating the configurations.
The referential domain :math:`R_{\boldsymbol\eta}` is mapped into the material and spatial domains by :math:`{\boldsymbol\Psi}`  and :math:`{\boldsymbol\Phi}`
respectively. The particle motion :math:`{\boldsymbol{\varphi}}` may then be expressed as :math:`{\boldsymbol{\varphi}}={\boldsymbol\Phi}\circ{\boldsymbol\Psi}^{-1}`, clearly showing
that, of course, the three mappings :math:`{\boldsymbol\Psi}`, :math:`{\boldsymbol\Phi}`, and :math:`{\boldsymbol{\varphi}}` are not independent.

.. figure:: images/ale2.png
   :width: 600
   :align: center
   
   The motion of the ALE computational mesh is independent of the material motion.


 
Some Example  
--------------------
  
Example 1:

.. math::
  \begin{align}
  x & = \varphi(\xi,t) = \xi +t \\
  x & = \phi(\eta,t) = \eta +2t \\
  x & =\varphi(\xi,t)=\phi(\eta,t)
  \Rightarrow \xi +t =\eta +2t 
  \Rightarrow \xi =\eta +t\\
  \xi &=\psi(\eta,t)=\eta +t\\
  \eta &=\psi^{-1}(\xi,t)=\xi-t\\
  \end{align}
  
the fluid particle velocity is

.. math::
  v_{\text{fluid}}=\cfrac{\text{d} x}{\text{d} t}\Bigg|_{\xi}=\cfrac{\partial \varphi(\xi,t)}{\partial t}=1
  
the mesh particle velocity is  
  
.. math::
  v_{\text{mesh}}=\cfrac{\text{d} x}{\text{d} t}\Bigg|_{\eta}=\cfrac{\partial \phi(\eta ,t)}{\partial t}=2\\
  
  
Let 

.. math::
  \begin{align}
  f(x,t)&=x-t\\
  f(x,t)&=x-t=f(\varphi(\xi,t),t)=(\xi +t)-t=\xi =f^{*}(\xi,t)\\
  f(x,t)&=f^{*}(\xi,t)=f^{*}(\psi(\eta,t),t)=f^{**}(\eta,t)=\xi=\eta +t\\
  \end{align}
  
then there is

.. math::
  \begin{align}
  f^{*}(\xi,t)&=\xi\\
  \cfrac{\text{d}\bar{f} }{\text{d} t}&\equiv \cfrac{\partial f^{*}(\xi,t)}{\partial t}=0\\
  \end{align}
  
-  
  
.. math::  
  \begin{align}
  \cfrac{\text{d}f(x,t)}{\text{d} t}\Bigg|_{\xi}\equiv \cfrac{\partial f^{*}(\xi,t)}{\partial t}
  &=\cfrac{\partial f(x,t)}{\partial t}+\cfrac{\partial f(x,t)}{\partial x}\cfrac{\partial x(\xi,t)}{\partial t}\\
  &=\cfrac{\partial f(x,t)}{\partial t}+\cfrac{\partial f(x,t)}{\partial x}\cfrac{\partial \varphi(\xi,t)}{\partial t}\\
  &=\cfrac{\partial (x-t)}{\partial t}+\cfrac{\partial (x-t)}{\partial x}\cfrac{\partial (\xi +t) }{\partial t}\\
  &=-1+1=0\\  
  \end{align}

-  
  
.. math::  
  \begin{align}
  \cfrac{\text{d}f^{**}(\eta,t)}{\text{d} t}\Bigg|_{\xi}\equiv \cfrac{\partial f^{*}(\xi,t)}{\partial t}
  &=\cfrac{\partial f^{**}(\eta,t)}{\partial t}+\cfrac{\partial f^{**}(\eta,t)}{\partial \eta}\cfrac{\partial \eta(\xi,t)}{\partial t}\\
  &=\cfrac{\partial f^{**}(\eta,t)}{\partial t}+\cfrac{\partial f^{**}(\eta,t)}{\partial \eta}\cfrac{\partial \psi^{-1}(\xi,t)}{\partial t}\\
  &=\cfrac{\partial (\eta +t)}{\partial t}+\cfrac{\partial (\eta +t)}{\partial \eta}\cfrac{\partial (\xi-t) }{\partial t}\\
  &=1-1=0\\  
  \end{align}
  
Example 2:  

Let 

.. math::
  \begin{align}
  f(x,t)&=x\\
  f(x,t)&=x=f(\varphi(\xi,t),t)=\xi +t =f^{*}(\xi,t)\\
  f(x,t)&=f^{*}(\xi,t)=f^{*}(\psi(\eta,t),t)=f^{**}(\eta,t)=\xi +t=\eta +2t\\
  \end{align}
  
then there is

.. math::
  \begin{align}
  f^{*}(\xi,t)&=\xi +t\\
  \cfrac{\text{d}\bar{f} }{\text{d} t}&\equiv \cfrac{\partial f^{*}(\xi,t)}{\partial t}=1\\
  \end{align}
  
-  
  
.. math::  
  \begin{align}
  \cfrac{\text{d}f(x,t)}{\text{d} t}\Bigg|_{\xi}\equiv \cfrac{\partial f^{*}(\xi,t)}{\partial t}
  &=\cfrac{\partial f(x,t)}{\partial t}+\cfrac{\partial f(x,t)}{\partial x}\cfrac{\partial x(\xi,t)}{\partial t}\\
  &=\cfrac{\partial f(x,t)}{\partial t}+\cfrac{\partial f(x,t)}{\partial x}\cfrac{\partial \varphi(\xi,t)}{\partial t}\\
  &=\cfrac{\partial (x)}{\partial t}+\cfrac{\partial (x)}{\partial x}\cfrac{\partial (\xi +t) }{\partial t}\\
  &=0+1=1\\  
  \end{align}

-  
  
.. math::  
  \begin{align}
  \cfrac{\text{d}f^{**}(\eta,t)}{\text{d} t}\Bigg|_{\xi}\equiv \cfrac{\partial f^{*}(\xi,t)}{\partial t}
  &=\cfrac{\partial f^{**}(\eta,t)}{\partial t}+\cfrac{\partial f^{**}(\eta,t)}{\partial \eta}\cfrac{\partial \eta(\xi,t)}{\partial t}\\
  &=\cfrac{\partial f^{**}(\eta,t)}{\partial t}+\cfrac{\partial f^{**}(\eta,t)}{\partial \eta}\cfrac{\partial \psi^{-1}(\xi,t)}{\partial t}\\
  &=\cfrac{\partial (\eta +2t)}{\partial t}+\cfrac{\partial (\eta +2t)}{\partial \eta}\cfrac{\partial (\xi-t) }{\partial t}\\
  &=2-1=1\\  
  \end{align}  