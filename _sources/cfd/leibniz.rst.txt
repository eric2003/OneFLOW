Leibniz integral
==================================

Integration by substitution
--------------------------------------
.. math:: 
  \int_{a}^{b} f(x)dx=\int_{\alpha }^{\beta } f(\varphi(t)){\varphi}' (t)dt

Leibniz integral rule
--------------------------------------

#. `Leibniz integral rule <https://en.wikipedia.org/wiki/Leibniz_integral_rule>`_
#. `How to take the derivative of the integral? <https://zhuanlan.zhihu.com/p/547766124>`_


General form: differentiation under the integral sign

.. math:: 

  \int_{a(x)}^{b(x)} f\left ( x,t \right ) dt
  
Theorem — Let f(x,t) be a function such that both :math:`f\left ( x,t \right )` are continuous it t and x in some 
region of the xt-plane, including :math:`a(x)\le t\le b(x),x_{0} \le x\le x_{1}`. Also suppose that the functions :math:`a(x)` and
functions :math:`b(x)` are both continuous and both have continuous derivatives for :math:`x_{0} \le x\le x_{1}`.Then for :math:`x_{0} \le x\le x_{1}`,

.. math::
  {\displaystyle {\frac {d}{dx}}\left(\int _{a(x)}^{b(x)}f(x,t)\,dt\right)=f{\big (}x,b(x){\big )}\cdot {\frac {d}{dx}}b(x)-f{\big (}x,a(x){\big )}\cdot {\frac {d}{dx}}a(x)+\int _{a(x)}^{b(x)}{\frac {\partial }{\partial x}}f(x,t)\,dt.}  
  
The right hand side may also be written using Lagrange's notation as: 

.. math::
  {\textstyle f(x,b(x))\,b^{\prime }(x)-f(x,a(x))\,a^{\prime }(x)+\displaystyle \int _{a(x)}^{b(x)}f_{x}(x,t)\,dt.}
  
another form

.. math::
  \Phi(x)=\int_{\alpha(x) }^{\beta(x)} f(x, y) dy
  
.. math::
  \Phi(t)=\int_{\alpha(t) }^{\beta(t)} f(t, x) dx
  
.. math::
  \Phi(t)=\int_{\alpha(t) }^{\beta(t)} F(x, t) dx

.. math::
  \begin{align}
    \frac{\mathrm{d} \Phi(t)}{\mathrm{d} t}&=\frac{\Phi(t+\Delta t)-\Phi(t)}{\Delta t}\\
  &=\int_{\alpha(t) }^{\beta(t)}\frac{ F(x, t+\Delta t)-F(x, t)}{\Delta t} dx\\
   &+\frac{1}{\Delta t}\int_{\beta(t) }^{\beta(t+\Delta t)}{ F(x, t+\Delta t)} dx \\
   &-\frac{1}{\Delta t}\int_{\alpha(t) }^{\alpha(t+\Delta t)}{ F(x, t+\Delta t)} dx 
  \end{align}

.. math::
  \begin{align}
    \frac{\mathrm{d} \Phi(t)}{\mathrm{d} t}&=\frac{\mathrm{d} }{\mathrm{d} t}\int_{\alpha(t) }^{\beta(t)}F(x, t) dx\\\\
  &=\int_{\alpha(t) }^{\beta(t)}\frac{\partial F(x, t)}{\partial t} dx
   +[ F(\beta(t), t)\dot{\beta}(t) - F(\alpha(t), t)\dot{\alpha}(t) ]\\
  \end{align}  
  
This can also be written as  

.. math::
  \cfrac{\mathrm{d}}{\mathrm{d}t}{\int_{g(t)}^{h(t)}F(x,t)\mathrm{d}x }={\int_{g(t)}^{h(t)}\frac{\partial F(x,t)}{\partial t}\mathrm{d}x }+\left \{ F[h(t),t]\dot{h}(t)-F[g(t),t]\dot{g}(t) \right \}
  
One proof runs as follows, modulo precisely stated hypotheses and some analytic 
details. Set 

.. math::
  \Phi(u, v, t)=\int_{u}^{v} F(x, t) d x
  
:math:`u = g(t)`, and :math:`v = h(t)`. By the chain rule

.. math::
  \frac{d}{d t} \Phi[g(t), h(t), t]=\left(\frac{\partial \Phi}{\partial u} \dot{g}+\frac{\partial \Phi}{\partial v} h\right)+\frac{\partial \Phi}{\partial t}  
  
The first two terms are bracketed because they measure all changes due to variation 
of the interval of integration [g(t), h(t)], and they are evaluated by applying the 
Fundamental Theorem to :math:`\Phi(u, v, t)=\int_{u}^{v} F(x, t) d x`. The third term measures change due to variation of 
the integrand. If enough smoothness is assumed to justify interchange of the integration and differentiation operators, then

.. math::
  \frac{\partial \Phi}{\partial t}=\frac{\partial}{\partial t} \int_{u}^{v} F(x, t) d x=\int_{u}^{v} \frac{\partial F(x, t)}{\partial t} d x .
  
Another proof

.. math::
  \begin{aligned}
    \phi_{t}(u) & =x(u, t), \\
    \phi_{t}:[a, b] & \rightarrow\left[\phi_{t}(a), \phi_{t}(b)\right]=[g(t), h(t)]=C_{t} .
  \end{aligned}
  
By the formula for change of variable in a simple integral

.. math::
  \int_{g(t)}^{h(t)} F(x) d x=\int_{\phi_{t}(a)}^{\phi_{t}(b)} F(x) d x=\int_{a}^{b} F[x(u, t)] \frac{\partial x}{\partial u} d u .
  
This transition is excellent, because it has changed the integral over a moving domain 
to one over a fixed domain. We pay for this fixed domain with a time-varying integrand. No matter, we like it; we thrive on differentiation under the integral sign:

.. math::
  \begin{aligned}
    \frac{\mathrm{d} }{\mathrm{d} t} \int_{g(t)}^{h(t)} F(x) d x 
    &=\frac{\mathrm{d} }{\mathrm{d} t} \int_{a}^{b} F[x(u, t)] \frac{\partial x(u, t)}{\partial u} d u \\
    &=\int_{a}^{b}\frac{\partial}{\partial t}\left \{  F[x(u, t)] \frac{\partial x(u, t)}{\partial u}\right \}d u \\
    &=\int_{a}^{b}\left \{\frac{\partial F[x(u, t)]}{\partial x}\frac{\partial x(u, t)}{\partial t}  \frac{\partial x(u, t)}{\partial u}+F[x(u, t)]\frac{\partial x^{2} (u, t)}{\partial u\partial t}\right \}d u \\
  \end{aligned}
  
The fixed domain has done its job, and we return to the moving domain. The instantaneous velocity is :math:`v = v(u, t) = \partial x(u, t)/ {\partial t}`, which we also consider as a function of :math:`x` 
and :math:`t` via the transformation :math:`(u, t) \leftrightarrow (x, t)`. When :math:`t` is fixed,

.. math::
  \begin{aligned}
    \frac{\partial x^{2} (u, t)}{\partial u\partial t}&=\frac{\partial v(u, t)}{\partial u}\\
    &=\cfrac{\cfrac{\partial v(u, t)}{\partial u}}{\cfrac{\partial x(u, t)}{\partial u}} {\frac{\partial x(u, t)}{\partial u}} \\
    &={\cfrac{\partial v(u, t)}{\partial x(u, t)}}{\cfrac{\partial x(u, t)}{\partial u}}
  \end{aligned}
  
hence

.. math::
  \begin{aligned}
    &\int_{a}^{b}\left \{\frac{\partial F[x(u, t)]}{\partial x}\frac{\partial x(u, t)}{\partial t}  \frac{\partial x(u, t)}{\partial u}+F[x(u, t)]\frac{\partial x^{2} (u, t)}{\partial u\partial t}\right \}d u \\
    =&\int_{a}^{b}\left \{\frac{\partial F[x(u, t)]}{\partial x}\frac{\partial x(u, t)}{\partial t}  \frac{\partial x(u, t)}{\partial u}+F[x(u, t)]{\cfrac{\partial v(u, t)}{\partial x(u, t)}}{\cfrac{\partial x(u, t)}{\partial u}}\right \}d u \\
    =&\int_{a}^{b}\left \{\frac{\partial F[x(u, t)]}{\partial x}\frac{\partial x(u, t)}{\partial t}  +F[x(u, t)]{\cfrac{\partial v(u, t)}{\partial x(u, t)}}\right \}{\cfrac{\partial x(u, t)}{\partial u}}d u \\
    =&\int_{\phi_{t}(a)}^{\phi_{t}(b)}\left \{\frac{\partial F[x(u, t)]}{\partial x}\frac{\partial x(u, t)}{\partial t}  +F[x(u, t)]{\cfrac{\partial v(u, t)}{\partial x(u, t)}}\right \}d x \\  
  \end{aligned}
  
that is  

.. math::
  \begin{aligned}
    &\int_{\phi_{t}(a)}^{\phi_{t}(b)}\left \{\frac{\partial F[x(u, t)]}{\partial x}\frac{\partial x(u, t)}{\partial t}  +F[x(u, t)]{\cfrac{\partial v(u, t)}{\partial x(u, t)}}\right \}d x \\
    =&\int_{\phi_{t}(a)}^{\phi_{t}(b)}\frac{\partial }{\partial x} [F[x(u, t)]v(u, t)]d x \\
    =&\int_{g(t)}^{h(t)}\frac{\partial }{\partial x} [F[x(u, t)]v(u, t)]d x 
  \end{aligned}
  
Two-dimensional, time-dependent case

We are also given a function :math:`F(x, y, t)`. The problem is to find

.. math::
  \frac{\mathrm{d} }{\mathrm{d} t}\iint_{D(t)}^{} F(x,y,t)dxdy
  
Certainly our first move should be separation of boundary variation from integrand variation. This is easy enough by the chain rule device in the first section 
and results in

.. math::
  \begin{aligned}
    &\frac{\mathrm{d} }{\mathrm{d} t}\iint_{D(t)}^{} F(x,y,t)dxdy{\Bigg|}_{t=t_{0} }\\
    =&\frac{\mathrm{d} }{\mathrm{d} t}\iint_{D(t)}^{} F(x,y,t_{0})dxdy{\Bigg|}_{t=t_{0} }+
    \iint_{D(t_{0})}^{} \frac{\partial F(x,y,t)}{\partial t}{\Bigg|}_{t=t_{0} }dxdy
  \end{aligned}
  
This is routine. The essence of the problem is to find

.. math::
  \frac{\mathrm{d} }{\mathrm{d} t}\iint_{D(t)}^{} F(x,y)dxdy  
  
Let :math:`\mathbf{v}=\mathbf{v}(x, y, t)` denote the velocity vector at a boundary point :math:`(x, y)` of :math:`Dt` and let :math:`n` 
denote the outward unit normal. In the difference

.. math::
  \iint_{D(t+\Delta t )}^{} F(x,y)dxdy-\iint_{D(t)}^{} F(x,y)dxdy
  
everything in the overlap of :math:`D(t)` and :math:`D(t+\Delta t)` cancels; only the thin boundary strip makes 
a contribution. From the detail, this contribution is  

.. math::
  F(x,y)(\mathbf{v}\Delta t)\cdot(\mathbf{n}ds)
  
up to higher order differentials, where :math:`ds`  is the element of arc length. 
Hence  

.. math::
  \begin{aligned}
    &\lim_{\Delta t \to 0} \cfrac{1}{\Delta t}\left \{  \iint_{D(t+\Delta t )}^{} F(x,y)dxdy-\iint_{D(t)}^{} F(x,y)dxdy\right \}\\
    =&\lim_{\Delta t \to 0}\cfrac{1}{\Delta t}\int_{\partial D(t)}^{} F(x,y)(\mathbf{v}\Delta t)\cdot(\mathbf{n}ds)
  \end{aligned}

.. math::
  \cfrac{\mathrm{d} }{\mathrm{d} t}\iint_{D(t)}^{} F(x,y)dxdy  =
  \int_{\partial D(t)}^{} F(x,y){\mathbf{v}}\cdot{\mathbf{n}}ds
  
where :math:`\partial` a denotes boundary. Before taking limits, we compute :math:`{\mathbf{v}}\cdot{\mathbf{n}}ds`. We rotate the
unit tangent :math:`(dy/ds,-dx/ds)` backwards through a right angle to obtain :math:`\mathbf{n}=(dy/ds,-dx/ds)`, hence  

.. math::
  \mathbf{n}=(dy/ds,-dx/ds)=(u,v)\cdot(dy,-dx)=udy-vdx
  
Therefore

.. math::
  \begin{aligned}
    \cfrac{\mathrm{d} }{\mathrm{d} t}\iint_{D(t)}^{} F(x,y)dxdy  =&
    \int_{\partial D(t)}^{} F(x,y){\mathbf{v}}\cdot{\mathbf{n}}ds \\
    =&\int_{\partial D(t)}^{} F(x,y)(u,v)\cdot(dy,-dx) \\
    =&\int_{\partial D(t)}^{} F(x,y)(udy-vdx)\\
    =&\int_{\partial D(t)}^{} (-F(x,y)vdx+F(x,y)udy)
  \end{aligned}
  
We can transform the boundary integral into an integral over D, by Green's Theorem. 

.. math::
   \begin{aligned}
     P(x,y)&\equiv -F(x,y)v(x,y)\\
     Q(x,y)&\equiv F(x,y)u(x,y)
   \end{aligned}
   
.. math::
   \begin{aligned}
     (\cfrac{\partial Q(x,y)}{\partial x}-\cfrac{\partial P(x,y)}{\partial y})=
     &(\cfrac{\partial [F(x,y)u(x,y)]}{\partial x}-\cfrac{\partial [-F(x,y)v(x,y)]}{\partial y})\\
     =&(\cfrac{\partial [F(x,y)u(x,y)]}{\partial x}+\cfrac{\partial [F(x,y)v(x,y)]}{\partial y})\\
     =&(\cfrac{\partial [F(x,y)u(x,y)]}{\partial x}+\cfrac{\partial [F(x,y)v(x,y)]}{\partial y})\\
   \end{aligned}   
   
.. math::
  \begin{aligned}
    \cfrac{\partial [F(x,y)u(x,y)]}{\partial x}-\cfrac{\partial [-F(x,y)v(x,y)]}{\partial y}=&\\
    \cfrac{\partial [F(x,y)u(x,y)]}{\partial x}+\cfrac{\partial [F(x,y)v(x,y)]}{\partial y}=\mathrm{div} \left \{ F(x,y)\mathbf{v(x,y)}  \right \} \\  
    \cfrac{\partial (F u)}{\partial x}+\cfrac{\partial (F v)}{\partial y}=\mathrm{div} ( F\mathbf{v} )
  \end{aligned}
  
.. math::
  \cfrac{\partial (F u)}{\partial x}+\cfrac{\partial (F v)}{\partial y}=\mathrm{div} ( F\mathbf{v} )= \nabla \cdot( F\mathbf{v} )
  
Here

.. math::
  \begin{aligned}
    \mathrm{div} ( F\mathbf{v} )\equiv\nabla \cdot( F\mathbf{v} )=&\cfrac{\partial (F u)}{\partial x}+\cfrac{\partial (F v)}{\partial y}\\
    =&F \nabla \cdot \mathbf{v}+\mathbf{v} \cdot \nabla F  \\
    =&\nabla F \cdot \mathbf{v} + F \nabla \cdot \mathbf{v}
  \end{aligned}
  
.. math::
  \nabla \cdot(\phi  \mathbf{v})=\phi \nabla \cdot \mathbf{v}+\mathbf{v} \cdot \nabla \phi  

.. math::
  \nabla \cdot(F  \mathbf{v})=F \nabla \cdot \mathbf{v}+\mathbf{v} \cdot \nabla F  

.. math::
  \begin{aligned}
    &\frac{\mathrm{d} }{\mathrm{d} t}\iint_{D(t)}^{} F(x,y,t)dxdy\\
  \end{aligned} 
  
  
.. math::
  \begin{aligned}
    \frac{\mathrm{d} }{\mathrm{d} t}\iint_{D(t)}^{} F(x,y,t)dxdy
    =\int_{\partial D(t)}^{} F(x,y){\mathbf{v}}\cdot{\mathbf{n}}ds+
    \iint_{D(t)}^{} \frac{\partial F(x,y,t)}{\partial t}dxdy
  \end{aligned}
  
.. math::
  \begin{aligned}
    \frac{\mathrm{d} }{\mathrm{d} t}\iint_{D(t)}^{} F(x,y,t)dxdy
    =&\iint_{D(t)}^{} \nabla \cdot(F  \mathbf{v})dxdy+
    \iint_{D(t)}^{} \frac{\partial F(x,y,t)}{\partial t}dxdy\\
    =&\iint_{D(t)}^{} \left \{\nabla \cdot(F  \mathbf{v})+\cfrac{\partial F(x,y,t)}{\partial t}\right \}dxdy
  \end{aligned}   

A space formula.

.. math::
  \begin{aligned}
    \frac{d}{d t} \iiint_{D_{t}} F(\mathbf{x}, t) d x d y d z & =\iint_{\partial D_{t}} F \mathbf{v} \cdot \mathbf{n} dS+\iiint_{D_{t}} \frac{\partial F}{\partial t} d x d y d z \\
  \end{aligned}

.. math::
  \begin{aligned}
    \frac{d}{d t} \iiint_{D_{t}} F(\mathbf{x}, t) d x d y d z & =\iint_{\partial D_{t}} F \mathbf{v} \cdot d\mathbf{S} +\iiint_{D_{t}} \frac{\partial F}{\partial t} d x d y d z \\
  \end{aligned}

.. math::
  \begin{aligned}
    \frac{d}{d t} \iiint_{D_{t}} F(\mathbf{x}, t) d x d y d z & =
    \iiint_{D_{t}}\left[\operatorname{div}(F \mathbf{v})+\frac{\partial F}{\partial t}\right] d x d y d z\\
    &=\iiint_{D_{t}}\left[\nabla \cdot(F  \mathbf{v})+\frac{\partial F}{\partial t}\right] d x d y d z\\
  \end{aligned}
  
Green's theorem
--------------------------------------
#. `Green's theorem <https://en.wikipedia.org/wiki/Green%27s_theorem>`_

Let C be a positively oriented, piecewise smooth, simple closed curve in a plane, and let D be the region bounded by C. If L and M are functions of (x, y) defined on an open region containing D and have continuous partial derivatives there, then

.. math::
  \oint_{C}^{} (Ldx+Mdy)=\iint_{D}^{} (\cfrac{\partial M}{\partial x}-\cfrac{\partial L}{\partial y})dxdy
  
where the path of integration along C is anticlockwise.

This formula can also be written as

.. math::
  \oint_{C}^{} (P(x,y)dx+Q(x,y)dy)=\iint_{D}^{} (\cfrac{\partial Q}{\partial x}-\cfrac{\partial P}{\partial y})dxdy\\

  
A Leibniz integral rule for a two dimensional surface moving in three dimensional space is[3][4]

.. math::
  {\displaystyle {\frac {d}{dt}}\iint _{\Sigma (t)}\mathbf {F} (\mathbf {r} ,t)\cdot d\mathbf {A} =\iint _{\Sigma (t)}\left(\mathbf {F} _{t}(\mathbf {r} ,t)+\left[\nabla \cdot \mathbf {F} (\mathbf {r} ,t)\right]\mathbf {v} \right)\cdot d\mathbf {A} -\oint _{\partial \Sigma (t)}\left[\mathbf {v} \times \mathbf {F} (\mathbf {r} ,t)\right]\cdot d\mathbf {s} ,}
