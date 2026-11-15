Calculus
==================================

- `《高等数学》同济版 全程教学视频（宋浩老师） <https://www.bilibili.com/video/BV1Eb411u7Fw/>`_

一元复合函数的求导法则
--------------------------

定理 1
```````````````
如果函数 :math:`u=g(x)` 在点 :math:`x` 可导，而函数 :math:`y=f(u)` 在对应点
:math:`u=g(x)` 可导，那么复合函数 :math:`y=f(g(x))` 在点 :math:`x` 可导，且有

.. math::
  \frac{\mathrm{d} y}{\mathrm{d} x}={f}'(u) \cdot {g}'(x) \quad  or \quad
  \frac{\mathrm{d} y}{\mathrm{d} x}=\frac{\mathrm{d} y}{\mathrm{d} u}\cdot \frac{\mathrm{d} u}{\mathrm{d} x}
  
严格来说函数 :math:`y=f(u)` 实际上只是 :math:`u` 的函数，并不是 :math:`x` 的函数。
而 :math:`h(x)=(f\circ g)(x)` 才是 :math:`x` 的函数。这种记号上的随意会在复杂问题上造成困扰。

.. math::
  \begin{align}
    \cfrac{h(x+\triangle x)-h(x)}{\triangle x } &=\cfrac{f(g(x+\triangle x))-f(g(x))}{\triangle x }\\
    &=\cfrac{[f(g(x+\triangle x))-f(g(x))]}{[g(x+\triangle x)-g(x)]}\cfrac{[g(x+\triangle x)-g(x)]}{\triangle x }
  \end{align}
  
.. math::
  \begin{align}
    \lim_{\triangle x \to 0}\cfrac{h(x+\triangle x)-h(x)}{\triangle x } &=\lim_{\triangle x \to 0}\cfrac{f(g(x+\triangle x))-f(g(x))}{\triangle x }\\
    &=\lim_{\triangle x \to 0}\cfrac{[f(g(x+\triangle x))-f(g(x))]}{[g(x+\triangle x)-g(x)]}\cfrac{[g(x+\triangle x)-g(x)]}{\triangle x }\\
    &=\frac{\mathrm{d} f(u)}{\mathrm{d} u} \Bigg|_{u=g(x)} \cdot \frac{\mathrm{d} g(x)}{\mathrm{d} x} \Bigg|_{x=x} 
  \end{align}

也就是说，如果把函数 :math:`y` 定义成 :math:`y=f(u)` ，那么复合函数就不要记作 :math:`y` ，这有歧义。这里不妨记作 :math:`h` ,
一般来说，:math:`h` 和 :math:`f` 的表达式显然不是一回事，这是两个不同的函数。
一般有：

.. math::
  \begin{align}
    y_{0}=f(u_{0})= f(u(x_{0} ))=  f(u(x ))\Bigg|_{x=x_{0}} 
  \end{align}
  
举个例子：

.. math::
  \begin{align}
    y(u)=f(u)=2u\\
    u(x)=g(x)=sin(x)\\
    h(x)=f(g(x))=(f \circ g)(x)=2sin(x)
  \end{align}
  
显然，这里的函数 :math:`y` 和 函数 :math:`h` 不是一回事。
我们可以将变量换为任何记号，比如 :math:`y(u)=f(u)=2u` 或者 :math:`y(t)=f(t)=2t` 或者 :math:`y(x)=f(x)=2x`
对于函数 :math:`h` 有 :math:`h(x)=2sin(x)` ,如果此时也记做 :math:`y(x)=2sin(x)` 则会造成歧义。

.. math::
  \begin{align}
    y(u)=f(u)=2u\\
    u(x)=g(x)=sin(x)\\
   \hat{y}(x)=h(x)=f(g(x))=(f \circ g)(x)=2sin(x)
  \end{align}
  
这里的函数 :math:`y` 和 函数 :math:`\hat{y}` 应该加以区分。

.. math::
  y(sin(x))=f(sin(x))=\hat{y}(x)=2sin(x)

从这个角度来说，定理更合理的表述应该是：

如果函数 :math:`u=g(x)` 在点 :math:`x` 可导，而函数 :math:`y=f(u)` 在对应点
:math:`u=g(x)` 可导，那么复合函数 :math:`\hat{y}=f(g(x))` 在点 :math:`x` 可导，且有

.. math::
  \frac{\mathrm{d} \hat{y}(x)}{\mathrm{d} x}={f}'(u) \cdot {g}'(x) \quad  or \quad
  \frac{\mathrm{d} \hat{y}(x)}{\mathrm{d} x}=\frac{\mathrm{d} y(u)}{\mathrm{d} u}\cdot \frac{\mathrm{d} u(x)}{\mathrm{d} x}
  

多元复合函数的求导法则
--------------------------

定理 1
```````````````
如果函数 :math:`u=\phi(t)` 及 :math:`v=\psi(t)` 都在点 :math:`t` 可导，函数 :math:`z=f(u,v)` 在对应点
:math:`(u,v)` 具有连续偏导数，那么复合函数 :math:`z=f(\phi(t),\psi(t))` 在点 :math:`t` 可导，且有

.. math::
  \frac{\mathrm{d} z}{\mathrm{d} t}=\frac{\partial z}{\partial u} \frac{\mathrm{d} u}{\mathrm{d} t}+
  \frac{\partial z}{\partial v} \frac{\mathrm{d} v}{\mathrm{d} t}
  
更准确的叙述：
如果函数 :math:`u=\phi(t)` 及 :math:`v=\psi(t)` 都在点 :math:`t` 可导，函数 :math:`z(u,v)=f(u,v)` 在对应点
:math:`(u,v)` 具有连续偏导数，那么复合函数 :math:`\hat{z}(t)=f(\phi(t),\psi(t))` 在点 :math:`t` 可导，且有

.. math::
  \frac{\mathrm{d} \hat{z}(t)}{\mathrm{d} t}=\frac{\partial z(u,v)}{\partial u} \frac{\mathrm{d} u(t)}{\mathrm{d} t}+
  \frac{\partial z(u,v)}{\partial v} \frac{\mathrm{d} v(t)}{\mathrm{d} t}  
  
即：

.. math::
  \frac{\mathrm{d} \hat{z}(t)}{\mathrm{d} t}=\frac{\partial z(u,v)}{\partial u} \frac{\mathrm{d} u(t)}{\mathrm{d} t}+
  \frac{\partial z(u,v)}{\partial v} \frac{\mathrm{d} v(t)}{\mathrm{d} t}  

定理 2
```````````````

如果函数 :math:`u=\phi(x,y)` 及 :math:`v=\psi(x,y)` 都在点 :math:`(x,y)` 具有对 :math:`x` 及对 :math:`y` 的偏导数，函数 :math:`z=f(u,v)` 在对应点
:math:`(u,v)` 具有连续偏导数，那么复合函数 :math:`z=f(\phi(x,y),\psi(x,y))` 在点 :math:`(x,y)` 的两个偏导数都存在，且有

.. math::
  \begin{align}
    \frac{\partial z}{\partial x} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial x}+
    \frac{\partial z}{\partial v} \frac{\partial v}{\partial x}\\
    \frac{\partial z}{\partial y} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial y}+
    \frac{\partial z}{\partial v} \frac{\partial v}{\partial y}\\
  \end{align}
  
类似的，设 :math:`u=\phi(x,y)` 、:math:`v=\psi(x,y)` 及 :math:`w=\omega(x,y)` 都在点 :math:`(x,y)` 具有对 :math:`x` 及对 :math:`y` 的偏导数，函数 :math:`z=f(u,v,w)` 在对应点
:math:`(u,v,w)` 具有连续偏导数，那么复合函数 

.. math::
  z=f(\phi(x,y),\psi(x,y),\omega(x,y))

在点 :math:`(x,y)` 的两个偏导数都存在，且有

.. math::
  \begin{align}
    \frac{\partial z}{\partial x} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial x}+
    \frac{\partial z}{\partial v} \frac{\partial v}{\partial x}+
    \frac{\partial z}{\partial w} \frac{\partial w}{\partial x}\\
    \frac{\partial z}{\partial y} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial y}+
    \frac{\partial z}{\partial v} \frac{\partial v}{\partial y}+
    \frac{\partial z}{\partial w} \frac{\partial w}{\partial y}
  \end{align}
  
更准确的叙述：  
如果函数 :math:`u=\phi(x,y)` 及 :math:`v=\psi(x,y)` 都在点 :math:`(x,y)` 具有对 :math:`x` 及对 :math:`y` 的偏导数，函数 :math:`z(u,v)=f(u,v)` 在对应点
:math:`(u,v)` 具有连续偏导数，那么复合函数 :math:`\hat{z}(x,y)=f(\phi(x,y),\psi(x,y))` 在点 :math:`(x,y)` 的两个偏导数都存在，且有

.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial z(u,v)}{\partial u} \frac{\partial u(x,y)}{\partial x}+
    \frac{\partial z(u,v)}{\partial v} \frac{\partial v(x,y)}{\partial x}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial z(u,v)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial z(u,v)}{\partial v} \frac{\partial v(x,y)}{\partial y}\\
  \end{align}
  
类似的，设 :math:`u=\phi(x,y)` 、:math:`v=\psi(x,y)` 及 :math:`w=\omega(x,y)` 都在点 :math:`(x,y)` 具有对 :math:`x` 及对 :math:`y` 的偏导数，函数 :math:`z(u,v,w)=f(u,v,w)` 在对应点
:math:`(u,v,w)` 具有连续偏导数，那么复合函数 

.. math::
  \hat{z}(x,y)=f(\phi(x,y),\psi(x,y),\omega(x,y))

在点 :math:`(x,y)` 的两个偏导数都存在，且有

.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial x}+
    \frac{\partial z(u,v,w)}{\partial v} \frac{\partial v(x,y)}{\partial x}+
    \frac{\partial z(u,v,w)}{\partial w} \frac{\partial w(x,y)}{\partial x}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial z(u,v,w)}{\partial v} \frac{\partial v(x,y)}{\partial y}+
    \frac{\partial z(u,v,w)}{\partial w} \frac{\partial w(x,y)}{\partial y}
  \end{align}  
  
定理 3
```````````````

如果函数 :math:`u=\phi(x,y)` 在点 :math:`(x,y)` 具有对 :math:`x` 及对 :math:`y`的偏导数，函数 :math:`v=\psi(y)` 在点 :math:`y` 可导，函数 :math:`z=f(u,v)` 在对应点
:math:`(u,v)` 具有连续偏导数，那么复合函数 :math:`z=f(\phi(x,y),\psi(y))` 在点 :math:`(x,y)` 的两个偏导数都存在，且有

.. math::
  \begin{align}
    \frac{\partial z}{\partial x} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial x}\\
    \frac{\partial z}{\partial y} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial y}+
    \frac{\partial z}{\partial v} \frac{\mathrm{d} v}{\mathrm{d} y}\\
  \end{align}  
  
如果复合函数的某些中间变量本身又是复合函数的自变量
设 :math:`u=\phi(x,y)` 、:math:`v=x` 及 :math:`w=y` 都在点 :math:`(x,y)` 具有对 :math:`x` 及对 :math:`y` 的偏导数，函数 :math:`z=f(u,v,w)` 即 :math:`z=f(u,x,y)` 在对应点
:math:`(u,v,w)` 具有连续偏导数，那么复合函数 

.. math::
  z=f(\phi(x,y),x,y)

有  

.. math::
  \begin{align}
    \frac{\partial v}{\partial x}&=1, \quad \frac{\partial v}{\partial y}=0\\
    \frac{\partial w}{\partial x}&=0, \quad \frac{\partial w}{\partial y}=1
  \end{align}
  
继续，有 

.. math::
  \begin{align}
    \frac{\partial z}{\partial x} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial x}+
    \frac{\partial z}{\partial v} \cdot 1+
    \frac{\partial z}{\partial w} \cdot 0\\
    \frac{\partial z}{\partial y} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial y}+
    \frac{\partial z}{\partial v} \cdot 0+
    \frac{\partial z}{\partial w} \cdot 1
  \end{align}
  
继续，有 

.. math::
  \begin{align}
    \frac{\partial z}{\partial x} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial x}+
    \frac{\partial z}{\partial v}\\
    \frac{\partial z}{\partial y} & = \frac{\partial z}{\partial u} \frac{\partial u}{\partial y}+
    \frac{\partial z}{\partial w}
  \end{align}
  
更准确的叙述：  
如果函数 :math:`u=\phi(x,y)` 在点 :math:`(x,y)` 具有对 :math:`x` 及对 :math:`y` 的偏导数，函数 :math:`v=\psi(y)` 在点 :math:`y` 可导，函数 :math:`z(u,v)=f(u,v)` 在对应点
:math:`(u,v)` 具有连续偏导数，那么复合函数 :math:`\hat{z}(x,y)=f(\phi(x,y),\psi(y))=\hat{f}(x,y)` 在点 :math:`(x,y)` 的两个偏导数都存在，且有

.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial z(u,v)}{\partial u} \frac{\partial u(x,y)}{\partial x}+
    \frac{\partial z(u,v)}{\partial v} \cancelto{0}{\frac{\partial v(x,y)}{\partial x}}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial z(u,v)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial z(u,v)}{\partial v} \cancelto{\frac{\mathrm{d} v(y)}{\mathrm{d} y}}{\frac{\partial v(x,y)}{\partial y}}
  \end{align}    

.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial z(u,v)}{\partial u} \frac{\partial u(x,y)}{\partial x}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial z(u,v)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial z(u,v)}{\partial v} {\frac{\mathrm{d} v(y)}{\mathrm{d} y}}
  \end{align}    

  
如果复合函数的某些中间变量本身又是复合函数的自变量
设 :math:`u=\phi(x,y)` 、:math:`v=x` 及 :math:`w=y` 都在点 :math:`(x,y)` 具有对 :math:`x` 及对 :math:`y` 的偏导数，函数 :math:`z(u,v,w)=f(u,v,w)` 即 :math:`z=f(u,x,y)` 在对应点
:math:`(u,v,w)` 具有连续偏导数，那么复合函数 

.. math::
  \hat{z}(x,y)=f(\phi(x,y),x,y)=\hat{f}(x,y)
  
有：  
  
.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial x}+
    \frac{\partial z(u,v,w)}{\partial v} \frac{\partial v(x,y)}{\partial x}+
    \frac{\partial z(u,v,w)}{\partial w} \frac{\partial w(x,y)}{\partial x}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial z(u,v,w)}{\partial v} \frac{\partial v(x,y)}{\partial y}+
    \frac{\partial z(u,v,w)}{\partial w} \frac{\partial w(x,y)}{\partial y}
  \end{align} 

继续，有：  
  
.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial x}+
    \frac{\partial z(u,v,w)}{\partial v} \cancelto{1}{\frac{\partial v(x,y)}{\partial x}}+
    \frac{\partial z(u,v,w)}{\partial w} \cancelto{0}{\frac{\partial w(x,y)}{\partial x}}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial z(u,v,w)}{\partial v} \cancelto{0}{\frac{\partial v(x,y)}{\partial y}}+
    \frac{\partial z(u,v,w)}{\partial w} \cancelto{1}{\frac{\partial w(x,y)}{\partial y}}
  \end{align}

即：

.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial x}+
    \frac{\partial z(u,v,w)}{\partial v}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial z(u,v,w)}{\partial w}
  \end{align}
  
即：

.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial x}+
    \frac{\partial z(u,v,w)}{\partial x}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial z(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial z(u,v,w)}{\partial y}
  \end{align}  
  
也可以写成：

.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial f(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial x}+
    \frac{\partial f(u,v,w)}{\partial x}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial f(u,v,w)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial f(u,v,w)}{\partial y}
  \end{align}  

继续，也可以写成：

.. math::
  \begin{align}
    \frac{\partial \hat{z}(x,y)}{\partial x} & = \frac{\partial f(u,x,y)}{\partial u} \frac{\partial u(x,y)}{\partial x}+
    \frac{\partial f(u,x,y)}{\partial x}\\
    \frac{\partial \hat{z}(x,y)}{\partial y} & = \frac{\partial f(u,x,y)}{\partial u} \frac{\partial u(x,y)}{\partial y}+
    \frac{\partial f(u,x,y)}{\partial y}
  \end{align} 
  
需要注意的是，虽然数值上有：

.. math::
  \begin{align}
    \hat{z}(x,y)&= \hat{f}(x,y)=z(u(x,y),x,y)=f(u(x,y),x,y)\\
  \end{align}
  
但是，函数形式上  

.. math::
  \begin{align}
    \hat{z}&\ne z\\
    \hat{f}&\ne f\\
    \hat{z}&\equiv \hat{f}\\
    {z}&\equiv {f}\\
  \end{align}
