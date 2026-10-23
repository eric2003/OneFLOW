Julia
==================================

Julia 
---------------------------------
#. `Julia 教程 <https://www.runoob.com/julia/julia-tutorial.html>`_
#. `Julia中文文档 <https://docs.juliacn.com/>`_
#. `Julia Documentation <https://docs.julialang.org/>`_
#. `Package to call Python functions from the Julia language <https://github.com/JuliaPy/PyCall.jl>`_
#. `CFD_Julia系列链接整理 <https://zhuanlan.zhihu.com/p/523584688/>`_
#. `Calling C and Fortran Code <https://julia-cn.readthedocs.io/pt-br/latest/manual/calling-c-and-fortran-code.html>`_


import Pkg
::

    import Pkg

ERROR: LoadError: ArgumentError: Package CPUTime not found in current path.
- Run `import Pkg; Pkg.add("CPUTime")` to install the CPUTime package.

::

    import Pkg; Pkg.add("CPUTime")

ERROR: LoadError: ArgumentError: Package Plots not found in current path.
- Run `import Pkg; Pkg.add("Plots")` to install the Plots package.

::

    import Pkg; Pkg.add("Plots")


ERROR: LoadError: ArgumentError: Package PyPlot not found in current path.
- Run `import Pkg; Pkg.add("PyPlot")` to install the PyPlot package.
::

    import Pkg; Pkg.add("PyPlot")

ERROR: Error building `PyCall`:
┌ Info: Using the Python distribution in the Conda package by default.
└ To use a different Python version, set ENV["PYTHON"]="pythoncommand" and re-run Pkg.build("PyCall").
[ Info: Downloading miniconda installer ...
ERROR: LoadError: RequestError: schannel: failed to receive handshake, SSL/TLS connection failed while requesting https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe	

ENV["PYTHON"] = raw"c:\Users\eric\AppData\Local\Programs\Python\Python311\python.exe"
julia> Pkg.build("PyCall")









