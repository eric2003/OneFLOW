Environment Modules 
==================================

Environment Modules is a tool for managing applications and environment variables, allowing users to dynamically load and unload software packages as needed, and configure the necessary environment variables for applications as required. Using Environment Modules, users can run multiple versions of applications on the same system without version conflicts or dependency issues. In addition, Environment Modules can help administrators simplify the installation and management of software packages, thereby improving the maintainability and reliability of the system.

#. `Environment Modules Official Website <http://modules.sourceforge.net/>`_
#. `Environment Modules Wikipedia page <https://en.wikipedia.org/wiki/Environment_Modules>`_
#. `Environment Modules usage tutorial:  <https://modules.readthedocs.io/>`_
#. `Environment Modules GitHub page <https://github.com/cea-hpc/modules/>`_
#. `A Compilation of Links Related to Environment Modules, an Environment Variable Management Tool <https://zhuanlan.zhihu.com/p/559136017/>`_

::

  module avail

cmake-3.28.0 modulefile
::

  #%Module1.0
  ##
  ##
  module-whatis "cmake-3.28.0 modulefile"
  set topdir "/usr/local/cmake-3.28.0-linux-x86_64"
  prepend-path PATH "${topdir}/bin"  
  
modulefiles location
::

  /home/eric/software/Modules/modulefiles

load module
::

  module load eric/cmake-3.28.0
 



