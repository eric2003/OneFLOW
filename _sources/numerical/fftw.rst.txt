FFTW
==================================

#. `Discrete Fourier Transform - Simple Step by Step <https://www.youtube.com/watch?v=mkGsMWi_j4Q/>`_
#. `The Two-Dimensional Discrete Fourier Transform <https://www.youtube.com/watch?v=Iz6C1ny-F2Q/>`_
#. `How the 2D FFT works <https://www.youtube.com/watch?v=v743U7gvLq0/>`_
#. `2D Fourier Transform Explained with Examples <https://www.youtube.com/watch?v=3gAZ0U66AEA/>`_
#. `2 Dimensional Discrete Fourier Transform <https://www.youtube.com/watch?v=7AddavtPWqM/>`_
#. `具体学习并实现快速傅里叶变换（FFT） <https://www.bilibili.com/video/BV1Y7411W73U/>`_
#. `fftpack <https://www.netlib.org/fftpack/>`_
#. `Netlib Repository <https://www.netlib.org/>`_
#. `Python Programming and Numerical Methods - A Guide for Engineers and Scientists <https://pythonnumericalmethods.berkeley.edu/notebooks/Index.html>`_
#. `fftw3_doc <https://www.fftw.org/fftw3_doc/index.html>`_
#. `The basic usage of FFTW <http://wuhongyi.cn/c_cppNote/ExternalLibrary/fftw3.html>`_
 

code

::

  /* determine precision and name-mangling scheme */
  #define CONCAT(prefix, name) prefix ## name
  #if defined(FFTW_SINGLE)
    typedef float R;
  # define X(name) CONCAT(fftwf_, name)
  #elif defined(FFTW_LDOUBLE)
    typedef long double R;
  # define X(name) CONCAT(fftwl_, name)
  # define TRIGREAL_IS_LONG_DOUBLE
  #elif defined(FFTW_QUAD)
    typedef __float128 R;
  # define X(name) CONCAT(fftwq_, name)
  # define TRIGREAL_IS_QUAD
  #else
    typedef double R;
  # define X(name) CONCAT(fftw_, name)
  #endif

  struct solvtab_s { void (*reg)(planner *); const char *reg_nam; };
  typedef struct solvtab_s solvtab[];
  void X(solvtab_exec)(const solvtab tbl, planner *p);