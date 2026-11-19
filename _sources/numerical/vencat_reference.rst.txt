Vencat analytic reference profile
==================================

Purpose
-------

The Vencat profile is retained as a small analytic-data generator for
numerical experiments. The previous FunctionTest runtime task was not
referenced by any checked-in simu.txt case, and its implementation mixed
formula definitions, output formatting, and production task registration.
The formulas and profile writer now live in
tools/verification/vencat_reference.py, outside the production C++ source
set.

Formulas
--------

The three historical functions are:

.. math::

   V(x)=\frac{1+2x}{1+x+2x^2}

.. math::

   V_c(x)=\frac{1+cx}{1+x+cx^2}

.. math::

   V_{\mathrm{Eric}}(x)=\frac{1+3x}{1+x+2x^2}

The reference-data writer uses V_c(x) with c=2,3,...,11.
For each value of c, it writes N=101 points over x in [0,2] and stores
x, V_c(x), and xV_c(x) in a separate zone.

Generation
----------

Choose an explicit output path:

.. code-block:: console

   python -m tools.verification.vencat_reference build/verification/vencat.dat

The output path may be relative to the current working directory or absolute.
Parent directories are created when needed. The generator does not write into
the repository source tree unless the caller explicitly chooses such a path.

Validation
----------

The focused tests are:

.. code-block:: console

   python -m unittest tools.verification.test_vencat_reference

They check representative formula values, the number of zones and points,
and the first/last coordinates and profile values. These checks validate the
analytic generator and output format only. They do not establish CFD solver
accuracy, discretization convergence, or agreement with experimental data.
