.. Plonk documentation master file, created by
   sphinx-quickstart on Mon Mar 18 15:51:00 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====
Plonk
=====

Smoothed particle hydrodynamics analysis and visualization with Python.

* Docs: https://plonk.readthedocs.io/
* Repo: https://www.github.com/dmentipl/plonk

.. warning:: Plonk is under development. It may contain bugs. The authors do not take responsibility for use of Plonk.

--------
Overview
--------

Plonk is a Python tool for analysis and visualization of smoothed particle hydrodynamics data.

Plonk aims to provide an alternative to the Phantom analysis tools and Splash. Both are compiled codes, and as such they are less flexible and scriptable than interpreted Python code.

.. note:: Plonk contains Fortran interpolation subroutines from Splash as part of the visualization module.

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   installation
   getting-started
   examples

.. role:: math(raw)
   :format: html latex
..

----------
References
----------

* Plonk:
   * https://github.com/dmentipl/plonk
   * Mentiplay, D., 2019, ASCL, ascl:1907.009 [`ADS <https://ui.adsabs.harvard.edu/abs/2019ascl.soft07009M/>`_, `DOI <http://ascl.net/1907.009>`_]
* Splash:
   * http://users.monash.edu.au/~dprice/splash/
   * Price D. J., 2007, PASA, 24, 159 [`ADS <https://ui.adsabs.harvard.edu/abs/2007PASA...24..159P>`_, `DOI <http://dx.doi.org/10.1071/AS07022>`_]
* Phantom:
   * https://phantomsph.bitbucket.io/
   * Price D. J., et al., 2018, PASA, 35, 31P [`ADS <https://ui.adsabs.harvard.edu/abs/2018PASA...35...31P>`_, `DOI <https://doi.org/10.1017/pasa.2018.25>`_]

============
Contributors
============

Author:

* `Daniel Mentiplay <https://github.com/dmentipl>`_

=======
License
=======

Copyright 2019 Daniel Mentiplay and contributors.

Plonk is available under the MIT license. For details see the `LICENSE <https://github.com/dmentipl/plonk/blob/master/LICENSE>`_ file.

Splash is available at https://github.com/danieljprice/splash.
