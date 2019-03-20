.. Plonk documentation master file, created by
   sphinx-quickstart on Mon Mar 18 15:51:00 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====
Plonk
=====

**Phantom analysis and visualization but with Python.**

* Docs: https://plonk.readthedocs.io/
* Repo: https://www.github.com/dmentipl/plonk

--------
Overview
--------

Plonk is a Python tool for analysis and visualization of Phantom data.

Plonk aims to provide an alternative to PhantomAnalysis and Splash. PhantomAnalysis, which is distributed with Phantom itself, and Splash are both written in Fortran, and thus, as compiled binary programs, are somewhat less flexible and scriptable than an interpreted Python program.

Plonk contains Fortran interpolation subroutines from Splash as part of the visualization module.

.. warning:: Plonk is under development. It may contain bugs. The authors do not take responsibility for use of Plonk.

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   installation
   getting-started
   examples
   contributions

.. role:: math(raw)
   :format: html latex
..

----------
References
----------

* Splash:
   * http://users.monash.edu.au/~dprice/splash/
   * Price D. J., 2007, PASA, 24, 159 [`ADS <http://adsabs.harvard.edu/abs/2007PASA...24..159P>`_, `DOI <http://dx.doi.org/10.1071/AS07022>`_]
* Phantom:
   * https://phantomsph.bitbucket.io/
   * Price D. J., et al., 2018, PASA, 35, 31P [`ADS <http://adsabs.harvard.edu/abs/2018PASA...35...31P>`_, `DOI <https://doi.org/10.1017/pasa.2018.25>`_]

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
