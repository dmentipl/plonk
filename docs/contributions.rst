=============
Contributions
=============

.. note:: Contributions are welcome.

If you want to contribute to Plonk you can fork the repository, clone it, and use Conda to link to your local copy of the code.

::

 git clone https://github.com/<user>/plonk
 cd plonk && conda develop .

There is a compiled Fortran component to Plonk which is derived from Splash. You must compile this before development. This requires a Fortran compiler, e.g. gfortran. The following compiles Splash into a shared object library, and then uses Cython to build a Python interface to that library.

::

 make install
 python setup.py build_ext --inplace

You need to make sure the required dependencies are installed (via Conda). To satisfy these requirements there is a ``environment.yml`` file. You can set up a Conda environment for development and install Plonk in it.

::

 git clone https://github.com/<user>/plonk && cd plonk
 conda env create --name plonk_dev --file environment.yml
 conda activate plonk_dev

and then follow the instructions above. (To leave the development environment: ``conda deactivate``.)

After you have committed and pushed your changes to your forked repository you can issue a `pull request <https://github.com/dmentipl/plonk/pull/new/master>`_.
