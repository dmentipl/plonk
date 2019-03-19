==========
Contribute
==========

*Any contributions are welcome.*

If you want to contribute to Plonk you can fork the repository, clone it, and use Conda to link to your local copy of the code::

 git clone https://github.com/<user>/plonk
 cd plonk && conda develop .

There is also a compiled Fortran component to Plonk which is derived from Splash. This needs to be compiled during development::

 python setup.py build_ext --inplace

You need to make sure the required dependencies are installed (via Conda). To satisfy these requirements there is a ``environment.yml`` file. You can set up a Conda environment for development and install Plonk in it::

 git clone https://github.com/<user>/plonk && cd plonk
 conda env create --name plonk_dev --file environment.yml
 conda activate plonk_dev
 conda develop .
 python setup.py build_ext --inplace

After you have committed and pushed your changes to your forked repository you
can issue a pull request: https://github.com/dmentipl/plonk/pull/new/master.
