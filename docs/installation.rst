============
Installation
============

The easiest and recommended way to install Plonk is via the package manager Conda::

 conda install -c dmentipl plonk

This will install the required dependencies. For details on Conda, see https://docs.conda.io/.

------------
Requirements
------------

Plonk has Python requirements listed in `environment.yml <https://github.com/dmentipl/plonk/blob/master/environment.yml>`_ located in the git repository. These requirements are satisfied by Conda.

.. note:: Plonk relies on Phantom output in `HDF <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ form (as opposed to the earlier Fortran binary output).

HDF output was added to Phantom as an option with git commit `9b22ded <https://bitbucket.org/danielprice/phantom/commits/9b22ded9e7b4d512966f2b2e4b84d693b1afc9e6>`_ on the 14th of March 2019. See the `Phantom wiki <https://bitbucket.org/danielprice/phantom/wiki/Home>`_ for instructions on how to compile with HDF output and to convert earlier output.
