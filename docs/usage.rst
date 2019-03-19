===========
Using Plonk
===========

Example scripts and IPython notebooks are found in the examples folder.

-------
Phantom
-------

Plonk relies on Phantom output in `HDF <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ form (as opposed to the earlier Fortran binary output). HDF output was added to Phantom as an option with git commit `9b22ded <https://bitbucket.org/danielprice/phantom/commits/9b22ded9e7b4d512966f2b2e4b84d693b1afc9e6>`_ on the 14th of March 2019. See the `Phantom wiki <https://bitbucket.org/danielprice/phantom/wiki/Home>`_ for instructions on how to compile with HDF output and to convert earlier output.

-----------
Basic usage
-----------

Say you have 50 dump files labelled ``disc_00000.h5, ...`` in the current directory. The following code loads data from each file into a Plonk Dump object and puts each object into a list::

 from plonk.dump import Dump

 files = [f'disc_{idx:05}.h5' for idx in range(50)]

 dumps = list()
 for file in files:
     dumps.append(Dump(dump_file_name))
