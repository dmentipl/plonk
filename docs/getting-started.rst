===============
Getting started
===============

-----------------
Data file formats
-----------------

Plonk supports the following SPH file formats:

* Phantom output in `HDF <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ form (as opposed to the sphNG-based Fortran binary format).

.. note:: HDF output was added to Phantom as an option with git commit `9b22ded <https://bitbucket.org/danielprice/phantom/commits/9b22ded9e7b4d512966f2b2e4b84d693b1afc9e6>`_ on the 14th of March 2019. See the `Phantom wiki <https://bitbucket.org/danielprice/phantom/wiki/Home>`_ for instructions on how to compile with HDF output and to convert from the Fortran binary output.

---------------------
Working with SPH data
---------------------

First import the Plonk package.

.. code-block:: pycon

 >>> import plonk

SPH dump files are represented by the `Dump` class. This object contains a header dictionary, particle arrays, which are represented by an `Arrays` class.  Here we demonstrate instantiating a `Dump` object, and accessing the header and particle arrays.

.. code-block:: pycon

 >>> filename = 'dump_99999.h5'
 >>> dump = plonk.Dump(filename)


In a single SPH simulation the data can be spread over multiple files of multiple types, even though, logically, a simulation is a singular "object". In Plonk we have the `Simulation` class to represent the totality of the simulation data. It is an aggregation of the `Dump` and `Evolution` objects, plus metadata, such as the directory on the file system.

.. code-block:: pycon

 >>> prefix = 'mysim'
 >>> sim = plonk.Simulation(prefix=prefix).

In the above code example we assume the files are like `mysim_00000.h5`, and `mysim01.ev`, and so on.

The `Visualization` class provides methods to visualise SPH data. Plonk uses Splash for interpolation to a pixel array. Instantiation of a `Visualization` object produces a figure.

.. code-block:: pycon

 >>> viz = plonk.Visualization(
 ...    dump=dump,
 ...    render='density',
 ... )
