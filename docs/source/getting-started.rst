===============
Getting started
===============

-----------------
Data file formats
-----------------

Plonk supports the following SPH file formats:

* Phantom output in
  `HDF <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_
  form (as opposed to the sphNG-based Fortran binary format).

.. note::
    HDF5 output was added to Phantom as an option with git commit
    `9b22ded <https://bitbucket.org/danielprice/phantom/commits/9b22ded9e7b4d512966f2b2e4b84d693b1afc9e6>`_
    on the 14th of March 2019. See the `Phantom documentation
    <https://phantomsph.readthedocs.io/>`_ for instructions on
    how to compile with HDF5 output and to convert from the sphNG-based
    output.

---------------------
Working with SPH data
---------------------

First import the Plonk package.

.. code-block:: pycon

    >>> import plonk

SPH snapshot files are represented by the `Snap` class. This object contains a
properties dictionary, particle arrays, which are lazily loaded from file. Here
we demonstrate instantiating a `Snap` object, and accessing some properties and
particle arrays.

.. code-block:: pycon

    >>> filename = 'snap_99999.h5'
    >>> snap = plonk.load_snap(filename)
    >>> snap['position']
    >>> snap.properties['time']


In a single SPH simulation the data can be spread over multiple files of
multiple types, even though, logically, a simulation is a singular "object". In
Plonk we have the `Simulation` class to represent the totality of the simulation
data. It is an aggregation of the `Snap` and `Evolution` objects, plus metadata,
such as the directory on the file system.

.. code-block:: pycon

    >>> prefix = 'mysim'
    >>> sim = plonk.load_sim(prefix=prefix)

In the above code example we assume the files are like `mysim_00000.h5`, and
`mysim01.ev`, and so on.

The `Visualization` class provides methods to visualise SPH data. Plonk uses
KDEpy for interpolation to a pixel array using kernel density estimation. This
is equivalent to the interpolation routines in Splash. Instantiation of a
`Visualization` object produces a figure.

.. code-block:: pycon

    >>> viz = plonk.visualize.render(
    ...    snap=snap,
    ...    quantity='density',
    ... )
