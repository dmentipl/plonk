==============
Plonk overview
==============

Plonk is a Python tool for analysis and visualization of Phantom data. It also contains Fortran code from Splash as part of the visualization module.

Plonk aims to provide an alternative to PhantomAnalysis and Splash. PhantomAnalysis, which is distributed with Phantom itself, and Splash are both written in Fortran, and thus, as compiled binary programs, are somewhat less flexible and scriptable than an interpreted Python program.

-------
Phantom
-------

Plonk relies on Phantom output in `HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ form (as opposed to the earlier Fortran binary output). HDF5 output was added to Phantom as an option with git commit `9b22ded <https://bitbucket.org/danielprice/phantom/commits/9b22ded9e7b4d512966f2b2e4b84d693b1afc9e6>`_ on the 14th of March 2019. See the `Phantom wiki <https://bitbucket.org/danielprice/phantom/wiki/Home>`_ for instructions on how to compile with HDF5 output and to convert earlier output.

------
Splash
------

Plonk *does not* rely on having Splash installed. The required Fortran code is included in the Plonk code repository.

-----
Links
-----

* Phantom: https://phantomsph.bitbucket.io/
* Splash: http://users.monash.edu.au/~dprice/splash/
