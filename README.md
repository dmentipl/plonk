Plonk
=====

Phantom analysis and visualization but with Python.

+ Docs: https://plonk.readthedocs.io/
+ Repo: https://www.github.com/dmentipl/plonk

[![Build Status](https://travis-ci.com/dmentipl/plonk.svg?token=AL8sPDxCNprS78nBjSQh&branch=master)](https://travis-ci.com/dmentipl/plonk)

Usage
-----

Plonk requires Phantom output to be in HDF5 format. See the Phantom wiki for more information: https://bitbucket.org/danielprice/phantom/wiki/Home.

### Basic usage

To read in a collection of Phantom dump files:

```python
from plonk.dump import Dump

files = [f'data/disc_{idx:05}.h5' for idx in range(50)]

dumps = list()
for file in files:
    dumps.append(Dump(dump_file_name))
```

For further usage examples see `examples` folder.

Install
-------

To install Plonk:

```
python setup.py install
```

### Requirements

Plonk has Python requirements listed in `requirements.txt`. The recommended way to satisfy these requirements is to use Anaconda (https://anaconda.org/).

In addition, Plonk requires a recent (2003+) version of the GCC Fortran compiler gfortran (https://gcc.gnu.org/wiki/GFortran) to compile the Splash subroutines.
