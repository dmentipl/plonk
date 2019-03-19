Plonk
=====

Phantom analysis and visualization but with Python.

+ Docs: https://plonk.readthedocs.io/
+ Repo: https://www.github.com/dmentipl/plonk

[![Build Status](https://travis-ci.com/dmentipl/plonk.svg?token=AL8sPDxCNprS78nBjSQh&branch=master)](https://travis-ci.com/dmentipl/plonk)

Usage
-----

Plonk requires Phantom output to be in HDF format.

*You can convert old dumps to HDF.* See the Phantom wiki for more information (https://bitbucket.org/danielprice/phantom/wiki).

### Basic usage

To read in a collection of Phantom dump files with names like `disc_00000.h5`, ...

```python
from plonk.dump import Dump

n_files = 100
files = [f'disc_{idx:05}.h5' for idx in range(n_files)]

dumps = list()
for file in files:
    dumps.append(Dump(dump_file_name))
```

For further usage, see `examples` folder.

Install
-------

*Note that Plonk is a Python 3 only package.*

### Conda

The easiest and recommended way to install Plonk is via the package manager Conda

```
conda install -c dmentipl plonk
```

This will install the required dependencies. For details on Conda, see https://docs.conda.io/.

Contribute
----------

*Any contributions are welcome.*

If you want to contribute to Plonk you can fork the repository, clone it, and use Conda to link to your local copy of the code.

```
git clone https://github.com/<user>/plonk
cd plonk && conda develop .
```

You need to make sure the required dependencies are installed (via Conda). To satisfy these requirements there is a `environment.yml` file. You can set up a Conda environment for development and install Plonk in it:

```
git clone https://github.com/<user>/plonk && cd plonk
conda env create --name plonk_dev --file environment.yml
conda activate plonk_dev
conda develop .
```

After you have committed and pushed your changes to your forked repository you
can issue a pull request: https://github.com/dmentipl/plonk/pull/new/master.
