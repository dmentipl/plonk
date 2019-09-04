Plonk
=====

Smoothed particle hydrodynamics analysis and visualization with Python.

+ Docs: https://plonk.readthedocs.io/
+ Repo: https://www.github.com/dmentipl/plonk

[![Build Status](https://travis-ci.org/dmentipl/plonk.svg?branch=master)](https://travis-ci.org/dmentipl/plonk)
[![Coverage Status](https://coveralls.io/repos/github/dmentipl/plonk/badge.svg?branch=master)](https://coveralls.io/github/dmentipl/plonk?branch=master)
[![Documentation Status](https://readthedocs.org/projects/plonk/badge/?version=latest)](https://plonk.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/dmentipl/plonk/blob/master/LICENSE)

[![Anaconda Version](https://img.shields.io/conda/v/dmentipl/plonk.svg)](https://anaconda.org/dmentipl/plonk)
[![Anaconda Platform](https://img.shields.io/conda/pn/dmentipl/plonk.svg)](https://anaconda.org/dmentipl/plonk)
[![ASCL](https://img.shields.io/badge/ascl-1907.009-blue.svg?colorB=262255)](http://ascl.net/1907.009)

Usage
-----

Plonk supports the following SPH file formats:

* [Phantom](https://phantomsph.bitbucket.io/) output in [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) format.

*Note: you can convert Phantom standard dumps to HDF5. See the [Phantom docs](https://phantomsph.readthedocs.io).*

### Accessing data

To read in a simulation with dump files like `disc_00000.h5`, ..., and evolution files like `disc01.ev`, ..., in the current directory, and see what dumps there are:

```python
>>> import plonk

>>> simulation = plonk.Simulation(prefix='disc')
>>> simulation.dumps
[<plonk.Dump: "disc_00000.h5">,
 ...
 <plonk.Dump: "disc_01000.h5">]
```

The Dump objects contain the particle and sinks arrays, lazily loaded from the HDF5 file, as well as the dump header stored as a dictionary.

### Visualization

To visualize a single dump file:

```python
>>> dump = plonk.Dump('disc_00000.h5')

>>> plonk.Visualization(
...     dump=dump,
...     render='density',
...     extent=[-200, 200, -200, 200]
...     )
```

For example, here is the deviation from Keplerian velocity around a planet embedded in a protoplanetary disc.

![](image.svg)

*Deviation from Keplerian velocity around a planet: at the disc midplane (left), and 10 (middle) and 20 au (right) above the disc midplane. Data from a Phantom simulation.*

### More

For further usage, see `examples` folder and documentation. The code is internally documented with docstrings. Try, for example, `help(plonk.Dump)` or `help(plonk.Visualization)`.

Install
-------

*Plonk is a Python 3 only package.*

### Conda

The easiest and recommended way to install Plonk is via the package manager Conda

```bash
conda install plonk --channel dmentipl
```

*Note*: Using this method you don't need to have this repository on your machine.

This will install the required dependencies. For details on Conda, see https://docs.conda.io/.

Getting help
------------

If you need help, try the following, in order:

1. Check the documentation.
2. Ask questions on Stack Overflow using the [plonk](https://stackoverflow.com/questions/tagged/plonk) tag.
3. File an issue, as a [bug report](https://github.com/dmentipl/plonk/issues/new?assignees=&labels=&template=bug_report.md&title=) or [feature request](https://github.com/dmentipl/plonk/issues/new?assignees=&labels=&template=feature_request.md&title=), using the issue tracker.

If you don't get an immediate response, please be patient. Plonk is run by one person, [@dmentipl](https://github.com/dmentipl).

Contributing
------------

Thank you for considering contributing to Plonk. *Contributions are welcome.*

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to contribute.

Citation
--------

If you use Plonk in a scientific publication, please cite

> [Plonk: Smoothed particle hydrodynamics data analysis and visualization](https://ui.adsabs.harvard.edu/abs/2019ascl.soft07009M/abstract)

You should also consider citing any other scientific software packages that you use.
