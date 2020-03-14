Plonk
=====

Smoothed particle hydrodynamics analysis and visualization with Python.

+ Docs: <https://plonk.readthedocs.io/>
+ Repo: <https://www.github.com/dmentipl/plonk>

[![Build Status](https://travis-ci.org/dmentipl/plonk.svg?branch=master)](https://travis-ci.org/dmentipl/plonk)
[![Coverage Status](https://coveralls.io/repos/github/dmentipl/plonk/badge.svg?branch=master)](https://coveralls.io/github/dmentipl/plonk?branch=master)
[![Documentation Status](https://readthedocs.org/projects/plonk/badge/?version=latest)](https://plonk.readthedocs.io/en/latest/?badge=latest)

[![PyPI](https://img.shields.io/pypi/v/plonk)](https://pypi.org/project/plonk/)
[![Anaconda Version](https://img.shields.io/conda/v/conda-forge/plonk.svg)](https://anaconda.org/conda-forge/plonk)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/dmentipl/plonk/blob/master/LICENSE)

[![JOSS](https://joss.theoj.org/papers/10.21105/joss.01884/status.svg)](https://doi.org/10.21105/joss.01884)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.3698382.svg)](https://doi.org/10.5281/zenodo.3698382)


Description
-----------

Plonk is a Python tool for analysis and visualization of smoothed particle hydrodynamics data with a focus on astrophysical fluid dynamics.

With Plonk we aim to integrate the high quality SPH visualisation of Splash into the modern Python astronomer workflow, and provide a framework for analysis of smoothed particle hydrodynamics simulation data.

Usage
-----

Plonk supports the following SPH file formats:

+ [Phantom](https://phantomsph.bitbucket.io/) output in [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) format.

*Note: you can convert Phantom non-HDF5 snapshots to HDF5. See the [Phantom docs](https://phantomsph.readthedocs.io).*

### Accessing data

To read in a simulation with snapshot files like `disc_00000.h5`, and global quantity time series files like `disc01.ev`, in the current directory, and see what snapshots there are:

```python
>>> import plonk

>>> simulation = plonk.load_sim(prefix='disc')
>>> simulation.snaps
[<plonk.Snap: "disc_00000.h5">,
 ...
 <plonk.Snap: "disc_01000.h5">]
```

The Snap objects contain the particle arrays, lazily loaded from the HDF5 file, as well as simulation metadata properties stored as a dictionary.

### Visualization

To render the density on a snapshot:

```python
>>> snap = plonk.load_snap('disc_00000.h5')
>>> plonk.visualize.plot(snap=snap, quantity='density')
```

For a more complicated example, here is the deviation from Keplerian velocity around a planet embedded in a protoplanetary disc.

![Planet embedded in protoplanetary disc](https://raw.githubusercontent.com/dmentipl/plonk/master/image.png)

*Deviation from Keplerian velocity around a planet: at the disc midplane (left), and 10 (middle) and 20 au (right) above the disc midplane. Data from a Phantom simulation.*

### More

For further usage, see documentation. The code is internally documented with docstrings. Try, for example, `help(plonk.Snap)` or `help(plonk.load_snap)`.

Install
-------

### Conda

You can install Plonk via the package manager Conda from conda-forge.

```bash
conda install plonk
```

This will install the required dependencies. Note: you may need to first add the `conda-forge` channel with `conda config --add channels conda-forge`. I also recommend strictly using conda-forge which you can do with `conda config --set channel_priority true`. For details on Conda, see <https://docs.conda.io/>.

### pip

You can also install Plonk via pip.

```bash
pip install plonk
```

This should install the required dependencies. For details on pip, see <https://pip.pypa.io/>.

Getting help
------------

If you need help, please try the following, in order:

1. Check the [documentation](https://plonk.readthedocs.io/).
2. Check the built-in help, e.g. `help(plonk.load_snap)`.
3. Raise an issue, as a [bug report](https://github.com/dmentipl/plonk/issues/new?assignees=&labels=&template=bug_report.md&title=) or [feature request](https://github.com/dmentipl/plonk/issues/new?assignees=&labels=&template=feature_request.md&title=), using the issue tracker.

Please don't be afraid to raise an issue. Even if your issue is not a bug, it's nice to have the question and answer available in a public forum so we can all learn from it together.

If you don't get an immediate response, please be patient. Plonk is maintained by one person, [@dmentipl](https://github.com/dmentipl).

Contributing
------------

Thank you for considering contributing to Plonk. *All types of contributions are welcome from all types of people with different skill levels.*

See [CONTRIBUTING.md](https://github.com/dmentipl/plonk/blob/master/CONTRIBUTING.md) for guidelines on how to contribute.

Citation
--------

If you use Plonk in a scientific publication, please cite the paper published in JOSS.

> [Plonk: Smoothed particle hydrodynamics analysis and visualization with Python](https://joss.theoj.org/papers/10.21105/joss.01884#)

You should also consider citing any other scientific software packages that you use.
