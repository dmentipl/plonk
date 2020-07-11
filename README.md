Plonk
=====

Smoothed particle hydrodynamics analysis and visualization with Python.

+ Docs: <https://plonk.readthedocs.io/>
+ Repo: <https://www.github.com/dmentipl/plonk>

[![Build Status](https://travis-ci.org/dmentipl/plonk.svg?branch=master)](https://travis-ci.org/dmentipl/plonk)
[![Coverage Status](https://coveralls.io/repos/github/dmentipl/plonk/badge.svg?branch=master)](https://coveralls.io/github/dmentipl/plonk?branch=master)
[![Documentation Status](https://readthedocs.org/projects/plonk/badge/?version=stable)](https://plonk.readthedocs.io/en/stable/?badge=stable)

[![PyPI](https://img.shields.io/pypi/v/plonk)](https://pypi.org/project/plonk/)
[![Anaconda Version](https://img.shields.io/conda/v/conda-forge/plonk.svg)](https://anaconda.org/conda-forge/plonk)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/dmentipl/plonk/blob/master/LICENSE)

[![JOSS](https://joss.theoj.org/papers/10.21105/joss.01884/status.svg)](https://doi.org/10.21105/joss.01884)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.3698382.svg)](https://doi.org/10.5281/zenodo.3698382)

Description
-----------

Plonk is a Python tool for analysis and visualization of smoothed particle hydrodynamics data with a focus on astrophysical fluid dynamics.

With Plonk we aim to integrate the high quality SPH visualisation of [Splash](https://github.com/danieljprice/splash) into the modern Python astronomer workflow, and provide a framework for analysis of smoothed particle hydrodynamics simulation data.

Usage
-----

Plonk supports the following SPH file formats:

+ [Phantom](https://github.com/danieljprice/phantom) output in [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) format.

*Note: you can convert Phantom non-HDF5 snapshots to HDF5. See the [Phantom docs](https://phantomsph.readthedocs.io).*

### Accessing data

To read in a simulation with snapshot files like `disc_00000.h5`, and global quantity time series files like `disc01.ev`, in the current directory, and see what snapshots there are:

```python
>>> import plonk

>>> simulation = plonk.load_sim(prefix='disc')
>>> simulation.snaps
[<plonk.Snap "disc_00000.h5">,
 ...
 <plonk.Snap "disc_00030.h5">]
```

You can load individual snapshots and access the particle arrays:

```python
>>> snap = plonk.load_snap('disc_00030.h5')
>>> snap['position']
array([[ -24.69953214,   49.60113417,   -4.98059478],
       [-108.99243136,   77.74663833,   12.89299546],
       [ -51.22218782,  108.64454019,    1.56619644],
       ...,
       [  93.296599  ,  -77.66042087,    5.40835798],
       [  63.75108128,   66.7446782 ,    3.30169363],
       [   8.11639008,  139.45117413,    7.55340187]])
```

The Snap objects contain the particle arrays, lazily loaded from the HDF5 file, as well as simulation metadata properties stored as a dictionary.

### Visualization

To visualize the column density on a snapshot:

```python
>>> plonk.visualize.plot(snap=snap, quantity='density')
```

For a more complicated example, here is the deviation from Keplerian velocity around a planet embedded in a protoplanetary disc.

![Planet embedded in protoplanetary disc](https://raw.githubusercontent.com/dmentipl/plonk/master/image.png)

*Deviation from Keplerian velocity around a planet: at the disc midplane (left), and 10 (middle) and 20 au (right) above the disc midplane. Data from a Phantom simulation.*

### Analysis

Extra quantities not written to the snapshot file are available:

```python
>>> snap.extra_quantities()
<plonk.Snap "disc_00030.h5">

>>> snap['angular_momentum']
array([ ... ])
```

You can generate radial profiles on the snapshot. For example, to calculate the scale height in a disc:

```python
>>> prof = plonk.load_profile(snap)

>>> prof['scale_height']
array([ ... ])
```

Physical units for array quantities and other properties are available.

```python
>>> snap['position'][0]
array([-24.69953214,  49.60113417,  -4.98059478])

>>> snap.physical_units()
<plonk.Snap "disc_00030.h5">

>>> snap['position'][0]
array([-3.69505001e+14,  7.42032967e+14, -7.45096980e+13]) <Unit('centimeter')>

>>> snap['position'][0].to('au')
array([-24.6998837 ,  49.60184016,  -4.98066567]) <Unit('astronomical_unit')>
```

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

This should install the required dependencies. For details on pip, see
<https://pip.pypa.io/>.

Requirements
------------

Python 3.6+ with h5py, matplotlib, numba, numpy, pandas, pint, scikit-image,
scipy, tqdm. Installing Plonk with conda or pip will install these dependencies.

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
