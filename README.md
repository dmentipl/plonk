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
array([[-3.69505001e+14,  7.42032967e+14, -7.45096980e+13],
       [-1.63052677e+15,  1.16308971e+15,  1.92879212e+14],
       [-7.66283930e+14,  1.62532232e+15,  2.34302988e+13],
       ...,
       [ 1.39571712e+15, -1.16179990e+15,  8.09090354e+13],
       [ 9.53716176e+14,  9.98500386e+14,  4.93933367e+13],
       [ 1.21421196e+14,  2.08618956e+15,  1.12998892e+14]]) <Unit('centimeter')>
```

The Snap objects contain the particle arrays, lazily loaded from the HDF5 file, as well as simulation metadata properties stored as a dictionary.

### Visualization

To visualize the column density on a snapshot:

```python
>>> snap.plot(quantity='density')
```

For a more complicated example, here is the deviation from Keplerian velocity around a planet embedded in a protoplanetary disc.

![Planet embedded in protoplanetary disc](https://raw.githubusercontent.com/dmentipl/plonk/master/image.png)

*Deviation from Keplerian velocity around a planet: at the disc midplane (left), and 10 (middle) and 20 au (right) above the disc midplane. See [here](https://plonk.readthedocs.io/en/latest/examples/deviation-from-keplerian.html) for details.*

### Analysis

Extra quantities not written to the snapshot file are available:

```python
>>> snap['angular_momentum']
array([ ... ])
```

You can generate radial profiles on the snapshot. For example, to calculate the scale height in a disc:

```python
>>> prof = plonk.load_profile(snap)

>>> prof['scale_height']
array([ ... ])
```

Physical units of array quantities and other properties allow for unit conversion:

```python
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

You can install Plonk via the package manager [Conda](https://docs.conda.io/) from [conda-forge](https://conda-forge.org/).

```bash
conda install plonk --channel conda-forge
```

This will install the required dependencies.

Note: You can simply use `conda install plonk` if you add the `conda-forge` channel with `conda config --add channels conda-forge`. I also recommend strictly using conda-forge which you can do with `conda config --set channel_priority true`. Both of these commands modify the Conda configuration file `~/.condarc`.

### pip

You can also install Plonk from [PyPI](https://pypi.org/) via [pip](https://pip.pypa.io/).

```bash
python -m pip install plonk
```

This should install the required dependencies.

### Source

You can install Plonk from source as follows.

```bash
# clone via HTTPS
git clone https://github.com/dmentipl/plonk.git

# or clone via SSH
git clone git@github.com:dmentipl/plonk

cd plonk
python -m pip install -e .
```

This assumes you have already installed the dependencies. One way to do this is by setting up a conda environment. The [environment.yml](https://github.com/dmentipl/plonk/blob/master/environment.yml) file provided sets up a conda environment "plonk" for using or developing Plonk.

```bash
conda env create --file environment.yml
conda activate plonk
```

Requirements
------------

Python 3.6+ with [h5py](https://www.h5py.org/), [matplotlib](https://www.matplotlib.org/), [numba](http://numba.pydata.org/), [numpy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [pint](https://pint.readthedocs.io/), [scipy](https://www.scipy.org/). Installing Plonk with conda or pip will install these dependencies.

Getting help
------------

If you need help, please try the following:

1. Check the [documentation](https://plonk.readthedocs.io/).
2. Check the built-in help, e.g. `help(plonk.load_snap)`.
3. Raise an issue, as a [bug report](https://github.com/dmentipl/plonk/issues/new?assignees=&labels=&template=bug_report.md&title=) or [feature request](https://github.com/dmentipl/plonk/issues/new?assignees=&labels=&template=feature_request.md&title=), using the issue tracker.

Please don't be afraid to raise an issue. Even if your issue is not a bug, it's nice to have the question and answer available in a public forum so we can all learn from it together.

If you don't get an immediate response, please be patient. Plonk is maintained by one person, [@dmentipl](https://github.com/dmentipl).

Contributing
------------

*All types of contributions are welcome from all types of people with different skill levels.*

Thank you for considering contributing to Plonk. There are many ways to contribute:

1. If you find any bugs or cannot work out how to do something, please file a [bug report](https://github.com/dmentipl/plonk/issues/new?assignees=&labels=&template=bug_report.md&title=) in the issue tracker. Even if the issue is not a bug it may be that there is a lack of documentation.
2. If you have any suggestions for new features, please raise a [feature request](https://github.com/dmentipl/plonk/issues/new?assignees=&labels=&template=feature_request.md&title=) in the issue tracker.
3. If you use Plonk to do anything please consider contributing to the [examples](https://plonk.readthedocs.io/en/stable/examples.html) section in particular, or any other section, of the documentation.
4. If you would like to contribute code, firstly thank you! We take code contributions via [pull request](https://github.com/dmentipl/plonk/pull/new/master).

See [CONTRIBUTING.md](https://github.com/dmentipl/plonk/blob/master/CONTRIBUTING.md) for detailed guidelines on how to contribute.

Citation
--------

If you use Plonk in a scientific publication, please cite the paper published in JOSS.

> [Plonk: Smoothed particle hydrodynamics analysis and visualization with Python](https://doi.org/10.21105/joss.01884)

A BibTeX entry is available in [CITATION.bib](https://github.com/dmentipl/plonk/blob/master/CITATION.bib)

If you use the interpolation to pixel grid component of Plonk please cite the [Splash paper](https://doi.org/10.1071/AS07022). You should also consider citing any other scientific software packages that you use.

Change log
----------

The change log is available in [CHANGELOG.md](https://github.com/dmentipl/plonk/blob/master/CHANGELOG.md)
