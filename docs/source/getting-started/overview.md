# Overview

```{eval-rst}
.. currentmodule:: plonk
```

This document gives an overview of using Plonk for analysis and visualization of
smoothed particle hydrodynamics data. For a further guide see {doc}`../../user-guide/usage`.

```{toctree}
:maxdepth: 1

overview/data
overview/visualization
overview/analysis
```

```{important}
To follow along, download the sample data `plonk_example_data.tar`
from [figshare](https://figshare.com/articles/dataset/Plonk_example_dataset/12885587).
Then extract with `tar xvf plonk_example_data.tar` and change into the
`plonk_example_data` directory. This data set is from a Phantom
simulation of a dust and gas protoplanetary disc with an embedded
protoplanet.
```

## Data file formats

Plonk supports the following SPH file formats:

- Phantom output in
  [HDF](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)
  form (as opposed to the sphNG-based Fortran binary format).

```{note}
HDF5 output was added to Phantom as an option with git commit
[9b22ded](https://github.com/danieljprice/phantom/commit/9b22ded9e7b4d512966f2b2e4b84d693b1afc9e6)
on the 14th of March 2019. See the [Phantom documentation](https://phantomsph.readthedocs.io/) for instructions on
how to compile with HDF5 output and to convert from the sphNG-based
output.
```
