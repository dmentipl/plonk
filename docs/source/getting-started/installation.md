# Installation

## Conda

You can install Plonk via the package manager [Conda](https://docs.conda.io/)
from the [conda-forge](https://conda-forge.org/) channel.

```console
conda install plonk --channel conda-forge
```

This will install the required dependencies.

```{note}
You can simply use `conda install plonk` if you add the `conda-forge`
channel with `conda config --add channels conda-forge`. I also
recommend strictly using conda-forge which you can do with `conda
config --set channel_priority true`. Both of these commands modify the Conda
configuration file `~/.condarc`.
```

## pip

You can also install Plonk from [PyPI](https://pypi.org/) via [pip](https://pip.pypa.io/).

```console
python -m pip install plonk
```

This should install the required dependencies.

## Source

You can install Plonk from source as follows.

```console
# clone via HTTPS
$ git clone https://github.com/dmentipl/plonk.git

# or clone via SSH
$ git clone git@github.com:dmentipl/plonk

$ cd plonk
$ python -m pip install -e .
```

This assumes you have already installed the dependencies. One way to do this is
by setting up a conda environment. The [environment.yml](https://github.com/dmentipl/plonk/blob/main/environment.yml) file provided
sets up a conda environment "plonk" for using or developing Plonk.

```console
conda env create --file environment.yml
conda activate plonk
```

Plonk has Python runtime requirements listed in [setup.cfg](https://github.com/dmentipl/plonk/blob/main/setup.cfg) in the
`install_requires` variable.
