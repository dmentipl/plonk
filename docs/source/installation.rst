============
Installation
============

.. note::
    The preferred way to install Plonk is with Conda.

-----
Conda
-----

Plonk can be installed with Conda.

.. code-block:: console

    $ conda install plonk

The Conda package is available via conda-forge at
`<https://anaconda.org/conda-forge/plonk/>`_.

If you're unfamiliar with Conda, see `<https://docs.conda.io/>`_.

---
pip
---

Plonk can also be installed with pip.

.. code-block:: console

    $ pip install plonk

If you're unfamiliar with pip, see `<https://pip.pypa.io/>`_.

Both of these installation methods should install all required Python runtime
packages.

------
Source
------

You can install Plonk from source as follows.

.. code-block:: console

    $ git clone https://github.com/dmentipl/plonk.git
    $ cd plonk
    $ conda env create --file environment.yml
    $ conda activate plonk-dev
    $ python setup.py install

This method requires using Conda to install the Python dependencies. Otherwise,
if you have already satisfied the dependencies, just cloning the repository and
running the `setup.py` script should work.

------------
Requirements
------------

Plonk has runtime Python requirements listed in `setup.py
<https://github.com/dmentipl/plonk/blob/master/setup.py>`_.

To set up a Conda environment for development of Plonk see
`environment.yml
<https://github.com/dmentipl/plonk/blob/master/environment.yml>`_.
