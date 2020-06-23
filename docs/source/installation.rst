============
Installation
============

.. note::
    The preferred way to install Plonk is with Conda.

-----
Conda
-----

You can install Plonk via the package manager Conda from conda-forge.

.. code-block:: console

    $ conda install plonk

This will install the required dependencies. Note: you may need to first add the
`conda-forge` channel with :code:`conda config --add channels conda-forge`. I
also recommend strictly using conda-forge which you can do with :code:`conda
config --set channel_priority true`. For details on Conda, see
`<https://docs.conda.io/>`_.


---
pip
---

You can also install Plonk via pip.

.. code-block:: console

    $ pip install plonk

This should install the required dependencies. For details on pip, see
`<https://pip.pypa.io/>`_.

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
