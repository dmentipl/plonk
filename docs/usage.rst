===========
Using Plonk
===========

Example scripts and IPython notebooks are found in the examples folder.

-----------
Basic usage
-----------

Say you have 50 dump files in a directory data. The following code loads data from each file into a Dump object and puts each object into a list::

    from plonk.dump import Dump

    files = [f'data/disc_{idx:05}.h5' for idx in range(50)]

    dumps = list()
    for file in files:
        dumps.append(Dump(dump_file_name))
