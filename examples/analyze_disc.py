'''
analyze_disc.py

Daniel Mentiplay, 2019.
'''

from plonk.analysis.disc import disc_analysis
from plonk.dump import Dump

#--- Options

# TODO: get as input
number_radial_bins = 150
radius_in          = 10
radius_out         = 200

#--- Dump file names

# TODO: get dump file names as input
dump_file_names = ['data/disc_00000.h5']

# ---------------------------------------------------------------------------- #

dumps = list()

for dump_file_name in dump_file_names:

#--- Read dump file

    print('\nReading in data from dumpfile: ' + dump_file_name + '...')

    dump = Dump(dump_file_name)

    dumps.append(dump)

#--- Perform analysis

radial_averages = list()
particles       = list()
sinks           = list()

for dump in dumps:

    print('\nPerforming disc analysis...\n')
    radial_averages_ = disc_analysis( radius_in          = radius_in,
                                      radius_out         = radius_out,
                                      number_radial_bins = number_radial_bins,
                                      dump               = dump )

    radial_averages.append(radial_averages_)
