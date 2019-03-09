'''
analyze_disc.py

Daniel Mentiplay, 2019.
'''

from plonk.plonk.analysis.disc import disc_analysis
from plonk.plonk.dump import Dump

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

    print('\nPerforming disc analysis...')
    results = disc_analysis( radius_in          = radius_in,
                             radius_out         = radius_out,
                             number_radial_bins = number_radial_bins,
                             particles          = dump.particles,
                             sinks              = dump.sinks,
                             parameters         = dump.parameters,
                             units              = dump.units )

    radial_averages.append(results[0])
    particles.      append(results[1])
    sinks.          append(results[2])


message = '''
Variables available:
  dumps           : list of Dump objects from dump files read in
  radial_averages : list of DataFrames with radial averages over the disc from
                    each dump file
  particles       : list of DataFrames with quantities on particles
  sinks           : list of DataFrames with quantities on sinks
'''

print(message)
