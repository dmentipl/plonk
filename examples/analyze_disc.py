'''
analyze_disc.py

Daniel Mentiplay, 2019.
'''

from plonk.plonk.analysis.disc import disc_analysis
from plonk.plonk.PhantomDump import PhantomDump

#--- Options

# TODO: get as input
numberRadialBins = 150
radiusIn         = 10
radiusOut        = 200

#--- Dump file names

# TODO: get dump file names as input
dumpFileNames = ['data/disc_00000.ascii', 'data/disc_00006.ascii']

# ---------------------------------------------------------------------------- #

dumps = list()

for dumpFileName in dumpFileNames:

#--- Read dump file

    print('Reading in data from dumpfile: ' + dumpFileName + '... ', end='',
          flush=True)

    dump = PhantomDump(dumpFileName)

    print('done')

    dumps.append(dump)

#--- Perform analysis

radialAverages = list()
particleData = list()
sinkData = list()

for dump in dumps:

    print('Performing disc analysis')
    results = disc_analysis( radiusIn=radiusIn,
                             radiusOut=radiusOut,
                             numberRadialBins=numberRadialBins,
                             particleData=dump.ParticleData,
                             sinkData=dump.SinkData,
                             parameters=dump.Parameters,
                             units=dump.Units )

    radialAverages.append(results[0])
    particleData.append(results[1])
    sinkData.append(results[2])
