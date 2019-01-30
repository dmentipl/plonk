'''
analyze_disc.py

Daniel Mentiplay, 2019.
'''

from plonk.analysis.disc import disc_analysis
from plonk.dumps import Dump

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

    dump = Dump()
    dump.read_dump(dumpFileName)

    print('done')

    dumps.append(dump)

#--- Perform analysis

radialBins      = list()
surfaceDensity  = list()
midplaneDensity = list()
smoothingLength = list()
scaleHeight     = list()
angularMomentum = list()
tilt            = list()
twist           = list()
psi             = list()
eccentricity    = list()
Stokes          = list()

for dump in dumps:

    print('Performing disc analysis on dumpfile: ' + dump.filename)

    analysis = disc_analysis(dump, radiusIn, radiusOut, numberRadialBins)

    radialBins.     append( analysis[0] )
    surfaceDensity. append( analysis[1] )
    midplaneDensity.append( analysis[2] )
    smoothingLength.append( analysis[3] )
    scaleHeight.    append( analysis[4] )
    angularMomentum.append( analysis[5] )
    tilt.           append( analysis[6] )
    twist.          append( analysis[7] )
    psi.            append( analysis[8] )
    eccentricity.   append( analysis[9] )
    Stokes.         append( analysis[10] )

output = \
'''
The following variables are available:

  - radialBins
  - surfaceDensity
  - midplaneDensity
  - smoothingLength
  - scaleHeight
  - angularMomentum
  - tilt
  - twist
  - psi
  - eccentricity
  - Stokes

for both gas and dust (if available).

There is also the list of dump objects 'dumps' each of which contain all
available information about the dump.
'''

print(output)
