'''
do_analysis_disc.py

Daniel Mentiplay, 2019.
'''

from analysis_disc import disc_analysis
from dump import Dump

#--- Dump file name

dumpFilePrefix = 'disc_00000'  # TODO: get dump filename as input
                               # TODO: read multiple dumpfiles

# ---------------------------------------------------------------------------- #


#--- Read dump file

print('\nReading in data from dumpfile: ' + dumpFilePrefix + '\n')

dump = Dump(dumpFilePrefix)

#--- Perform analysis

print('Performing disc analysis on dumpfile: ' + dumpFilePrefix + '\n')

analysis = disc_analysis(dump)

radialBins      = analysis[0]
surfaceDensity  = analysis[1]
midplaneDensity = analysis[2]
smoothingLength = analysis[3]
scaleHeight     = analysis[4]
angularMomentum = analysis[5]
tilt            = analysis[6]
twist           = analysis[7]
psi             = analysis[8]
eccentricity    = analysis[9]
Stokes          = analysis[10]

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
'''

print(output)
