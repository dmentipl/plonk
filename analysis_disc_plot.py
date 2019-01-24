'''
analysis_disc_plot.py

Daniel Mentiplay, 2019.
'''

from analysis_disc import disc_analysis
from dump import Dump

#--- Dump file name

dumpFilePrefix = 'disc_00000'  # TODO: get dump filename as input
                               # TODO: read multiple dumpfiles

# ---------------------------------------------------------------------------- #


#--- Read dump file

dumpObj = Dump(dumpFilePrefix)

#--- Perform analysis

radialBinsDisc, \
surfaceDensityGas, \
midplaneDensityGas, \
meanSmoothingLengthGas, \
scaleHeightGas, \
meanAngularMomentumGas, \
tiltGas, \
twistGas, \
psiGas, \
meanEccentricityGas, \
surfaceDensityDust, \
midplaneDensityDust, \
meanSmoothingLengthDust, \
scaleHeightDust, \
meanAngularMomentumDust, \
tiltDust, \
twistDust, \
psiDust, \
eccentricityDust, \
meanEccentricityDust = disc_analysis(dumpObj)
