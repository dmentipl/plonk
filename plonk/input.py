'''
pysplash.py

Daniel Mentiplay, 2019.
'''

import matplotlib.pyplot as plt
import numpy as np

import splash

from .dump import Dump
from .utils import density_from_smoothing_length

#--- Options

dumpFilePrefix = 'disc_00000'

npixx = 1000
npixy = 1000

index = 9

#--- Splash subroutines

interpolate2d = splash.splash.interpolate2d
set_interpolation_weights = splash.splash.set_interpolation_weights

#--- Read dump file

print('Reading dumpfile...')

dump = Dump(dumpFilePrefix)

npart = dump.parameters['npartoftype'][0]
mass = dump.parameters['massoftype'][0]

dat = np.empty((npart,11), dtype=np.float32, order='F')

dat[:,  0] = dump.arrays.position['gas'][:, 0]
dat[:,  1] = dump.arrays.position['gas'][:, 1]
dat[:,  2] = dump.arrays.position['gas'][:, 2]
dat[:,  3] = np.array(npart*[dump.parameters['massoftype'][0]])
dat[:,  4] = dump.arrays.smoothingLength['gas']
dat[:,  5] = density_from_smoothing_length(dat[:, 4], mass)
dat[:,  6] = dump.arrays.velocity['gas'][:, 0]
dat[:,  7] = dump.arrays.velocity['gas'][:, 1]
dat[:,  8] = dump.arrays.velocity['gas'][:, 2]
dat[:,  9] = dump.arrays.dustFrac[:, 0]
dat[:, 10] = dump.arrays.dustFrac[:, 1]

itype = np.ones(npart)

ih = 1
irho = 1
ipmass = 1
ndim = 3
iamtypei = np.ones(1)
inormalise = False
usetype = True
ntypes = 1
masstype = mass
npartoftype = npart
ninterp = npart
ndataplots = 1
irescale = False
idensityweighted = False
inormalise = False
units = np.ones(1)
unit_interp = 1.
required = np.array(1*[False])
rendersinks = False

weight = np.empty(npart, dtype=np.float32, order='F')

#--- Interpolation weights

print('Set interpolation weights...')
set_interpolation_weights(weighti=weight,
                          dati=dat,
                          iamtypei=iamtypei,
                          usetype=usetype,
                          ninterp=ninterp,
                          npartoftype=npartoftype,
                          masstype=masstype,
                          ntypes=ntypes,
                          ndataplots=ndataplots,
                          irho=irho,
                          ipmass=ipmass,
                          ih=ih,
                          ndim=ndim,
                          irescale=irescale,
                          idensityweighted=idensityweighted,
                          inormalise=inormalise,
                          units=units,
                          unit_interp=unit_interp,
                          required=required,
                          rendersinks=rendersinks)

#--- Interpolate onto 2d pixel array

maxx      = np.max(dat[:, 0])
maxy      = np.max(dat[:, 1])
xmin      = - maxx
ymin      = - maxy
pixwidthx = 2*maxx / npixx
pixwidthy = 2*maxy / npixy

normalise = False
exact     = False
periodicx = False
periodicy = False

datsmooth = np.zeros((npixx,npixy), dtype=np.float32, order='F')

print('Interpolate to 2d pixel array...')
interpolate2d(x=dat[:, 0],
              y=dat[:, 1],
              hh=dat[:, 4],
              weight=weight,
              dat=dat[:, index],
              itype=itype,
              npart=npart,
              xmin=xmin,
              ymin=ymin,
              datsmooth=datsmooth,
              npixx=npixx,
              npixy=npixy,
              pixwidthx=pixwidthx,
              pixwidthy=pixwidthy,
              normalise=normalise,
              exact=exact,
              periodicx=periodicx,
              periodicy=periodicy)

print('Plotting...')
plt.imshow(datsmooth)
