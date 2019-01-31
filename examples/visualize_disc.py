'''
visualize_disc.py

Daniel Mentiplay, 2019.
'''

from plonk.dumps import Dump
from plonk.visualization.image import Image

#--- Dump file

filename = 'data/disc_00000.ascii'

#--- Read dump file

print('Reading dumpfile...')

dump = Dump()
dump.read_dump(filename)

#--- Plot image

print('Plotting...')
image = Image(dump)
image.plot()
