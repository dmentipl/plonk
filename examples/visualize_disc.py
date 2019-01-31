'''
visualize_disc.py

Daniel Mentiplay, 2019.
'''

from plonk.dumps import Dump
from plonk.visualization.image import Image

#--- Dump file names

# TODO: get dump file names as input
dumpFileNames = ['data/disc_00000.ascii', 'data/disc_00006.ascii']

dumps = list()

for dumpFileName in dumpFileNames:

#--- Read dump file

    print('Reading in data from dumpfile: ' + dumpFileName + '... ', end='',
          flush=True)

    dump = Dump()
    dump.read_dump(dumpFileName)

    print('done')

    dumps.append(dump)

#--- Plot image

print('Creating Image objects for each dump file...', end='')

images = list()
for dump in dumps:
    images.append(Image(dump))

print('done')
print('\n\nTo plot the ith dump:  images[i].plot()')
