'''
visualize_disc.py

Daniel Mentiplay, 2019.
'''

from plonk.plonk.dump import PhantomDump
from plonk.plonk.visualization.image import Image

#--- Dump file names

# TODO: get dump file names as input
dump_file_names = ['data/disc_00000.ascii', 'data/disc_00006.ascii']

dumps = list()

for dump_file_name in dump_file_names:

#--- Read dump file

    print('Reading in data from dumpfile: ' + dump_file_name + '... ', end='',
          flush=True)

    dump = PhantomDump(dump_file_name)

    print('done')

    dumps.append(dump)

#--- Plot image

print('Creating Image objects for each dump file...', end='')

images = list()
for dump in dumps:
    images.append(Image(dump))

print('done')
print('\n\nTo plot the ith dump:  images[i].plot()')
