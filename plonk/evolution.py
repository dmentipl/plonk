'''
evolution.py

Daniel Mentiplay, 2019.
'''

import numpy as np

class Evolution:
    '''
    Evolution class represents the data from Phantom .ev files.

    Arguments:
        filename : e.g. disc01.ev, discSink0001N01.ev
    '''

    def __init__(self, filename):

        self.filename = filename

        self.data = np.loadtxt(filename)

        with open(filename) as file:
            firstLine = file.readline().rstrip('\n').split('[')

        labels = list()
        for col in firstLine[1:]:
            labels.append(col[2:].rstrip().rstrip(']').lstrip())

        self.labels = labels
