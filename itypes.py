'''
itypes.py

Daniel Mentiplay, 2019.
'''

class ITypes:
    '''
    Dump class represents a dump file (ascii only).
    '''

    # pylint: disable=too-many-instance-attributes

    def __init__(self):

        self.iGas  = 1
        self.iSink = 3
        self.iDust = 8

        self.iGasLabel = 'gas'
        self.iSinkLabel = 'sink'
        self.iDustLabel = 'dust'

        self.maxDustTypes = 10

# Instantiate Types class to create global types object.
iTypes = ITypes()
