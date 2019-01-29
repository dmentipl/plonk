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
        self.iDust = 7

        self.iGasSplash  = 1
        self.iSinkSplash = 3
        self.iDustSplash = 8

        self.iGasLabel  = 'gas'
        self.iSinkLabel = 'sink'
        self.iDustLabel = 'dust'

# Instantiate Types class to create global types object.
iTypes = ITypes()
