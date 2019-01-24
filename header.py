'''
header.py

Daniel Mentiplay, 2019.
'''

class Header:
    '''
    Header class represents the header from a dump file (in the format produced
    by the phantom showheader utility program).
    '''

    def __init__(self, filename):

        self.filename = filename

        self._read_file()

    def _read_file(self):

        with open(self.filename, 'r') as file:

            keys = list()
            values = list()

            firstLine = True
            for line in file:
                if firstLine:
                    firstLine = False
                    continue
                keys.append(line.rstrip('\n').split()[0])
                value = line.rstrip('\n').split()[1]
                if '.' in value:
                    value = float(value)
                else:
                    value = int(value)
                values.append(value)

        newKeys = list()
        newValues = list()

        prevKey = None
        firstMultipleKey = True
        idxj = 0

        for idxi, key in enumerate(keys):

            if key == prevKey:

                if firstMultipleKey:
                    newValues[idxj-1] = list()
                    newValues[idxj-1].append(values[idxi-1])
                    newValues[idxj-1].append(values[idxi])
                    firstMultipleKey = False
                else:
                    newValues[idxj-1].append(values[idxi])

            else:

                firstMultipleKey = True
                idxj += 1
                newKeys.append(key)
                newValues.append(values[idxi])

            prevKey = key

        self.header = dict(zip(newKeys, newValues))
