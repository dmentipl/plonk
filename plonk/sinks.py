'''
sinks.py

Daniel Mentiplay, 2019.
'''

from evolution import Evolution

prefix = 'disc'

nSinks = 4

sink = list()
for i in range(1,nSinks+1):
    filename = prefix + 'Sink000' + str(i) + 'N01.ev'
    sink.append(Evolution(filename))
