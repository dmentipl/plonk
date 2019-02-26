'''
utils.py

Daniel Mentiplay, 2019.
'''

def print_warning(message):
    '''
    Print formatted warning message.
    '''
    l = len(message)
    l = min(l, 80)
    message = 'WARNING: ' + message
    print( '\n' \
         + l*'-' + '\n' \
         + message + '\n' \
         + l*'-' )

def print_error(message):
    '''
    Print formatted error message.
    '''
    l = len(message)
    l = min(l, 80)
    message = 'ERROR: ' + message
    print( '\n' \
         + l*'-' + '\n' \
         + message + '\n' \
         + l*'-' )
