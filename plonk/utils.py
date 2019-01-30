'''
utils.py

Daniel Mentiplay, 2019.
'''

def print_warning(message):
    '''
    Print formatted warning message.
    '''
    message = 'WARNING: ' + message
    print( '\n' \
         + len(message)*'-' + '\n' \
         + message + '\n' \
         + len(message)*'-' )

def print_error(message):
    '''
    Print formatted error message.
    '''
    message = 'ERROR: ' + message
    print( '\n' \
         + len(message)*'-' + '\n' \
         + message + '\n' \
         + len(message)*'-' )
