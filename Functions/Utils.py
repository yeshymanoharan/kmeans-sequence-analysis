import os, os.path
import errno
from sys import argv

def mkdir_p(path):
    try:
        #print(os.path)
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def safe_open(path, open_type):
    ''' Open "path" for writing, creating any parent directories as needed.
    '''
    if open_type == "w":
        mkdir_p(os.path.dirname(path))
    return open(path, open_type)

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts
