import eosWrap as eos
import numpy as np
from Wrapper import Wrapper
import matplotlib.pyplot as plt
import Models
from tabulate import tabulate
from os.path import join
import sys
import getopt

def main(argv):
    HYPER = 0
    try:
        opts, args = getopt.getopt(argv, 'hv', ['help', 'verbose'])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('Help me!')
            sys.exit()
        if opt == '--hyper':
            HYPER = 1
    
    
if __name__ == '__main__':
    main(sys.argv[1:])
    
        

