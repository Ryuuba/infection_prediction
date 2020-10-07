#!/usr/bin/env python3

import sys
import numpy as np

def main():
    rij_file = sys.argv[1]
    f_a = open('aij.txt', 'w')
    with open(rij_file) as f:
        for line in f.readlines():
            for token in line.split():
                print(token)
                val = float(token) 
                if val > 0.0:
                    f_a.write('1 ')
                else:
                    f_a.write('0 ')
            f_a.write('\n')
    f_a.close()
                

if __name__ == "__main__":
    main()