#!/usr/bin/python

import numpy as np
import sys

def main():
    filep = sys.argv[1]             # Path for the xyz structure where you want to add noise
    maxnoise = float(sys.argv[2])   # Maximum noise that can be added to an atom
    
    natom,atoms,x,y,z = read_xyz(filep)

    rx = x + random_array(natom,-maxnoise,maxnoise)
    ry = y + random_array(natom,-maxnoise,maxnoise)
    rz = z + random_array(natom,-maxnoise,maxnoise)

    print_xyz(atoms,rx,ry,rz)
    # print x,rx

    # import matplotlib.pyplot as plt
    # plt.plot(rx)
    # plt.plot(x)
    # plt.ylabel('some numbers')
    # plt.show()

def print_xyz(atoms,x,y,z):
    
    txt = '%5i\n' % len(atoms)
    txt += 'Randomized atoms by %s\n' % sys.argv[0]
    for i,at in enumerate(atoms):
        txt += '%3s %12.6f %12.6f %12.6f\n' % tuple([at,x[i],y[i],z[i]])
    
    print txt
    
def random_array(size,minr,maxr):
    return np.random.uniform(minr,maxr,size)


def read_xyz(filep):
    filec = open(filep,'r').readlines()
    natom = int(filec[0].strip())

    x = np.zeros(natom)
    y = np.zeros(natom)
    z = np.zeros(natom)
    atoms = []

    for i,line in enumerate(filec[2:]):
        atom,x[i],y[i],z[i] = line.split()
        atoms.append(atom)
        
    return natom,atoms,x,y,z

if __name__ == '__main__':
    main()
