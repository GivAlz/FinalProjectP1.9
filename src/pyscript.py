from ctypes import *
import os, sys

dso = CDLL('./ljmd_serial.so')


def read_file(inputfile):
    f = open(inputfile,"r")
    lines = [line.rstrip('\n').split(" ")[0] for line in f]

    natoms,mass,epsilon,sigma,rcut,box,\
    restfile,trajfile,ergfile,line,nstep,dt = lines

    c_natoms = c_int(int(natoms))
    c_sigma = c_double(float(sigma))
    c_restfile = c_char_p(restfile)

    dso.mklkop(c_)

    return c_natoms, c_sigma, c_restfile.value


if __name__ == "__main__":

    if(len(sys.argv)>1):
        print read_file(sys.argv[-1])
