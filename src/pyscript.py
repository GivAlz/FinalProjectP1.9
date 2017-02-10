from ctypes import *
import os, sys

dso = CDLL('../ljmd_serial.so')


def read_file(inputfile):
    f = open(inputfile,"r")
    lines = [line.rstrip('\n').split(" ")[0] for line in f]

    natoms,mass,epsilon,sigma,rcut,box,\
    restfile,trajfile,ergfile,line,nstep,dt = lines

    c_natoms = c_int(int(natoms))
    c_nstep = c_int(int(float(nstep)))
    c_dt = c_int(int(dt))
    #c_line = c_int(int(line))

    c_mass = c_double(float(mass))
    c_epsilon = c_double(float(epsilon))
    c_sigma = c_double(float(sigma))
    c_rcut = c_double(float(rcut))
    c_box = c_double(float(box))

    c_restfile = c_char_p(restfile)
    c_trajfile = c_char_p(trajfile)
    c_ergfile = c_char_p(ergfile)

    #dso.read_from_py.arg = [c_char_p]
    #dso.read_from_py(c_natoms, c)

    #dso.compute_ljmd(c_natoms, c_mass, c_epsilon, c_sigma, c_restfile, c_trajfile,\
    #c_ergfile, c_nstep, c_dt)

    out = [c_natoms, c_mass, c_epsilon, c_sigma, c_restfile, c_trajfile,\
    c_ergfile, c_nstep, c_dt]

    return out


if __name__ == "__main__":

    if(len(sys.argv)>1):
        print read_file(sys.argv[-1])
        read_file(sys.argv[-1])
        #return
