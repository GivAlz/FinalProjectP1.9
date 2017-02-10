from ctypes import *
import os, sys

dso = CDLL('../ljmd.so')


def read_file(inputfile):
    f = open(inputfile,"r")
    lines = [line.rstrip('\n').split(" ")[0] for line in f]

    natoms,mass,epsilon,sigma,rcut,box,\
    restfile,trajfile,ergfile,nstep,dt,nprint = lines

    c_natoms = c_int(int(natoms))
    c_nstep = c_int(int(float(nstep)))
    c_nprint = c_int(int(nprint))


    c_dt = c_double(float(dt))
    c_mass = c_double(float(mass))
    c_epsilon = c_double(float(epsilon))
    c_sigma = c_double(float(sigma))
    c_rcut = c_double(float(rcut))
    c_box = c_double(float(box))

    c_restfile = c_char_p(restfile.encode())
    c_trajfile = c_char_p(trajfile.encode())
    c_ergfile = c_char_p(ergfile.encode())

    dso.compute_ljmd(c_natoms, c_mass, c_rcut, c_box, c_epsilon, c_sigma, c_restfile, c_trajfile,\
    c_ergfile, c_nstep, c_dt, c_nprint)




if __name__ == "__main__":

    if(len(sys.argv)>1):
        read_file(sys.argv[-1])
