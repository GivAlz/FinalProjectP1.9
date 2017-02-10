#ifndef LJMD_H
#define LJMD_H

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#include "prototypes.h"
#include "utilities.h"
#include "force.h"
#include "output.h"
#include "verlet1.h"



int read_from_py(int atoms, double mass, double rcut, double box, double epsilon,
   double sigma, char * restfile, char * trajfile, char * ergfile,
   int nsteps, double dt, mdsys_t * sys);

int compute_ljmd(int atoms, double mass, double rcut,  double box, double epsilon,
      double sigma, char restfile[BLEN], char trajfile[BLEN], char ergfile[BLEN], int nsteps, double dt, int nprint);

#endif
