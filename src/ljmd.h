#ifndef LJMD_H
#define LJMD_H

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#include "prototypes.h"
#include "utilities.h"



int read_from_py(int atoms, double mass, double epsilon,
   double sigma, char * restfile, char * trajfile, char * ergfile,
   char * line, int nsteps, int dt);

int compute_ljmd(int atoms, double mass, double epsilon,
      double sigma, char * restfile, char * trajfile, char * ergfile,
      char * line, int nsteps, int dt);

#endif
