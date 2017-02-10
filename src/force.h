#ifndef FORCE_H
#define FORCE_H

#include "prototypes.h"
#include "utilities.h"

double pbc(double x, const double boxby2, const double twice_boxby2);

void ekin(mdsys_t *sys);
#ifndef __MPI_H__

void force(mdsys_t *sys);

#endif

#ifdef __MPI_H__
void force(mdsys_t *sys, temp_t* tmp);
#endif

#endif

