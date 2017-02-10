#ifndef UTIL_H
#define UTIL_H

void azzero(double *d, const int n);

#ifdef __MPI_H__
#include "prototypes.h"
void set_mpi(mdsys_t* sys,temp_t* tmp,int rank,int nprocs);
#endif


#endif
