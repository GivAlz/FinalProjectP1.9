#include "utilities.h"
/* helper function: zero out an array */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

void set_mpi(mdsys_t* sys,int* rank,int* nprocs)
{
	sys->nprocs = nprocs;
	sys->rank = rank;
}
