#include "utilities.h"
/* helper function: zero out an array */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

#ifdef __MPI_H__
#include <stdlib.h> 
void set_mpi(mdsys_t* sys,temp_t* tmp,int rank,int nprocs)
{
	sys->nprocs = nprocs;
	sys->rank = rank;
	sys->slice = sys->natoms%nprocs;
	sys->start = sys->slice*sys->rank;
	/*Handling reminders:
	 * if the number of atoms can't be distributed exactly
	 * the first cores get one extra atom
	 * thus start is shifted*/
	if(sys->natoms%nprocs != 0){
	    if(sys->rank < sys->natoms%nprocs){
			++sys->slice;
			sys->start += sys->rank;
		}
		else{
			sys->start += sys->nprocs;
		}
	}
	tmp->tmpx = (double *)malloc(sys->slice*sizeof(double));
	tmp->tmpy = (double *)malloc(sys->slice*sizeof(double));
	tmp->tmpz = (double *)malloc(sys->slice*sizeof(double));
}
#endif
