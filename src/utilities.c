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
	tmp->nprocs = nprocs;
	tmp->rank = rank;
	//tmp->slice = sys->natoms%nprocs;
	//tmp->start = tmp->slice*tmp->rank;
	/*Handling reminders:
	 * if the number of atoms can't be distributed exactly
	 * the first cores get one extra atom
	 * thus start is shifted*/
	//if(sys->natoms%nprocs != 0){
	    //if(tmp->rank < sys->natoms%nprocs){
			//++tmp->slice;
			//tmp->start += tmp->rank;
		//}
		//else{
			//tmp->start += tmp->nprocs;
		//}
	//}
	tmp->tmpx = (double *)malloc(sys->natoms*sizeof(double));
	tmp->tmpy = (double *)malloc(sys->natoms*sizeof(double));
	tmp->tmpz = (double *)malloc(sys->natoms*sizeof(double));
}
#endif
