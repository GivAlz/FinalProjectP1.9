#ifndef PROT_H
#define PROT_H

/* generic file- or pathname buffer length */
#define BLEN 200
/* a few physical constants */
extern const double kboltz;     /* boltzman constant in kcal/mol/K */
extern const double mvsq2e; /* m*v^2 in kcal/mol */
/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    
	#ifdef __MPI_H__
	//These are values needed to handle parallelization
	int rank, nprocs, slice,start;
	#endif
};
typedef struct _mdsys mdsys_t;

#ifdef __MPI_H__
//These temporary values are used to simplify the parallelization

struct _temp {
	double *tmpx,*tmpy,*tmpz;
};
typedef struct _temp temp_t;


#endif

#endif
