/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#include "prototypes.h"
#include "input.h"
#include "utilities.h"
#include "force.h"
#include "output.h"
#include "verlet1.h"

#ifdef __MPI_H__
#include <mpi.h>
#endif

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */


/* main */
int main(int argc, char **argv)
{
    int nprint;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *traj,*erg;
    mdsys_t sys;

	#ifndef __MPI_H__
	read_input(&sys, &nprint,restfile,trajfile,ergfile,line);
	#endif

    
	#ifdef __MPI_H__
	int rank, nprocs;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
	
		
	//This function sets the values needed for splitting the array
	temp_t tmp;
	set_mpi(&sys,&tmp,rank,nprocs);
	
	if(rank==0)
	{
		read_input(&sys, &nprint,restfile,trajfile,ergfile,line);
	}
	//Using MPI to Broadcast a structure is a nightmare
	//This is easier:
	MPI_Bcast(&sys.natoms,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&sys.mass,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sys.epsilon,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sys.sigma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sys.rcut,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sys.box,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.nsteps,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));
	
    azzero(sys.fx, sys.natoms);
    azzero(sys.fy, sys.natoms);
    azzero(sys.fz, sys.natoms);

    MPI_Bcast(sys.rx,sys.natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(sys.ry,sys.natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(sys.rz,sys.natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);

    MPI_Bcast(sys.vx,sys.natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(sys.vy,sys.natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(sys.vz,sys.natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);

	#endif
	
	
    /* initialize forces and energies.*/
    sys.nfi=0;
    #ifndef __MPI_H__
    force(&sys);
    #endif
    
    #ifdef __MPI_H__
    force(&sys,&tmp);
    #endif
    ekin(&sys);

    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
		#ifndef __MPI_H__
        velverlet(&sys);
        #endif
		#ifdef __MPI_H__
		velverlet(&sys,&tmp);
		#endif
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done.\n");
    fclose(erg);
    fclose(traj);

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
	

	#ifdef __MPI_H__
	free(tmp.tmpx);
	free(tmp.tmpy);
	free(tmp.tmpz);
	MPI_Finalize();
	#endif
    return 0;
}
