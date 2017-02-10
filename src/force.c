#include "force.h"

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#ifdef __MPI_H__
#include <mpi.h>
#include <string.h>
#endif

/* helper function: apply minimum image convention */
double pbc(double x, const double boxby2,const double twice_boxby2)
{
    while (x >  boxby2) x -= twice_boxby2;
    while (x < -boxby2) x += twice_boxby2;
    return x;
}


/* compute kinetic energy */
void ekin(mdsys_t *sys)
{   
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}
#ifndef __MPI_H__
/* compute forces */
void force(mdsys_t *sys) 
{
    double ffac;
    double rx,ry,rz;
    int i,j;
    const double half_box = 0.5*sys->box;
    const double twice_boxby2 = sys->box;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
	
	//Constants added to avoid computing pow,sqrt etc.
	
	double c12 = 4.0*sys->epsilon*pow(sys->sigma,12.0);
	double c6 = 4.0*sys->epsilon*pow(sys->sigma,6.0);
	const double rcsq = sys->rcut * sys->rcut;
	double rsq,rinv,r6;

	#ifdef _OPENMP
	double epot=0.0;
	#pragma omp parallel for private(i,j,rx,ry,rz,rsq,rinv,r6,ffac) shared(sys) reduction(+:epot)
	#endif

    for(i=0; i < (sys->natoms)-1; ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            //if (i==j) continue; //With 3rd law we can eliminate this
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], half_box,twice_boxby2);
            ry=pbc(sys->ry[i] - sys->ry[j], half_box,twice_boxby2);
            rz=pbc(sys->rz[i] - sys->rz[j], half_box,twice_boxby2);
            //r = sqrt(rx*rx + ry*ry + rz*rz);
      
            /* compute force and energy if within cutoff */
            /**if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);
                
                //We use Netwon third law: the following is used twice
                //So it's not multiplied by 0.5 anymore
                sys->epot += 4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));
			**/
			rsq = rx*rx + ry*ry + rz*rz;
			if (rsq < rcsq){
				
				rinv=1.0/rsq;
				r6 = rinv*rinv*rinv;
				
				ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
				#ifndef _OPENMP
				sys->epot += r6*(c12*r6 - c6);
				#endif
				
				#ifdef _OPENMP
				epot += r6*(c12*r6 - c6);
				#endif
                sys->fx[i] += rx*ffac;
                sys->fx[j] -= rx*ffac;
                
                sys->fy[i] += ry*ffac;
                sys->fy[j] -= ry*ffac;
                
                sys->fz[i] += rz*ffac;
                sys->fz[j] -= rz*ffac;
            }
        }
    }
    #ifdef _OPENMP
	sys->epot = epot;
	#endif
}
#endif

#ifdef __MPI_H__
void force(mdsys_t *sys, temp_t* tmp)
{
    double ffac,epot;
    double rx,ry,rz;
    int i,j,ii;
    const double half_box = 0.5*sys->box;
    const double twice_boxby2 = sys->box;

	
    /* zero energy and forces */
    epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    //azzero(tmp->tmpx,sys->natoms);
    //azzero(tmp->tmpy,sys->natoms);
    //azzero(tmp->tmpz,sys->natoms);
	
    MPI_Bcast(sys->fx,sys->natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(sys->fy,sys->natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(sys->fz,sys->natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	memcpy(tmp->tmpx,sys->fx,sys->natoms*sizeof(double));
	memcpy(tmp->tmpy,sys->fy,sys->natoms*sizeof(double));
	memcpy(tmp->tmpz,sys->fz,sys->natoms*sizeof(double));
	//Constants added to avoid computing pow,sqrt etc.
	
	double c12 = 4.0*sys->epsilon*pow(sys->sigma,12.0);
	double c6 = 4.0*sys->epsilon*pow(sys->sigma,6.0);
	const double rcsq = sys->rcut * sys->rcut;
	double rsq,rinv,r6;
	
    for(i=0; i < (sys->natoms)-1; i+=tmp->nprocs) {
		ii = i + tmp->rank;
		if(ii>=(sys->natoms-1)) break;
        for(j=i+1; j < (sys->natoms); ++j) {
            /* particles have no interactions with themselves */
            //if (i==j) continue; //With 3rd law we can eliminate this
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], half_box,twice_boxby2);
            ry=pbc(sys->ry[i] - sys->ry[j], half_box,twice_boxby2);
            rz=pbc(sys->rz[i] - sys->rz[j], half_box,twice_boxby2);
            //r = sqrt(rx*rx + ry*ry + rz*rz);
      
            /* compute force and energy if within cutoff */
            /**if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);
                
                //We use Netwon third law: the following is used twice
                //So it's not multiplied by 0.5 anymore
                sys->epot += 4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));
			**/
			rsq = rx*rx + ry*ry + rz*rz;
			if (rsq < rcsq){
				rinv=1.0/rsq;
				r6 = rinv*rinv*rinv;
				
				ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
				epot += r6*(c12*r6 - c6);
				
                tmp->tmpx[i] += rx*ffac;
                tmp->tmpx[j] -= rx*ffac;
                
                tmp->tmpy[i] += ry*ffac;
                tmp->tmpy[j] -= ry*ffac;
                
                tmp->tmpz[i] += rz*ffac;
                tmp->tmpz[j] -= rz*ffac;
            }
        }
    }
    
    MPI_Reduce(tmp->tmpx,sys->fx,sys->natoms,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(tmp->tmpy,sys->fy,sys->natoms,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(tmp->tmpz,sys->fz,sys->natoms,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&epot,&sys->epot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
}
#endif
