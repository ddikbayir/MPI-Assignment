/* 
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 */

#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include "time.h"
#include "apf.h"
#include "Plotting.h"
#ifdef _MPI_
#include <mpi.h>
#endif
#define NEIGHBOR 1
#define FINAL 2
using namespace std;

void repNorms(ofstream& logfile, double l2norm, double mx, double dt, int m,int n, int niter, int stats_freq);


// Reports statistics about the computation: the L2 Norm and the Infinity NOrm
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem


// The L2 norm of an array is computed by taking sum of the squares
// of each element, normalizing by dividing by the number of points
// and then taking the sequare root of the result
//
// The Linf norm is simply the maximum (absolute) value over
// all points in the array

 double stats(double **E, int m, int n, double *_mx){
     double mx = -1;
     double l2norm = 0;
     int i, j;
     for (j=1; j<=m+1; j++)
       for (i=1; i<=n+1; i++) {
	   l2norm += E[j][i]*E[j][i];
	   double fe = fabs(E[j][i]);
	   if (fe > mx)
	       mx = fe;
      }

// In the parallel version, you must sum all the local contributoins
// before dividing by (m+1)*(n+1)
     l2norm /= (double) ((m+1)*(n+1));
     l2norm = sqrt(l2norm);

     *_mx = mx;
     return l2norm;
 }
// Added px and py to solve
int solve(ofstream& logfile, double ***_E, double ***_E_prev, double **R, int m, int n, int niters, double alpha, double dt, int plot_freq, Plotter *plotter, int stats_freq, int px, int py){

 // Simulated time is different from the integer timestep number
 double t = 0.0;
 int siz = m*n*sizeof(double);
 char * buffer = new char[siz];
 int position = 0;
 double **E = *_E, **E_prev = *_E_prev;
 int niter;
 int rank =0, np=1;
 #ifdef _MPI_
 MPI_Comm_size(MPI_COMM_WORLD,&np);
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Request send, recv, send2, recv2;
 #endif
 cout << "Inside solve from rank " << rank << "/" << np << endl;
 // We continue to sweep over the mesh until the simulation has reached
 // the desired simulation Time
 // This is different from the number of iterations

/* Setting up our boundaries for i, not pretty but it works */
int i_s = (((n+1)*rank)/py)+1;
i_s += rank==0?0:1;
int i_e = rank==py-1?n+1:(((n+1)*(rank+1))/py)+1;

//   cout << "i_s = " << i_s << " rank = " << rank << " n = " << n+1 << endl;
//cout << "i_e = " << i_e << " rank = " << rank << " n = " << n+1 << endl;
  for (niter = 0; niter < niters; niter++){
  
#ifdef DEBUG
   double mx;
   double l2norm = stats(E_prev,m,n,&mx);
   repNorms(logfile,l2norm,mx,dt,m,n,niter, stats_freq);
   if (plot_freq)
	plotter->updatePlot(E,  niter, m+1, n+1, WAIT);
//    splot(E_prev,niter,m+1,n+1,WAIT);
#endif

   /* 
    * Copy data from boundary of the computational box to the
    * padding region, set up for differencing computational box's boundary
    *
    * These are physical boundary conditions, and are not to be confused
    * with ghost cells that we would use in an MPI implementation
    *
    * The reason why we copy boundary conditions is to avoid
    * computing single sided differences at the boundaries
    * which increase the running time of solve()
    *
    */
    
   int i,j;
/*
     for (j=1; j<=m+1; j++) {
       E_prev[j][0] = E_prev[j][2];
	E_prev[j][n+2] = E_prev[j][n];
     } 
   for (i=i_s; i<=i_e; i++) {
      E_prev[0][i] = E_prev[2][i];
      E_prev[m+2][i] = E_prev[m][i];
   }*/
   
   // Solve for the excitation, a PDE 

   
   for (j=1; j<=m+1; j++){
     E_prev[j][0] = E_prev[j][2];
     E_prev[j][n+2] = E_prev[j][n];
     for (i=i_s; i<=i_e; i++) {
	E_prev[0][i] = E_prev[2][i];
        E_prev[m+2][i] = E_prev[m][i];
	E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
     }
	/* 1 | 2 | 3	0 | 1 | 2 r/3= 0
	 *-----------	---------
	 * 4 | 5 | 6	3 | 4 | 5 = 1
	 *-----------	---------
	 * 7 | 8 | 9	6 | 7 | 8 = 2
	 * 1  >1   0 */
	// Left edge
	if ((rank+1)%py== 1) {
		/* This sends the boundary of the first thread (the items at the end of first thread) */
		MPI_Send(&E[j][i_e], 1, MPI_DOUBLE, rank+1, NEIGHBOR, MPI_COMM_WORLD);
		/* This recieves the boundary of the second thread (i+1) */
		MPI_Recv(&E[j][i_e+1], 1, MPI_DOUBLE, rank+1, NEIGHBOR, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
	// Middle
	} else if ((rank+1)%py>1) {
		/* This sends the boundary of the second thread (the items at the end of the second thread) */
		MPI_Send(&E[j][i_s], 1, MPI_DOUBLE, rank-1, NEIGHBOR, MPI_COMM_WORLD);
		/* This recieves the boundary of the first thread needed for computation (i-1) */
		MPI_Recv(&E[j][i_s-1], 1, MPI_DOUBLE, rank-1, NEIGHBOR, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
		MPI_Send(&E[j][i_e], 1, MPI_DOUBLE, rank+1, NEIGHBOR, MPI_COMM_WORLD);
		/* This recieves the boundary of the second thread (i+1) */
		MPI_Recv(&E[j][i_e+1], 1, MPI_DOUBLE, rank+1, NEIGHBOR, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

	// Right Edge
	} else if ((rank+1)%py==0) {
		/* This sends the boundary of the second thread (the items at the end of the second thread) */
		MPI_Send(&E[j][i_s], 1, MPI_DOUBLE, rank-1, NEIGHBOR, MPI_COMM_WORLD);
		/* This recieves the boundary of the first thread needed for computation (i-1) */
		MPI_Recv(&E[j][i_s-1], 1, MPI_DOUBLE, rank-1, NEIGHBOR, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
	}	
    }

/* Printing our table each iteration. Just for debugging 
cout << "------------------------n = " << niter << ", rank = " << rank << "-------------------------\n";
for (j=1; j<=m+1; j++){
 for (i=1; i<=n+1; i++) {
   cout << E[j][i] << " ";
  }
  cout << "\n";
}*/

   /* 
    * Solve the ODE, advancing excitation and recovery variables
    *     to the next timtestep
    */
   int new_i_s = i_s==1?i_s:i_s-1;
   int new_i_e = i_e==n+1?i_e:i_e+1;
   for (j=1; j<=m+1; j++){
     double *RR = &R[j][new_i_s];
     double *EE = &E[j][new_i_s];
     for (i=new_i_s; i<=new_i_e; i++, EE++,RR++) {
	EE[0] += -dt*(kk*EE[0]*(EE[0]-a)*(EE[0]-1)+EE[0]*RR[0]);
	RR[0] += dt*(epsilon+M1* RR[0]/( EE[0]+M2))*(-RR[0]-kk*EE[0]*(EE[0]-b-1));
     }
   }


   if (stats_freq){
	   double mx;
     double l2norm = stats(E_prev,m,n,&mx);
     repNorms(logfile,l2norm,mx,dt,m,n,niter, stats_freq);
  
   }

   if (plot_freq){
          if (!(niter % plot_freq)){
//            splot(E,niter,m+1,n+1,WAIT);
/* *************************** */
/* Gather matrix for plot      */
/* *************************** */

/* We're going to let thread 0 do the return value, so that means every other thread
 * needs to pack up it's own working segment and send it to thread 0 */
if (rank != 0) {
	position = 0;
	/* Packing up our rank so we know which part to fill in */
	MPI_Pack(&rank, 1, MPI_INT, buffer, siz, &position, MPI_COMM_WORLD);

/* Packing up our array */
for (int j=1; j<=m+1; j++){
  for (int i=i_s; i<=i_e; i++) {
    MPI_Pack(&E_prev[j][i],1, MPI_DOUBLE, buffer, siz, &position, MPI_COMM_WORLD);
  }
}
	/* Then we send to thread 0 */
	MPI_Send(buffer, siz, MPI_PACKED,0,FINAL, MPI_COMM_WORLD);
}else { 
// IN THREAD 0

int sofar = 0;
int in_rank;
int in_i_s;
int in_i_e;
double in_val;
	/* We loop until we've recieved communication from every thread */
	while (sofar < np-1) {
		position=0; // Reset buffer iterator each time we recieve a message
		/* We recieve our packed stuff from some other thread */
		MPI_Recv(buffer, siz, MPI_PACKED, MPI_ANY_SOURCE,FINAL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		/* We unpack the rank first */
		MPI_Unpack(buffer, siz, &position, &in_rank, 1, MPI_INT, MPI_COMM_WORLD);
		/* Use the rank to figure out the thread's boundary so we can fill our own array
		 * with it's array elements */
		in_i_s = (((n+1)*in_rank)/py)+1;
		in_i_s += in_rank==0?0:1;
		in_i_e = in_rank==py-1?n+1:(((n+1)*(in_rank+1))/py)+1;

		/* Now we unpack each value in the array and add it to our own */
		for (int j=1; j<=m+1; j++){
  			for (int i=in_i_s; i<=in_i_e; i++) {
				MPI_Unpack(buffer, siz, &position, &in_val, 1, MPI_DOUBLE, MPI_COMM_WORLD);
					E_prev[j][i] = in_val;
			}
		}
		sofar++; // Move to the next process
	}
		    plotter->updatePlot(E,  niter, m+1, n+1, WAIT);
}

        }
    }

   // Swap current and previous
   double **tmp = E; E = E_prev; E_prev = tmp;
 }

  // Store them into the pointers passed in

/* *************************** */
/* Gather matrix before return */
/* *************************** */

/* We're going to let thread 0 do the return value, so that means every other thread
 * needs to pack up it's own working segment and send it to thread 0 */
if (rank != 0) {
	position = 0;
	/* Packing up our rank so we know which part to fill in */
	MPI_Pack(&rank, 1, MPI_INT, buffer, siz, &position, MPI_COMM_WORLD);

/* Packing up our array */
for (int j=1; j<=m+1; j++){
  for (int i=i_s; i<=i_e; i++) {
    MPI_Pack(&E_prev[j][i],1, MPI_DOUBLE, buffer, siz, &position, MPI_COMM_WORLD);
  }
}
	/* Then we send to thread 0 */
	MPI_Send(buffer, siz, MPI_PACKED,0,FINAL, MPI_COMM_WORLD);
}else { 
// IN THREAD 0

int sofar = 0;
int in_rank;
int in_i_s;
int in_i_e;
double in_val;
	/* We loop until we've recieved communication from every thread */
	while (sofar < np-1) {
		position=0; // Reset buffer iterator each time we recieve a message
		/* We recieve our packed stuff from some other thread */
		MPI_Recv(buffer, siz, MPI_PACKED, MPI_ANY_SOURCE,FINAL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		/* We unpack the rank first */
		MPI_Unpack(buffer, siz, &position, &in_rank, 1, MPI_INT, MPI_COMM_WORLD);
		/* Use the rank to figure out the thread's boundary so we can fill our own array
		 * with it's array elements */
		in_i_s = (((n+1)*in_rank)/py)+1;
		in_i_s += in_rank==0?0:1;
		in_i_e = in_rank==py-1?n+1:(((n+1)*(in_rank+1))/py)+1;

		/* Now we unpack each value in the array and add it to our own */
		for (int j=1; j<=m+1; j++){
  			for (int i=in_i_s; i<=in_i_e; i++) {
				MPI_Unpack(buffer, siz, &position, &in_val, 1, MPI_DOUBLE, MPI_COMM_WORLD);
					E_prev[j][i] = in_val;
			}
		}
		sofar++; // Move to the next process
	}

	/* Print out final array for debugging 
	cout << "------------------------FINAL-------------------------\n";
	for (int j=1; j<=m+1; j++){
	 for (int i=1; i<=n+1; i++) {
	   cout << E_prev[j][i] << " ";
  	}
  	cout << "\n";
	}*/
}

  *_E = E;
  *_E_prev = E_prev;
  return niter;

}
