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
#define NEIGHBOR_C 2
#define PLOT 3
#define FINAL 5
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
     #ifdef _MPI_
     double partial_sum = l2norm;
     *_mx = mx;
     return partial_sum;
     #else
     l2norm /= (double) ((m+1)*(n+1));
     l2norm = sqrt(l2norm);

     *_mx = mx;
     return l2norm;
     #endif
 }

void coord (int px, int py, int procID, int m, int n, int &cbeg, int &cend, int &rbeg, int &rend) {
    int rinterval = (m+1) / py;
    int cinterval = (n+1) / px;
    int rowID = procID / px;
    int colID = procID % px;
    int rrem = (m+1) % py ;
    int crem = (n+1) % px ;
   //calculating the beginning and ending global indexes of the column 
    if ( colID < crem ){
	cbeg = colID*(cinterval + 1 ) + 1;//deadly one:it's starting from 1 !
	cend = cbeg + (cinterval + 1)-1;
    }
    else {
	cbeg = crem * (cinterval + 1)+((colID+1-crem)-1)*cinterval + 1;
	cend = cbeg + cinterval-1;
    }
   // calculating the beginning and ending global indexes of the row 
    if ( rowID < rrem ){
	rbeg = rowID*(rinterval + 1 ) + 1;//deadly one:it's starting from 1 !
	rend = rbeg + (rinterval + 1)-1;
    }
    else {
	rbeg = rrem * (rinterval + 1)+((rowID+1-rrem)-1)*rinterval + 1;
	rend = rbeg + rinterval-1;
    }
}
// Added px and py to solve
int solve(ofstream& logfile, double ***_Ew, double ***_E, double ***_E_prev, double **R, int m, int n, int niters, double alpha, double dt, int plot_freq, Plotter *plotter, int stats_freq, int px, int py,int allm = -1, int alln = -1, bool noComm= false){

 // Simulated time is different from the integer timestep number
 double t = 0.0;
 int root = 0;
 double **E = *_E, **E_prev = *_E_prev, **Ew=*_Ew;
 int niter;
 int rank =0, np=1;
 #ifdef _MPI_
 MPI_Comm_size(MPI_COMM_WORLD,&np);
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Request	send_request,recv_request;
 MPI_Request	send_request2,recv_request2;
 MPI_Request	send_requestc,recv_requestc;
 MPI_Request	send_request2c,recv_request2c;
 int plot_msgsiz = (m+1)*(n+1)*sizeof(double);
 char * msgbuffer = new char[plot_msgsiz];
 int col_msgsiz = (m+1)*sizeof(double);
 double * sendbuffer = new double [m+1];
 double * recvbuffer = new double [m+1];
 double * sendbuffer2 = new double [m+1];
 double * recvbuffer2 = new double [m+1];

 #else
 m=allm;
 n=alln;
 #endif
 // We continue to sweep over the mesh until the simulation has reached
 // the desired simulation Time
 // This is different from the number of iterations

/* Setting up our boundaries for i, not pretty but it works */
/*
int i_s = (((n+1)*rank)/py)+1;
i_s += rank==0?0:1;
int i_e = rank==py-1?n+1:(((n+1)*(rank+1))/py)+1;
*/
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
   
   // Solve for the excitation, a PDE 
   int colid = rank % px;
   int rowid = rank / px;
   
   for (j=1; j<=m+1; j++){
     E_prev[j][0] = E_prev[j][2];
     E_prev[j][n+2] = E_prev[j][n];
   }
     for (i=1; i<=n+1; i++) {
	E_prev[0][i] = E_prev[2][i];
        E_prev[m+2][i] = E_prev[m][i];
   }
	/* 0 | 1 | 2	0 | 1 | 2 
	 *-----------	---------
	 * 3 | 4 | 5	3 | 4 | 5 
	 *-----------	---------
	 * 6 | 7 | 8	6 | 7 | 8 = 2
     * col 0   1   2                    */

	#ifdef _MPI_
    if (!noComm) { // Turn off communication if this is true.
	if (py > 1) {
		/* Row Message Passing */
		// Top
		if (rowid == 0) {
			/* This sends the boundary of the first thread (the items at the end of first thread) */
			MPI_Isend(&E_prev[m+1][1], n+1 , MPI_DOUBLE, rank + px , NEIGHBOR, MPI_COMM_WORLD, &send_request);
			/* This recieves the boundary of the second thread (i+1) */
			MPI_Irecv(&E_prev[m+2][1], n+1, MPI_DOUBLE, rank + px, NEIGHBOR, MPI_COMM_WORLD, &recv_request);
			MPI_Wait (&send_request,MPI_STATUS_IGNORE);
			MPI_Wait (&recv_request,MPI_STATUS_IGNORE);
		}
		// bottom
		else if (rowid == py-1) {
			/* This sends the boundary of the second thread (the items at the end of the second thread) */
			MPI_Isend(&E_prev[1][1], n+1, MPI_DOUBLE, rank - px, NEIGHBOR, MPI_COMM_WORLD, &send_request);
			/* This recieves the boundary of the first thread needed for computation (i-1) */
			MPI_Irecv(&E_prev[0][1], n+1, MPI_DOUBLE, rank - px, NEIGHBOR, MPI_COMM_WORLD, &recv_request);
			MPI_Wait (&send_request,MPI_STATUS_IGNORE);
			MPI_Wait (&recv_request,MPI_STATUS_IGNORE);
		}
		 else {
		//middle
			/* This sends the boundary of the first thread (the items at the end of first thread) */
			MPI_Isend(&E_prev[m+1][1], n+1 , MPI_DOUBLE, rank + px , NEIGHBOR, MPI_COMM_WORLD, &send_request);
			/* This recieves the boundary of the second thread (i+1) */
			MPI_Irecv(&E_prev[m+2][1], n+1, MPI_DOUBLE, rank + px, NEIGHBOR, MPI_COMM_WORLD, &recv_request);
			MPI_Wait (&send_request,MPI_STATUS_IGNORE);
			MPI_Wait (&recv_request,MPI_STATUS_IGNORE);

			MPI_Isend(&E_prev[1][1], n+1, MPI_DOUBLE, rank - px, NEIGHBOR, MPI_COMM_WORLD, &send_request);
			/* This recieves the boundary of the first thread needed for computation (i-1) */
			MPI_Irecv(&E_prev[0][1], n+1, MPI_DOUBLE, rank - px, NEIGHBOR, MPI_COMM_WORLD, &recv_request);
			MPI_Wait (&send_request,MPI_STATUS_IGNORE);
			MPI_Wait (&recv_request,MPI_STATUS_IGNORE);
		} 
	}
	if (px > 1) {
	/* Col Message Passing */
		
	// Leftmost
		if (colid == 0) {
			///cout << "Rank " << rank << endl;
			/* Packing to send right */
			int position=0;
			double chk=0;
			for (int j=1; j<=m+1; j++){
			    sendbuffer[j-1] = E_prev[j][n+1];
			    //MPI_Pack(&E_prev[j][n+1],1, MPI_DOUBLE, sendbuffer, col_msgsiz, &position, MPI_COMM_WORLD);
			  //  chk += E_prev[j][n+1];
			}
			//cout << "Rank sending " << rank << "|| chk = " << chk << endl; 
			MPI_Isend(&sendbuffer[0], m+1, MPI_DOUBLE,rank+1,NEIGHBOR_C, MPI_COMM_WORLD, &send_requestc);
			MPI_Irecv(&recvbuffer[0], m+1, MPI_DOUBLE, rank+1,NEIGHBOR_C, MPI_COMM_WORLD, &recv_requestc);
			// Recieve packed array
			position=0; // Reset position
			// Unpack in ghost cells
			chk=0;
		//	cout << "Rank recieving " << rank << "|| chk = " << chk << endl;
			MPI_Wait (&recv_requestc,MPI_STATUS_IGNORE);
			for (int j=1;j<=m+1;j++,position++) {
				  E_prev[j][n+2] = recvbuffer[j-1];
				 // chk += recvbuffer[j-1];
				  //MPI_Unpack(recvbuffer, col_msgsiz, &position, &E_prev[j][n+2], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			MPI_Wait (&send_requestc,MPI_STATUS_IGNORE);

		}
		// Rightmost
		if (colid == px-1) {
			int position =0;
			/* Packing to send left */
			double chk = 0;
			for (int j=1; j<=m+1; j++, position++){
			    sendbuffer[j-1] = E_prev[j][1];
			  //  chk += E_prev[j][1];
			    //MPI_Pack(&E_prev[j][1],1, MPI_DOUBLE, sendbuffer, col_msgsiz, &position, MPI_COMM_WORLD);
			}
			//cout << "Rank sending " << rank << "|| chk = " << chk << endl;
			MPI_Isend(&sendbuffer[0], m+1, MPI_DOUBLE,rank-1,NEIGHBOR_C, MPI_COMM_WORLD, &send_requestc);
			MPI_Irecv(&recvbuffer[0], m+1, MPI_DOUBLE, rank-1,NEIGHBOR_C, MPI_COMM_WORLD, &recv_requestc);
			// Recieve packed array
			position=0; // Reset position
			// Unpack in ghost cells
			
			//cout << "Rank recieving " << rank << "|| chk = " << chk << endl;
			MPI_Wait (&recv_requestc,MPI_STATUS_IGNORE);
			chk =0;
			for (int j=1; j<=m+1; j++, position++){
			    E_prev[j][0] = recvbuffer[j-1];
			//    chk += recvbuffer[j-1];
			    //MPI_Unpack(recvbuffer, col_msgsiz, &position, &E_prev[j][0], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			MPI_Wait (&send_requestc,MPI_STATUS_IGNORE);
		}
		if (colid>0 && colid <px-1)  {
		//middle
			/* Packing to send right */
			for (int j=1; j<=m+1; j++){
			    sendbuffer[j-1] = E_prev[j][n+1];
			}
			MPI_Isend(&sendbuffer[0], m+1, MPI_DOUBLE,rank+1,NEIGHBOR_C, MPI_COMM_WORLD, &send_requestc);
			MPI_Irecv(&recvbuffer[0], m+1, MPI_DOUBLE, rank+1,NEIGHBOR_C, MPI_COMM_WORLD, &recv_requestc);
			// Recieve packed array
			MPI_Wait (&recv_requestc,MPI_STATUS_IGNORE);
			// Unpack in ghost cells
			for (int j=1;j<=m+1;j++) {
				  E_prev[j][n+2] = recvbuffer[j-1];
			}
			MPI_Wait (&send_requestc,MPI_STATUS_IGNORE);
			// Unpack in ghost cells
			
			/* Packing to send left */
			for (int j=1; j<=m+1; j++){
			    sendbuffer[j-1] = E_prev[j][1];
			}
			MPI_Isend(&sendbuffer[0], m+1, MPI_DOUBLE,rank-1,NEIGHBOR_C, MPI_COMM_WORLD, &send_request2c);
			// Unpack in ghost cells
			MPI_Irecv(&recvbuffer[0], m+1, MPI_DOUBLE, rank-1,NEIGHBOR_C, MPI_COMM_WORLD, &recv_request2c);
			// Unpack in ghost cells
			MPI_Wait (&recv_request2c,MPI_STATUS_IGNORE);
			// Unpack in ghost cells
			// Recieve packed array
			for (int j=1; j<=m+1; j++){
			    E_prev[j][0] = recvbuffer[j-1];
			}
			MPI_Wait (&send_request2c,MPI_STATUS_IGNORE);
			// Unpack in ghost cells
		}
		
	}
   }
	#endif
     for (j=1; j<=m+1; j++){
	     for (i=1; i<=n+1; i++) {;
		E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
	     }
     

    }

   /* 
    * Solve the ODE, advancing excitation and recovery variables
    *     to the next timtestep
    */
   for (j=1; j<=m+1; j++){
     double *RR = &R[j][1];
     double *EE = &E[j][1];
     for (i=1; i<=n+1; i++, EE++,RR++) {
	EE[0] += -dt*(kk*EE[0]*(EE[0]-a)*(EE[0]-1)+EE[0]*RR[0]);
	RR[0] += dt*(epsilon+M1* RR[0]/( EE[0]+M2))*(-RR[0]-kk*EE[0]*(EE[0]-b-1));
     }
   }


   if (stats_freq){
     #ifdef _MPI_
     
     double mx;
     double partial_sum = stats(E_prev,m,n,&mx); 
     double total;
     double max;
     double l2norm;
   
     MPI_Reduce(&mx,&max,1,MPI_DOUBLE,MPI_MAX,root,MPI_COMM_WORLD); 
     MPI_Reduce(&partial_sum, &total,1,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
	if ( rank == root ){
		l2norm = total / ((allm+1)*(alln+1));
		l2norm = sqrt(l2norm);
	}
	repNorms(logfile,l2norm,max,dt,allm,alln,niter, stats_freq);

     #else
     double mx;
     double l2norm = stats(E_prev,m,n,&mx);
     repNorms(logfile,l2norm,mx,dt,m,n,niter, stats_freq);
     #endif
     
  
   }

   if (plot_freq){
          if (!(niter % plot_freq)){
#ifdef _MPI_
		MPI_Request sendrequestp, recvrequestp;
		// splot(E,niter,m+1,n+1,WAIT);
		/* *************************** */
		/* Gather matrix for plot      */
		/* *************************** */

		/* We're going to let thread 0 do the return value, so that means every other thread
		* needs to pack up it's own working segment and send it to thread 0 */
		
		int plot_pos = 0;

		if (rank != root) {	
			double* submatrix = new double[(m+1)*(n+1)];
			int pmsgsize = (m+1)*(n+1);
			for (int j=1; j<=m+1; j++){
			    //first 0-n, seconde n+1-2n+1
			    for (int i = 1; i<=n+1;i++){
				     submatrix[(j-1)*(n+1)+(i-1)] = E[j][i];	
				}
			   // MPI_Pack(&E[j][1],n+1, MPI_DOUBLE, msgbuffer, plot_msgsiz, &plot_pos, MPI_COMM_WORLD);
			}
				/* Then we send to thread 0 */
			MPI_Isend(submatrix, pmsgsize, MPI_DOUBLE,root,PLOT, MPI_COMM_WORLD, &sendrequestp);
			MPI_Wait (&sendrequestp,MPI_STATUS_IGNORE);
			if (submatrix != NULL) delete [] submatrix;

		}
		else { 
			// IN THREAD 0
			int cbeg,cend,rbeg,rend;
			MPI_Status packstat;
			int proc_seen=0;
			// Copy the matrix in thread 0 first to Large Matrix.
			for (int j=1;j<=m+1;j++) 
			  for (int i=1;i<=n+1;i++) 
			    Ew[j][i]=E[j][i];

			/* We loop until we've recieved communication from every thread */
			while (proc_seen < np-1) {
				plot_pos=0; // Reset buffer iterator each time we recieve a message
				/* We recieve our packed stuff from some other thread */
				int revmsgsize = (m+2)*(n+2);
				double* revmatrix = new double[revmsgsize];
				MPI_Irecv(revmatrix, revmsgsize, MPI_DOUBLE, MPI_ANY_SOURCE,PLOT, MPI_COMM_WORLD, &recvrequestp);
				MPI_Wait (&recvrequestp, &packstat);
				int recv_procId = packstat.MPI_SOURCE;
				coord (px, py, recv_procId, allm, alln, cbeg, cend, rbeg, rend);
				/* We unpack the rank first */
				int colsize = cend-cbeg+1;
				for (int j=rbeg;j<=rend;j++) {
					  int localrow = j-rbeg;
					  for (int i = cbeg ; i<=cend;i++){
						int localcol = i-cbeg;
						Ew[j][i]=revmatrix[localrow*colsize+localcol];
					  }
				}
				proc_seen++; // Move to the next process
			}

			plotter->updatePlot(Ew,  niter, allm+1, alln+1, WAIT);
		}
#else
		plotter->updatePlot(E,  niter, m+1, n+1, WAIT);
#endif
	}
    }

   // Swap current and previous
   double **tmp = E; E = E_prev; E_prev = tmp;
 }

  // Store them into the pointers passed in
  #ifdef _MPI_
  delete [] msgbuffer;
  delete [] sendbuffer;
  delete [] recvbuffer;
  delete [] sendbuffer2;
  delete [] recvbuffer2;  
  #endif
  *_E = E;
  *_E_prev = E_prev;
  return niter;

}
