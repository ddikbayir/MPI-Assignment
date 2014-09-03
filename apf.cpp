/* 
 * Driver for a cardiac elecrophysioly simulatin that uses the
 * Aliev-Panfilov model
 * We use an explicit method
 *
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <iomanip>
#include <string>
#include <math.h>
#include "apf.h"
#include "Plotting.h"
#ifdef _MPI_
#include <mpi.h>
#endif
using namespace std;

// Utilities
// 

// Allocate a 2D array
/*double **alloc2D(int sizeX, int sizeY){
	double** ary = new double*[sizeX];
	for(int i = 0; i < sizeX; ++i)
	    ary[i] = new double[sizeY];
        for (int i = 0;i<sizeX;i++){
		for (int j=0;j<sizeY;j++){
			ary[i][j]=0;
		}
	}
	return ary;
}*/

double **alloc2D(int m,int n){

   double **E;
   int nx=n+1, ny=m+1;
   E = (double**)malloc(sizeof(double*)*ny + sizeof(double)*nx*ny);
   assert(E);
   int j;
   for(j=0;j<ny;j++) E[j] = (double*)(E+ny) + j*nx;
   return(E);
}


void init (double **E,double **E_prev,double **R,int m,int n,int procID = -1, int m_new = -1, int n_new= -1, int px = -1, int py = -1){
    int i,j;
    //----------original--------------
    // Initialization
    /*for (j=1; j<=m + 1; j++)
        for (i=1; i<= n+1; i++){
            E_prev[j][i] = R[j][i] = 0;
    }
    for (j=1; j<=m + 1; j++)
        for (i=n/2+2; i<= n+1 ; i++){
            E_prev[j][i] = 1.0;
    }

    for (j=m/2+2; j<=m+1; j++)
        for (i=1; i<=n+1; i++)
            R[j][i] = 1.0;
    */
    //-----------end original-------
    //-----------myversion--------
    #ifdef _MPI_
    int rinterval = (m+1) / py;
    int cinterval = (n+1) / px;
    int rowID = procID / px;
    int colID = procID % px;
    int rrem = (m+1) % py ;
    int crem = (n+1) % px ;
    for (j=1; j<=m_new + 1; j++)
        for (i=1; i<= n_new+1; i++){
            E_prev[j][i] = R[j][i] = 0;
    }
   //calculating the beginning and ending global indexes of the column 
    int cbeg =0 , cend = 0;
    if ( colID < crem ){
	cbeg = colID*(cinterval + 1 ) + 1;//deadly one:it's starting from 1 !
	cend = cbeg + (cinterval + 1)-1;
    }
    else {
	cbeg = crem * (cinterval + 1)+((colID+1-crem)-1)*cinterval + 1;
	cend = cbeg + cinterval-1;
    }
    for (j=1; j<=m_new + 1; j++){
        for (i=1; i<= n_new+1 ; i++){
	    int colindex = cbeg + i - 1;
	    if (colindex >= n/2+2 && colindex <= n+1)
                E_prev[j][i] = 1.0;
        }
    }
   // calculating the beginning and ending global indexes of the row 
    int rbeg =0 , rend = 0;
    if ( rowID < rrem ){
	rbeg = rowID*(rinterval + 1 ) + 1;//deadly one:it's starting from 1 !
	rend = rbeg + (rinterval + 1)-1;
    }
    else {
	rbeg = rrem * (rinterval + 1)+((rowID+1-rrem)-1)*rinterval + 1;
	rend = rbeg + rinterval-1;
    }
    for ( j = 1 ; j <= m_new + 1;j++){ 
	 //calculate the global index
   	 int rowindex = rbeg + j - 1;
	// if it falls into the range, change it to one
	 if (rowindex >= m/2+2 && rowindex <= m+1){
		for (int i = 1; i<= n+1;i++){
			R[j][i] = 1.0;
		}
	} 
   } 
    #else
    for (j=1; j<=m + 1; j++)
        for (i=1; i<= n+1; i++){
            E_prev[j][i] = R[j][i] = 0;
    }
    for (j=1; j<=m + 1; j++)
        for (i=n/2+2; i<= n+1 ; i++){
            E_prev[j][i] = 1.0;
    }

    for (j=m/2+2; j<=m+1; j++)
        for (i=1; i<=n+1; i++)
            R[j][i] = 1.0;
    #endif
//debug
/*
    string filename = "matrix form process ";
    string ID = to_string(procID); 
    filename = filename + ID ;
    ofstream ofs;
    ofs.open(filename);
    for (j=1; j<=m_new+1 ; j++){
        for (i=1; i<= n_new+1; i++){
           ofs<< E_prev[j][i]<<" " ;
	}
	   ofs<<"\n";
    }
    ofs<<"\n";
    for (j=1; j<=m_new+1 ; j++){
        for (i=1; i<= n_new+1; i++){
           ofs<< R[j][i]<<" " ;
	}
	   ofs<<"\n";
    }
    ofs.close();
*/
    
    //-----------end my version--------
}

// External functions
void cmdLine(int argc, char *argv[], int& n, int& stats_freq, int& plot_freq, int& px, int& py, bool &noComm, int &niters);
int solve(ofstream& logfile, double ***_Ew, double ***_E, double ***_E_prev, double **R, int m, int n, int niters, double alpha, double dt, int plot_freq, Plotter *plotter, int stats_freq, int px , int py,int allm, int alln,bool noComm);
void printTOD(ofstream& logfile, string mesg);
double stats(double **E, int m, int n, double *_mx);
void ReportStart(ofstream& logfile, double dt, int niters, int m, int n, int px, int py, bool noComm);
void ReportEnd(ofstream& logfile, int niters, double l2norm, double mx, int m,int n, double t0, int px, int py, bool noComm);
double getTime();

// Main program
int main(int argc, char** argv)
{
 /*
  *  Solution arrays
  *   E is the "Excitation" variable, a voltage
  *   R is the "Recovery" variable
  *   E_prev is the Excitation variable for the previous timestep,
  *      and is used in time integration
  */
 double **E, **R, **E_prev, **Ew=NULL;

 // Default values for the command line arguments
 int m=100,n=100;
 int stats_freq = 0;
 int plot_freq = 0;
 int px = 1, py = 1;
 int niters=100;
 bool noComm = false;
 int root = 0;

#ifdef _MPI_
 MPI_Init(&argc,&argv);
#endif

// Parse command line arguments
 cmdLine( argc, argv, n, stats_freq,  plot_freq, px, py, noComm, niters);
 if (n < 26){
    cout << "\n *** N must be larger than 25.  Exiting ... " << endl << endl;
    exit(-1);
 }
 m = n;
 int col =0, row = 0;
 int nprocs=1, myrank=0;
#ifdef _MPI_
 MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
 MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#endif

 // The log file
 // Do not change the file name or remove this call
    ofstream logfile;
    if (!myrank)
        logfile.open("Log.txt",ios::out);
//  ofstream logfile("Log.txt",ios::out);
 printTOD(logfile, "Simulation begins");

 // Allocate contiguous memory for solution arrays
 // The computational box is defined on [1:m+1,1:n+1]
 // We pad the arrays in order to facilitate differencing on the 
 // boundaries of the computation box
 // ------------original------------
/* E = alloc2D(m+2,n+2);
 E_prev = alloc2D(m+2,n+2);
 R = alloc2D(m+2,n+2);

 init(E,E_prev,R,m,n);
*/
//--------------original-----------
//
//
//--------my version------
#ifdef _MPI_
 int ny = py;
 int nx = px;
 int rrem = (m+1) % ny;
 int crem = (n+1) % nx;
 /*
 cout << "py is " <<py<<" rem is " <<rem<<endl;
 cout << "myrank "<<myrank <<endl;
 */
 if( myrank == 0 ) {
        
	Ew = alloc2D( m+2, n+2 );
}
 int rowID= myrank / nx;
 int colID= myrank % nx;
 row =0;
 if ( rowID < rrem )
     row = (m+1)/ny + 1  - 1;// -1 to fit in init
 //total points should be row + 1;
 else 
     row = (m+1)/ny - 1;
 //total points should be row
 
 col = 0;
 if ( colID < crem )
     col = (n+1)/nx + 1  - 1;// -1 to fit in init
 //total points should be col + 1;
 else 
     col = (n+1)/nx - 1;
 //total points should be col

 
 //cout<<"row "<<row <<" col "<<col << " rank " << myrank <<endl;
 
 E = alloc2D(row+100,col+100);//ghose cells are included in the 2!
 E_prev = alloc2D(row+100,col+100);
 R = alloc2D(row+100,col+100);
 init(E,E_prev,R,m,n,myrank,row,col,px,py); //myrank starts from 0
#else
 E = alloc2D(m+2,n+2);
 E_prev = alloc2D(m+2,n+2);
 R = alloc2D(m+2,n+2);

 init(E,E_prev,R,m,n);
#endif
/*debug
    string filename = "data of matrix form process ";
    string ID = to_string(myrank); 
    filename = filename + ID ;
    ofstream ofs;
    ofs.open(filename);
    ofs<< "myrank " << myrank<<endl;
    ofs<< "crem" <<crem<<endl;
    ofs<< "rrem" <<rrem<<endl;
    ofs<< "col " << col <<endl;
    ofs<< "row " << row <<endl;
    ofs<< "rowID " << rowID <<endl;
    ofs<< "colID " << colID <<endl;
    ofs.close();
*/

//-------------end my version

 //
 // Initization of various simulation variables
 // Do not change the code these assignments statements, as doing so
 // could cause your submission to be graded incorrectly
 //

 // We compute the timestep dt

 // We should always use double precision values for the folowing variables:
 //    rp, dte, dtr, ddt
 //
 // This ensures that the computation of dte and especially dt
 // will not lose precision (i.e. if computed as single precision values)

 double dx = 1.0/n;
 double rp= kk*(b+1)*(b+1)/4;
 double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
 double dtr=1/(epsilon+((M1/M2)*rp));
 double ddt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
 double dt = (double) ddt;
 double alpha = d*dt/(dx*dx);

 // End Initization of various simulation variables

 // Report various information
 // Do not remove this call, it is needed for grading
  ReportStart(logfile, dt, niters, m, n, px, py, noComm);

 Plotter *plotter = NULL;
 if (plot_freq){
     plotter = new Plotter();
     assert(plotter);
 }

 // Start the timer
#ifdef _MPI_
 double t0 = -MPI_Wtime();
#else
 double t0 = -getTime();
#endif
 int niter = solve(logfile, &Ew, &E, &E_prev, R, row, col, niters, alpha, dt, plot_freq, plotter, stats_freq, px, py ,m , n, noComm);

#ifdef _MPI_
 t0 += MPI_Wtime();
#else
 t0 += getTime();
#endif
/*    string filename = "data of matrix form process ";
    string ID = to_string(myrank); 
    filename = filename + ID ;
    ofstream ofs;
    ofs.open(filename);
    ofs<< "myrank " << myrank<<endl;
    ofs.close();
*/
if (niter != niters)
   cout << "*** niters should be equal to niters" << endl;
 // Report various information
 // Do not remove this call, it is needed for grading    
     #ifdef _MPI_
     
     double mx;
     double partial_sum = stats(E_prev,row,col,&mx); 
     double total;
     double max;
     double l2norm;
     
     MPI_Reduce(&mx,&max,1,MPI_DOUBLE,MPI_MAX,root,MPI_COMM_WORLD); 
     MPI_Reduce(&partial_sum, &total,1,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
	if ( myrank == root ){
		l2norm = total / ((m+1)*(n+1));
		l2norm = sqrt(l2norm);
	}
    ReportEnd(logfile,niters,l2norm,max,m,n,t0,px, py, noComm);
     #else

    double mx;
    double l2norm = stats(E_prev,m,n,&mx);
    ReportEnd(logfile,niters,l2norm,mx,m,n,t0,px, py, noComm);
     #endif

 if (plot_freq){
    cout << "\n\nEnter any input to close the program and the plot...";
    int resp;
    cin >> resp;
  }

 if (!myrank) {
    logfile.close();
    free (Ew);
 }
 free (E);
 free (E_prev);
 free (R);
 if (plot_freq)
     delete plotter;
#ifdef _MPI_
 MPI_Finalize();
#endif
}
