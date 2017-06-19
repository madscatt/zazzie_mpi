#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char* argv[]){	
	int global[8];   /* only task 0 has this */
	int local[2];    /* everyone has this */
	const int root = 0;   /* the processor with the initial global data */
	int rank;
	MPI_Init(&argc, &argv);
if (rank == root) {
   for (int i=0; i<7; i++) global[i] = i;


MPI_Scatter(global, 2, MPI_INT,      /* send everyone 2 ints from global */
            local,  2, MPI_INT,      /* each proc receives 2 ints into local */
            root, MPI_COMM_WORLD);   /* sending process is root, all procs in */
                                     /* MPI_COMM_WORLD participate */

printf("%s", str_local);
 }   
	

	

	
	MPI_Finalize;
}
