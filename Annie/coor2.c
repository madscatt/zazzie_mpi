
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

/*
 This is a version with integers, rather than char arrays, presented in this
 very good answer: http://stackoverflow.com/a/9271753/2411320
 It will initialize the 2D array, scatter it, increase every value by 1 and then gather it back.
*/

int malloc2D(double ***array, int n, int m) {
    int i;
    /* allocate the n*m contiguous items */
    double *p = malloc(n*m*sizeof(double));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = malloc(n*sizeof(double*));
    if (!(*array)) {
       free(p);
       return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (i=0; i<n; i++)
       (*array)[i] = &(p[i*m]);

    return 0;
}

int free2D(double ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}

int main(int argc, char **argv) {
    double **global, **local, **globaly, **localy, **globalz, **localz;
    const int xn_atoms=2; // number of rows or atoms for x coor
    const int xn_frames=4; // number of columns or frames for x coor
    const int xproc_frames=2;  // size of process grid, how many frames packaged at a time
    int rank, size;        // rank of current process and no. of processes
    int i, j, p;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (size != xproc_frames) {
        fprintf(stderr,"%s: Only works with np=%d for now\n", argv[0], xproc_frames);
        MPI_Abort(MPI_COMM_WORLD,1);
    }


    if (rank == 0) {
        /* fill in the array, and print it */
        malloc2D(&global, xn_atoms, xn_frames);
        double counter = 0.0;
        for (i=0; i< xn_atoms; i++) {
            for (j=0; j< xn_frames; j++)
                global[i][j] = ++counter;
        }


        printf("Global x array is:\n");
        for (i=0; i< xn_atoms; i++) {
            for (j=0; j< xn_frames; j++) {
                printf("%2g ", global[i][j]);
            }
            printf("\n");
        }

	malloc2D(&globaly, xn_atoms, xn_frames);
        counter = 0.0;
        for (i=0; i< xn_atoms; i++) {
            for (j=0; j< xn_frames; j++)
                globaly[i][j] = ++counter +1.0;
        }


        printf("Global y array is:\n");

        for (i=0; i< xn_atoms; i++) {
            for (j=0; j< xn_frames; j++) {
                printf("%2g ", globaly[i][j]);
            }
            printf("\n");
        }

	malloc2D(&globalz, xn_atoms, xn_frames);
        counter = 0.0;
        for (i=0; i< xn_atoms; i++) {
            for (j=0; j< xn_frames; j++)
                globalz[i][j] = ++counter +2.0;
        }


        printf("Global z array is:\n");
        for (i=0; i< xn_atoms; i++) {
            for (j=0; j< xn_frames; j++) {
                printf("%2g ", globalz[i][j]);
            }
            printf("\n");
        }
    }
    //return;

    /* create the local array which we'll process */
    malloc2D(&local, xn_atoms, xn_frames/xproc_frames);
    malloc2D(&localy, xn_atoms, xn_frames/xproc_frames);
    malloc2D(&localz, xn_atoms, xn_frames/xproc_frames);



    /* create a datatype to describe the subarrays of the global array */
    int sizes[2]    = {xn_atoms, xn_frames};         /* global size */
    int subsizes[2] = {xn_atoms, xn_frames/xproc_frames};     /* local size */
    int starts[2]   = {0,0};                        /* where this one starts */
    MPI_Datatype type, subarrtype;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
    MPI_Type_create_resized(type, 0, xn_frames/xproc_frames *sizeof(double), &subarrtype);
    MPI_Type_commit(&subarrtype);

    double *globalptr;
    if (rank == 0)
        globalptr = &(global[0][0]);

    double *globalptry;
    if (rank == 0)
        globalptry = &(globaly[0][0]);

    double *globalptrz;
    if (rank == 0)
        globalptrz = &(globalz[0][0]);

    /* scatter the array to all processors */
    int sendcounts[xn_atoms*xproc_frames];
    int displs[xn_atoms*xproc_frames];

    if (rank == 0) {
        for (i=0; i<xn_atoms*xproc_frames; i++)
            sendcounts[i] = 1;
        int disp = 0;
        for (i=0; i<xn_atoms; i++) {
            for (j=0; j<xproc_frames; j++) {
                displs[i*xproc_frames+j] = disp;
                disp += 1;
            }
            disp += ((xn_frames/xproc_frames)-1)*xproc_frames;
        }
    }


    MPI_Scatterv(globalptr, sendcounts, displs, subarrtype, &(local[0][0]),
                 xn_atoms*xn_frames/(xproc_frames), MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    MPI_Scatterv(globalptry, sendcounts, displs, subarrtype, &(localy[0][0]),
                 xn_atoms*xn_frames/(xproc_frames), MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    MPI_Scatterv(globalptrz, sendcounts, displs, subarrtype, &(localz[0][0]),
                 xn_atoms*xn_frames/(xproc_frames), MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    /* now all processors print their local data: */

    for (p=0; p<size; p++) {
        if (rank == p) {
            printf("Local process on rank %d is:\n", rank);
            for (i=0; i<xn_atoms; i++) {
                putchar('|');
                for (j=0; j<xn_frames/xproc_frames; j++) {
                    printf("%2g ", local[i][j]);
                }
                printf("|\n");
            }

	    for (i=0; i<xn_atoms; i++) {
                putchar('|');
                for (j=0; j<xn_frames/xproc_frames; j++) {
                    printf("%2g ", localy[i][j]);
                }
                printf("|\n");
            }

	    for (i=0; i<xn_atoms; i++) {
                putchar('|');
                for (j=0; j<xn_frames/xproc_frames; j++) {
                    printf("%2g ", localz[i][j]);
                }
                printf("|\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free2D(&local);

    
    MPI_Type_free(&subarrtype);

        free2D(&global);



    MPI_Finalize();

    return 0;
}
