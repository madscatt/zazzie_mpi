#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include <time.h>
 
void sleep(unsigned int mseconds)
{
    clock_t goal = mseconds + clock();
    while (goal > clock());
}

int main(int argc, char* argv[]){
	int my_rank;
	int number_of_processes;
	int source;
	int dest;
	int tag = 0;
	char message[100];
	MPI_Status status;
    char processor_name[MPI_MAX_PROCESSOR_NAME] ;
    int name_length ;

    int nframes = 4 ; 
    int natoms = 10 ; 
    double x[nframes][natoms] ;
    double lx[nframes/2][natoms] ;

    MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Get_processor_name(processor_name , &name_length) ;

    if(my_rank == 0){
        for(int i=0 ; i < nframes ; i++){
            for(int j = 0 ; j < natoms ; j++){
                x[i][j] = double(i*10) + double(j) ;
            }
        }
        printf("%s\n", "rank zero showing typical values") ;
	    printf("x[0][0] = %lf\n", x[0][0]);
	    printf("x[0][1] = %lf\n\n", x[0][1]);
	    printf("x[2][0] = %lf\n\n", x[2][0]);
	    printf("x[0][8] = %lf\n", x[0][8]);
	    printf("x[0][9] = %lf\n", x[0][9]);
	    printf("x[3][8] = %lf\n", x[3][8]);
	    printf("x[3][9] = %lf\n", x[3][9]);
    }
    else{
        for(int i=0 ; i < nframes/2 ; i++){
            for(int j = 0 ; j < natoms ; j++){
                lx[i][j] = double(0.0) ;
            }
        }
    }

    if(my_rank == 0){
        MPI_Send(&(x[2][0]), 20, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD) ;
    }
    else{
        MPI_Recv(&(lx[0][0]), 20, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status) ;
        
	    printf("my_rank = %i\tlx[0][0] = %lf\n", my_rank, lx[0][0]);
	    printf("my_rank = %i\tlx[0][1] = %lf\n", my_rank, lx[0][1]);
	    printf("my_rank = %i\tlx[0][2] = %lf\n", my_rank, lx[0][2]);
	    printf("my_rank = %i\tlx[0][3] = %lf\n", my_rank, lx[0][3]);
	    printf("my_rank = %i\tlx[1][0] = %lf\n", my_rank, lx[1][0]);
	    printf("my_rank = %i\tlx[1][1] = %lf\n", my_rank, lx[1][1]);
	    printf("my_rank = %i\tlx[1][2] = %lf\n", my_rank, lx[1][2]);
	    printf("my_rank = %i\tlx[1][3] = %lf\n", my_rank, lx[1][3]);
    }


	MPI_Finalize();
}



