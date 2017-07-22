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
    float x[nframes][natoms] ;


    MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Get_processor_name(processor_name , &name_length) ;

    if(my_rank == 0){
        for(int i=0 ; i < nframes ; i++){
            for(int j = 0 ; j < natoms ; j++){
                x[i][j] = float(i*10) + float(j) ;
            }
        }
        printf("%s\n", "rank zero showing typical values") ;
	    printf("x[0][0] = %f\n", x[0][0]);
	    printf("x[0][1] = %f\n\n", x[0][1]);
	    printf("x[2][0] = %f\n\n", x[2][0]);
	    printf("x[0][8] = %f\n", x[0][8]);
	    printf("x[0][9] = %f\n", x[0][9]);
	    printf("x[9][8] = %f\n", x[9][8]);
	    printf("x[9][9] = %f\n", x[9][9]);
    }

    if(my_rank == 0){
        MPI_Send(&(x[2][0]), 20, MPI_FLOAT, 1, 0, MPI_COMM_WORLD) ;
    }
    else{
        MPI_Recv(&(x[0][0]), 20, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status) ;
        
	    printf("my_rank = %i\tx[0][0] = %f\n", my_rank, x[0][0]);
	    printf("my_rank = %i\tx[0][1] = %f\n", my_rank, x[0][1]);
	    printf("my_rank = %i\tx[0][2] = %f\n", my_rank, x[0][2]);
	    printf("my_rank = %i\tx[0][3] = %f\n", my_rank, x[0][3]);
	    printf("my_rank = %i\tx[1][0] = %f\n", my_rank, x[1][0]);
	    printf("my_rank = %i\tx[1][1] = %f\n", my_rank, x[1][1]);
	    printf("my_rank = %i\tx[1][2] = %f\n", my_rank, x[1][2]);
	    printf("my_rank = %i\tx[1][3] = %f\n", my_rank, x[1][3]);
	    printf("my_rank = %i\tx[2][0] = %f\n", my_rank, x[2][0]);
	    printf("my_rank = %i\tx[2][1] = %f\n", my_rank, x[2][1]);
	    printf("my_rank = %i\tx[2][2] = %f\n", my_rank, x[2][2]);
	    printf("my_rank = %i\tx[2][3] = %f\n", my_rank, x[2][3]);
    }


	MPI_Finalize();
}



