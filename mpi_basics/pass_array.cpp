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

    int nframes = 10 ; 
    int natoms = 10 ; 
    float x[nframes][natoms] ;

    for(int i=0 ; i < nframes ; i++){
        for(int j = 0 ; j < natoms ; j++){
            x[i][j] = float(i*j) ;
        }
    }

    printf("%s\n", "outside mpi") ;
	printf("%f\n", x[0][1]);
	printf("%f\n", x[1][1]);

    MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Get_processor_name(processor_name , &name_length) ;

    printf("my_rank = %d\t number_of_processes = %d\n", my_rank, number_of_processes) ;
    printf("processor_name = %s\n", processor_name) ;

    //MPI_Datatype column_mpi_t ;

    //MPI_Type_vector(10, 1, np, MPI_FLOAT, &column_mpi_t) ;
    //MPI_Type_commit(&column_mpi_t);

    if(my_rank == 0){
	    //sprintf(message, "Hello from process %d!", my_rank);
	    dest = 0;
	    //MPI_Send(message, strlen(message) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
	    for(source = 1; source < number_of_processes; source ++) {
            MPI_Send(&(x[source][0]), natoms, MPI_FLOAT, source, 0, MPI_COMM_WORLD) ;
        };
	    //MPI_Send(&(x[0][0]), 1, column_mpi_t, 1, 0, MPI_COMM_WORLD) ;
    }
    else{
	    //for(source = 1; source < number_of_processes; source ++) {
	        //MPI_Recv(message, 100, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
	    for(source = 1; source < number_of_processes; source ++) {
	        MPI_Recv(&(x[source][0]), natoms, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status) ;
	        //MPI_Recv(&(x[0][0]), 1, column_mpi_t, 0, 0, MPI_COMM_WORLD, &status) ;
	        //printf("%s\n", message);
	        printf("%f\n", x[source][1]);
        };
	//} ; 
	} ; 

	MPI_Finalize();
}



