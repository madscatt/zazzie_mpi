#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char* argv[]){
    
	int my_rank;
	int number_of_processes;
	int source;
	int destination;
	int tag = 0;
	char message[100];
	MPI_Status status;

    char processor_name[MPI_MAX_PROCESSOR_NAME] ;
    int name_length ;

    MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Get_processor_name(processor_name , &name_length) ;

    if(my_rank != 0){
	    sprintf(message, "Hello from process rank = %d : running on %s", my_rank, processor_name);
	    destination = 0;
	    MPI_Send(message, strlen(message) + 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD);
    }
    else{
	    sprintf(message, "RANK 0 : Hello from process rank = %d : running on %s", my_rank, processor_name);
	    printf("%s\n", message);
	    for(source = 1; source < number_of_processes; source ++) {
	        MPI_Recv(message, 100, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
	        printf("%s\n", message);
	} ; 
	} ; 

	MPI_Finalize();
}



