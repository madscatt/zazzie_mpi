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

//    printf("hello, I am running and ready for mpi stuff\n") ;

    MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Get_processor_name(processor_name , &name_length) ;

    printf("my_rank = %d\t number_of_processes = %d\n", my_rank, number_of_processes) ;
    printf("processor_name = %s\n", processor_name) ;

    if(my_rank != 0){
	    sprintf(message, "Hello from process %d!", my_rank);
	    dest = 0;
	    MPI_Send(message, strlen(message) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    }
    else{
	    for(source = 1; source < number_of_processes; source ++) {
	        MPI_Recv(message, 100, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
	        printf("%s\n", message);
	} ; 
	} ; 

	MPI_Finalize();
}



