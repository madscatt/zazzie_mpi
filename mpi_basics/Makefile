DIRINC=/usr/include/mpich/
DIRLIB=/usr/lib/mpich/

CC=mpicc

CCFLAGS=-std=c++11 -c -fPIC -w -g -I$(DIRINC) 
LDFLAGS=-std=c++11 -g -w -I$(DIRINC) 

all: pass_array.o
	$(CC) pass_array.o -o pass_array.exe

pass_array.o: pass_array.cpp
	$(CC) -c pass_array.cpp $(CCFLAGS) -I$(DIRINC)  

clean:
	rm -Rf *.o *.exe
