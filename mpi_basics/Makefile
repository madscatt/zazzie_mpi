
# onsager
DIRINC=/usr/include/mpich/
DIRLIB=/usr/lib/mpich/

#entropy
DIRINC=/usr/lib64/openmpi/include
DIRLIB=/usr/lib64/openmpi/

#mac
DIRINC=/opt/pgi/osx86-64/2016/mpi/mpich/include
DIRLIB=/opt/pgi/osx86-64/2016/mpi/mpich/lib

CC=mpic++
CC=/opt/pgi/osx86-64/2016/mpi/mpich/bin/mpic++

CCFLAGS=-std=c++11 -c -fPIC -w -g -I$(DIRINC) 
LDFLAGS=-std=c++11 -g -w -I$(DIRINC) 

all: pass_array.o
	$(CC) pass_array.o -o pass_array.exe

pass_array.o: pass_array.cpp
	$(CC) -c pass_array.cpp $(CCFLAGS) -I$(DIRINC)  

clean:
	rm -Rf *.o *.exe
