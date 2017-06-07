# Annie's zazzie_mpi programs


### NOTES from JC

mpicc -o hello.exe hello.c

mpirun -n 3 --hostfile hosts hello.exe

where the file "hosts" has three lines

onsager
compute-0-0
compute-0-1

