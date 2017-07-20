multicore is working fine (pre mpi)
gpu runs fine, but only one p(r) is recorded, there is something wrong with the copy data not working correctly, note that parallel_hist was defined in order to rule out <vector> hist ... once the data management is worked out go back and see if we can do this with <vector> hist or not.


