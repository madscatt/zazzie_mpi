
# cluster (untested)

DIRINC_BOOST=/share/apps/local/include/boost/
DIRLIB_BOOST=/share/apps/local/lib/

DIRINC=/share/apps/local/core_libraries/include/
DIRLIB=/share/apps/local/core_libraries/lib/
DIRLIB_OACC=./openacc/

CC=/share/apps/local/bin/g++

# mac (tested)

#DIRINC_BOOST=/usr/local/include/
#DIRLIB_BOOST=/usr/local/lib/

#DIRINC=/usr/local/core_libraries/include/
#DIRLIB=/usr/local/core_libraries/lib/
#DIRLIB_OACC=./openacc/
#CC=g++

# set flags

CCFLAGS=-std=c++11 -c -fPIC -w -g -I$(DIRINC) -I$(DIRINC_BOOST)
LDFLAGS=-std=c++11 -g -w -I$(DIRINC) -I$(DIRINC_BOOST)

# set sasmol library
#
LIB=sasmol
LIBRARY=$(DIRLIB)lib$(LIB).a

# set openacc library
#
LIB_OACC=oacc_pr
LIBRARY_OACC=$(DIRLIB_OACC)lib$(LIB_OACC).a

# compile pr_mpi

pr_mpi.out: pr_mpi.cpp $(LIBRARY) $(LIBRARY_OACC)
	$(CC) $< $(LDFLAGS) -L$(DIRLIB_BOOST) -L$(DIRLIB) -l$(LIB) $(DIRLIB_BOOST)/libboost_regex.a -o $@

clean:
	rm -Rf pr_mpi.out *.o *.so pgprof*

