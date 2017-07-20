
# cluster (untested)

DIRINC_BOOST=/share/apps/local/include/boost/
DIRLIB_BOOST=/share/apps/local/lib/

DIRINC=/share/apps/local/core_libraries/include/
DIRLIB=/share/apps/local/core_libraries/lib/
DIRLIB_OACC=./openacc/

CC=/share/apps/local/bin/g++
CC=/usr/lib64/openmpi/bin/mpic++
#CC=/usr/lib64/mpich/bin/mpic++
#CC=/share/apps/local/pgi/linux86-64/16.10/bin/pgc++

# mac (tested)

#DIRINC_BOOST=/usr/local/include/
#DIRLIB_BOOST=/usr/local/lib/

#DIRINC=/usr/local/core_libraries/include/
#DIRLIB=/usr/local/core_libraries/lib/
#DIRLIB_OACC=./openacc/
#CC=g++

# set flags

CPPFLAGS=-std=c++11 -fPIC -w -g -I$(DIRINC) -I$(DIRINC_BOOST) -I$(DIRLIB_OACC)
#CPPFLAGS=-c11 -v -w -g -I$(DIRINC) -I$(DIRINC_BOOST) -I$(DIRLIB_OACC)
#CCFLAGS=-std=c++11 -c -fPIC -fPIE -w -g -I$(DIRINC) -I$(DIRINC_BOOST) -I$(DIRLIB_OACC)
#LDFLAGS=-std=c++11 -g -w -I$(DIRINC) -I$(DIRINC_BOOST) -I$(DIRLIB_OACC)
LDFLAGS=-std=c++11 -g -w -I$(DIRINC) -I$(DIRINC_BOOST) -I$(DIRLIB_OACC)
#LDFLAGS=-c11 -g -w -I$(DIRINC) -I$(DIRINC_BOOST) -I$(DIRLIB_OACC)

# set sasmol library
#
LIB=sasmol
LIBRARY=$(DIRLIB)lib$(LIB).a

# set openacc library
#
LIB_OACC=oacc_pr
LIBRARY_OACC=$(DIRLIB_OACC)lib$(LIB_OACC).a

# set pgi dependencies
#

DIR_PGILIB_1 = /share/apps/local/pgi/linux86-64/16.10/lib/
## warning the next two may not be on slave nodes
## warning the next two may not be on slave nodes
DIR_PGILIB_2 = /usr/lib64/
DIR_PGILIB_3 = /usr/lib/gcc/x86_64-redhat-linux/4.4.7/

PGI_LIB_1 = $(DIR_PGILIB_1)libaccapi.a
PGI_LIB_2 = $(DIR_PGILIB_1)libaccg.a
PGI_LIB_3 = $(DIR_PGILIB_1)libaccn.a
PGI_LIB_4 = $(DIR_PGILIB_1)libaccg2.a
PGI_LIB_5 = $(DIR_PGILIB_2)libdl.so
PGI_LIB_6 = $(DIR_PGILIB_1)libcudadevice.a

#PGI_LIB_7 = $(DIR_PGILIB)libpg.a
PGI_LIB_8 = $(DIR_PGILIB_2)libmp.a
PGI_LIB_9 = $(DIR_PGILIB_1)libnuma.so

PGI_LIB_10 = $(DIR_PGILIB_2)libpthread.so
PGI_LIB_11 = $(DIR_PGILIB_1)libnspgc.a

PGI_LIB_12 = $(DIR_PGILIB_1)libpgc.a

PGI_LIB_13 = $(DIR_PGILIB_2)libm.a
PGI_LIB_14 = $(DIR_PGILIB_3)libgcc.a

PGI_LIB_15 = $(DIR_PGILIB_2)libc.so

PGI_LIBS = $(PGI_LIB_1) $(PGI_LIB_2) $(PGI_LIB_3) $(PGI_LIB_4) $(PGI_LIB_5) $(PGI_LIB_6) $(PGI_LIB8)
PGI_LIBS += $(PGI_LIB_9) $(PGI_LIB_10) $(PGI_LIB_11) $(PGI_LIB_12) $(PGI_LIB_13) $(PGI_LIB_14) $(PGI_LIB_15) 
#PGI_LIBS += $(PGI_LIB_9) $(PGI_LIB_10) $(PGI_LIB_11) $(PGI_LIB_12) $(PGI_LIB_13) $(PGI_LIB_14) 

# compile pr_mpi

#$(CC) $< $(LDFLAGS) -L$(DIRLIB_BOOST) -L$(DIRLIB) -L$(DIRLIB_OACC) -l$(LIB_OACC)  -l$(LIB) $(DIRLIB_BOOST)/libboost_regex.a $(PGI_LIBS) -o $@

pr_mpi.out: pr_mpi.cpp $(LIBRARY) $(LIBRARY_OACC) 
	$(CC) $< $(CPPFLAGS) -L$(DIRLIB_BOOST) -L$(DIRLIB) -L$(DIRLIB_OACC) -L$(LIBRARY_OACC) -l$(LIB_OACC)  -l$(LIB) $(DIRLIB_BOOST)/libboost_regex.a $(PGI_LIBS) -o $@

clean:
	rm -Rf pr_mpi.out *.o *.so pgprof*

