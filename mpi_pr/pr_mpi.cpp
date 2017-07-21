#include <sasmol.h>
#include <sasio.h>
#include <sascalc.h>
#include <sasop.h>
#include <sassubset.h>
#include <dcdio.h>
#include <util.h>
#include <vector>
#include <boost/dynamic_bitset.hpp>

#include <vector>
#include "openacc/oacc_pr.h"

#include "mpi.h"

//extern void get_distances();

/********* methods        ******************/

/********* main           ******************/

int main(int argc, char **argv){

    int i, j, p ;
	int frame = 0 ;
//	std::string dcd_input_file_name = "new_ten_mer_1000.dcd" ;
	std::string dcd_input_file_name = "ten_mer.dcd" ;
	int input_natoms, input_nset, input_reverseEndian ;
	int nframes, natoms ;
    int read_frames = 400 ;
    int nbins = 200 ; 
	double bin_width = 1.0 ;  

    int xproc_frames = 2 ;

    double **x, **y, **z ; 	
    std::cout << "\n\n\n" ;
	const std::string pdb_filename = "ten_mer.pdb" ;
	
    int rank, size ;
    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &size) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
	
	util::compile_time(__FILE__, __func__, __DATE__, __TIME__ );

	std::cout << "mpi / multicore / gpu p[r] calculator" << " " << pdb_filename << std::endl ; 
    if(rank == 0){

	    sasmol::SasMol mol ;
	    mol.read_pdb(pdb_filename);

	    std::cout << "number of atoms = " << mol._natoms() << std::endl ;
	    std::cout << "total mass = " << mol._total_mass() << std::endl ;

	    std::cout << "number_of_frames = " << mol.number_of_frames << std::endl ;

	    FILE *dcd_file_pointer ;

	    dcd_file_pointer = sasio::open_dcd_read(dcd_input_file_name, input_natoms, input_nset, input_reverseEndian) ;

	    std::cout << "input_natoms = " << input_natoms << std::endl ;
	    std::cout << "input_nset = " << input_nset << std::endl ;
	    std::cout << "input_reverseEndian = " << input_reverseEndian << std::endl ;

        input_nset = read_frames ; 

        mol.number_of_frames = 0 ;
        mol.x.setZero(mol.natoms, input_nset);
        mol.y.setZero(mol.natoms, input_nset);
        mol.z.setZero(mol.natoms, input_nset);

//        mol.read_dcd(dcd_input_file_name) ;

//        for(frame=0 ; frame < input_nset ; frame++){ 
        for(frame=0 ; frame < 200 ; frame++){ 
            mol.read_dcd_step(dcd_file_pointer, frame, input_natoms, input_reverseEndian) ;	
        } ;	
	    int close_read_result = close_dcd_read(dcd_file_pointer) ;

	    std::cout << "number_of_frames = " << mol.number_of_frames << std::endl ;

	    nframes = mol.number_of_frames ;
	    natoms  = mol.natoms ;

        x = new double*[nframes] ;
        y = new double*[nframes] ;
        z = new double*[nframes] ;

        for(int i=0 ; i < nframes ; i++){
            x[i] = new double[natoms] ;
            y[i] = new double[natoms] ;
            z[i] = new double[natoms] ;
        }

        std::cout << "hello here is some eigen stuff" << "\n" ;
        std::cout << "atom 0: frame 0" << std::endl ;
        std::cout << "x(0,0) = " << mol.x(0,0) << "\n" ;
        std::cout << "y(0,0) = " << mol.y(0,0) << "\n" ;
        std::cout << "z(0,0) = " << mol.z(0,0) << "\n" ;
        std::cout << "atom 1: frame 0" << std::endl ;
        std::cout << "x(1,0) = " << mol.x(1,0) << "\n" ;
        std::cout << "y(1,0) = " << mol.y(1,0) << "\n" ;
        std::cout << "z(1,0) = " << mol.z(1,0) << "\n" ;
        std::cout << "atom 0: frame 1" << std::endl ;
        std::cout << "x(0,1) = " << mol.x(0,1) << "\n" ;
        std::cout << "y(0,1) = " << mol.y(0,1) << "\n" ;
        std::cout << "z(0,1) = " << mol.z(0,1) << "\n" ;
        std::cout << "atom 1: frame 1" << std::endl ;
        std::cout << "x(1,1) = " << mol.x(1,1) << "\n" ;
        std::cout << "y(1,1) = " << mol.y(1,1) << "\n" ;
        std::cout << "z(1,1) = " << mol.z(1,1) << "\n" ;

        std::cout<< "size x = " << mol.x.size() << "\n" ;

        for(frame=0 ; frame < mol.number_of_frames ; frame++){
            for(int i=0 ; i < mol.natoms; i++){
                x[frame][i] = mol.x(i, frame) ; 
                y[frame][i] = mol.y(i, frame) ; 
                z[frame][i] = mol.z(i, frame) ; 

            } 
        }

    } ;

    const int xn_atoms = 148 ; // natoms ; // number of rows or atoms for x coor
    //const int xn_frames = 1000 ; // nframes ; // number of columns or frames for x coor
    const int xn_frames = read_frames ; // nframes ; // number of columns or frames for x coor

    double **localx = new double*[xn_frames/xproc_frames] ;
    double **localy = new double*[xn_frames/xproc_frames] ;
    double **localz = new double*[xn_frames/xproc_frames] ;

    for(int i=0 ; i < nframes ; i++){
        localx[i] = new double[natoms] ;
        localy[i] = new double[natoms] ;
        localz[i] = new double[natoms] ;
    }
        
    int sizes[2] = {xn_atoms, xn_frames} ;
    int subsizes[2] = {xn_atoms, xn_frames/xproc_frames} ;
    int starts[2] = {0, 0} ;

    MPI_Datatype type, subarrtype ;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type) ;
    MPI_Type_create_resized(type, 0, xn_frames/xproc_frames *sizeof(double), &subarrtype) ;
    MPI_Type_commit(&subarrtype) ;

    int sendcounts[xn_atoms*xproc_frames] ;
    int displs[xn_atoms*xproc_frames] ;
    double *globalptrx, *globalptry, *globalptrz ;
    
    if(rank == 0) {
        globalptrx = &(x[0][0]) ;
        globalptry = &(y[0][0]) ;
        globalptrz = &(z[0][0]) ;
    
    }

    if (rank == 0) {
        for (i=0; i<xn_atoms*xproc_frames; i++)
            sendcounts[i] = 1;
        int disp = 0;
        for (i=0; i<xn_atoms; i++) {
            for (j=0; j<xproc_frames; j++) {
                displs[i*xproc_frames+j] = disp;
                disp += 1;
            }
            disp += ((xn_frames/xproc_frames)-1)*xproc_frames;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD) ;
/*
    MPI_Scatterv(globalptrx, sendcounts, displs, subarrtype, &(localx[0][0]),
                 xn_atoms*xn_frames/(xproc_frames), MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD) ;

    MPI_Scatterv(globalptry, sendcounts, displs, subarrtype, &(localy[0][0]),
        xn_atoms*xn_frames/(xproc_frames), MPI_DOUBLE,
        0, MPI_COMM_WORLD);

    MPI_Scatterv(globalptrz, sendcounts, displs, subarrtype, &(localz[0][0]),
        xn_atoms*xn_frames/(xproc_frames), MPI_DOUBLE,
        0, MPI_COMM_WORLD);
*/

//    std::vector<std::vector<int> > local_hist(xn_frames/xproc_frames, std::vector<int>(nbins, 0));
 //   get_distances(localx, localy, localz, xn_frames/xproc_frames, natoms, local_hist, nbins, bin_width) ; 
 
    //get_distances(localx, localy, localz, xn_frames/xproc_frames, xn_atoms, nbins, bin_width) ; 
   // for(p = 0 ; p < size ; p++){
    //    if(rank == p){
     //       std::cout << "Local process on rank : " << rank << std::endl ;
      //      //get_distances(localx, localy, localz, xn_frames/xproc_frames, xn_atoms, nbins, bin_width) ; 
       // }
        //MPI_Barrier(MPI_COMM_WORLD) ;
    //}

    std::cout << "done with MPI stuff" << std::endl ;  
    MPI_Barrier(MPI_COMM_WORLD) ;

//    delete [] localx ; 
 //   delete [] localy ; 
  //  delete [] localz ; 

    MPI_Type_free(&subarrtype) ;
    
    //delete [] x ; 
   // delete [] y ; 
   // delete [] z ; 
    
    MPI_Finalize() ;
    return 0 ;


}
