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

double **contiguous_2d_array(int rows, int cols){

    double *b =  (double *)malloc(rows * cols * sizeof(double)) ;
    double **array = (double **) malloc(rows * sizeof(double*)) ;
    for (int i = 0 ; i < rows ; i++){
        array[i] = &b[i*cols] ;
    }

    for (int i = 0 ; i < rows ; i++){
        for (int j = 0 ; j < cols ; j++){
            array[i][j] = 0.0 ;
        }
    }

    return array ;
}

double **alloc_2d_double(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    for (int i = 0 ; i < rows ; i++){
        for (int j = 0 ; j < cols ; j++){
            array[i][j] = 0.0 ;
        }
    }

    return array;
}


template <typename T>
T** create2DArray(unsigned nrows, unsigned ncols)
{
   T** ptr = new T*[nrows];  // allocate pointers
   T* pool = new T[nrows*ncols];  // allocate pool
   for (unsigned i = 0; i < nrows; ++i, pool += ncols )
       ptr[i] = pool;
   return ptr;
}

template <typename T>
void delete2DArray(T** arr)
{
   delete [] arr[0];  // remove the pool
   delete [] arr;     // remove the pointers
}


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
    double *localx, *localy, *localz ; 	
    std::cout << "\n\n\n" ;
	const std::string pdb_filename = "ten_mer.pdb" ;
	
    int rank, size ;
    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &size) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
    MPI_Status status ;
	
	util::compile_time(__FILE__, __func__, __DATE__, __TIME__ );

	std::cout << "mpi / multicore / gpu p[r] calculator" << " " << pdb_filename << std::endl ; 
    if(rank == 0){

	    sasmol::SasMol mol ;
	    mol.read_pdb(pdb_filename);

	    std::cout << "number of atoms = " << mol._natoms() << std::endl ;
	    std::cout << "total mass = " << mol._total_mass() << std::endl ;

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
        for(frame=0 ; frame < read_frames ; frame++){ 
            mol.read_dcd_step(dcd_file_pointer, frame, input_natoms, input_reverseEndian) ;	
        } ;	
	    int close_read_result = close_dcd_read(dcd_file_pointer) ;

	    std::cout << "number_of_frames = " << mol.number_of_frames << std::endl ;

	    nframes = mol.number_of_frames ;
	    natoms  = mol.natoms ;

        x = contiguous_2d_array(nframes, natoms) ;
        y = contiguous_2d_array(nframes, natoms) ;
        z = contiguous_2d_array(nframes, natoms) ;

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

        std::cout << "x.rows = " << mol.x.rows() << std::endl ;
        std::cout << "x.cols = " << mol.x.cols() << std::endl ;
    } ;


    const int xn_atoms = natoms ; // number of rows or atoms for x coor
    const int xn_frames = read_frames ; // nframes ; // number of columns or frames for x coor

    localx = new double[xn_frames/xproc_frames * natoms] ;
    localy = new double[xn_frames/xproc_frames * natoms] ;
    localz = new double[xn_frames/xproc_frames * natoms] ;

    int send_value = natoms * xn_frames / xproc_frames ;

    if(rank == 0){
        std::cout << "natoms = " << natoms << std::endl ; 
        std::cout << "send value = " << send_value << std::endl ; 
        MPI_Send(&(x[0][0]), send_value, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD) ;
        MPI_Send(&(y[0][0]), send_value, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD) ;
        MPI_Send(&(z[0][0]), send_value, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD) ;
    }
    else{
        MPI_Recv(&(localx[0]), send_value, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status) ;
        MPI_Recv(&(localy[0]), send_value, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status) ;
        MPI_Recv(&(localz[0]), send_value, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status) ;
    }

    for(p = 1 ; p < size ; p++){
        if(rank == p){
            std::cout << "Local process on rank : " << rank << std::endl ;
            get_distances(localx, localy, localz, xn_frames/xproc_frames, xn_atoms, nbins, bin_width) ; 
       }
    } 

    std::cout << "\ndone with MPI stuff" << std::endl ;  
/*
    delete [] localx ; 
    delete [] localy ; 
    delete [] localz ; 

    delete [] x ; 
    delete [] y ; 
    delete [] z ; 
*/    
    MPI_Finalize() ;
    return 0 ;


}
