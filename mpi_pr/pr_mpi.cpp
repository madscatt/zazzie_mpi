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

//extern "C" {
#include "openacc/oacc_pr.h"
//}


//extern void get_distances();

/********* methods        ******************/

/********* main           ******************/

int main(){

	std::cout << "\n\n\n" ;
	
	sasmol::SasMol mol ;
	
	util::compile_time(__FILE__, __func__, __DATE__, __TIME__ );

//	const std::string filename="ubq.xyz" ;
//	mol.read_xyz(filename) ;
      
//	const std::string pdb_filename = "min3.pdb" ;
	const std::string pdb_filename = "ten_mer.pdb" ;
	const std::string output_filename = "new_min3.pdb" ;
	const std::string dcd_output_filename="new_min3.dcd" ;

	std::cout << "testing sasmol " << " " << pdb_filename << std::endl ; 

	mol.read_pdb(pdb_filename);

	int frame = 0 ;

	std::cout << "number of atoms = " << mol._natoms() << std::endl ;
	std::cout << "total mass = " << mol._total_mass() << std::endl ;

	std::vector<float> com = mol.calc_com(frame) ;
	util::print(com) ;	

	std::cout << "number_of_frames = " << mol.number_of_frames << std::endl ;

	std::cout << "TESTING DCD READ" << std::endl << std::endl ;

	std::string dcd_input_file_name = "ten_mer.dcd" ;
	int input_natoms, input_nset, input_reverseEndian ;

	FILE *dcd_file_pointer ;

	dcd_file_pointer = sasio::open_dcd_read(dcd_input_file_name, input_natoms, input_nset, input_reverseEndian) ;

	std::cout << "BACK FROM OPEN DCD READ" << std::endl ;

	std::cout << "input_natoms = " << input_natoms << std::endl ;
	std::cout << "input_nset = " << input_nset << std::endl ;
	std::cout << "input_reverseEndian = " << input_reverseEndian << std::endl ;

	std::cout << "TESTING DCD READ STEP" << std::endl << std::endl ;
	int input_frame = 0 ;

    mol.read_dcd(dcd_input_file_name) ;

/*    int f = 0 ;
    for(frame=0 ; frame < input_nset ; frame++){ 
//	    std::cout << "frame = " << frame << std::endl ; 
        mol.read_dcd_step( dcd_file_pointer, f , input_natoms, input_reverseEndian) ;	
    } ;	
	int close_read_result = close_dcd_read(dcd_file_pointer) ;
*/

	std::cout << "number_of_frames = " << mol.number_of_frames << std::endl ;

	int nframes, natoms, nbins ;
	double bin_width ;  

	nframes = mol.number_of_frames ;
	natoms  = mol.natoms ;
	nbins  = 200 ;
	bin_width = 1.0 ;

    std::vector<std::vector<int> > hist(nframes, std::vector<int>(nbins,       0));

    double **x = new double*[nframes] ;
    double **y = new double*[nframes] ;
    double **z = new double*[nframes] ;

    for(int i=0 ; i < nframes ; i++){
        x[i] = new double[natoms] ;
        y[i] = new double[natoms] ;
        z[i] = new double[natoms] ;
    }
    

/*    malloc2D(&x, nframes, natoms);
    malloc2D(&y, nframes, natoms);
    malloc2D(&z, nframes, natoms);
 */ 
    std::cout << "hello here is some eigen stuff" << "\n" ;
    std::cout << "x(0,0) = " << mol.x(0,0) << "\n" ;
    std::cout << "y(0,0) = " << mol.y(0,0) << "\n" ;
    std::cout << "z(0,0) = " << mol.z(0,0) << "\n" ;

    std::cout<< "size x = " << mol.x.size() << "\n" ;

   // for(frame=0 ; frame < 2 ; frame++){
    for(frame=0 ; frame < mol.number_of_frames ; frame++){
        for(int i=0 ; i < mol.natoms; i++){
            x[frame][i] = mol.x(i, frame) ; 
            y[frame][i] = mol.y(i, frame) ; 
            z[frame][i] = mol.z(i, frame) ; 

        } 
    }

    get_distances(x, y, z, nframes, natoms, hist, nbins, bin_width) ; 
    //get_distances(c_array, nframes, natoms, hist, nbins, bin_width) ; 
  
/* 
    free2D(&x) ; 
    free2D(&y) ; 
    free2D(&z) ;
*/ 
    return 0 ;


}
