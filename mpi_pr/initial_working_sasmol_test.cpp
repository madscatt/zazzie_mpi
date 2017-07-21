#include <sasmol.h>
#include <sasio.h>
#include <sascalc.h>
#include <sasop.h>
#include <sassubset.h>
#include <dcdio.h>
#include <util.h>
#include <vector>
#include <boost/dynamic_bitset.hpp>

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

	float rg = mol.calc_rg(frame) ;
	std::cout << "radius of gyration = " << rg << std::endl ;
	std::cout << "sasmol.py : radius of gyration = 64.043168998442439 " << std::endl ;

	mol.center(frame) ;
	
	mol.move_to(frame,com) ;

	mol.calc_pmi(frame) ;

	std::cout << "mol.uk = " << mol.uk << std::endl ;
	std::cout << "mol.ak = " << mol.ak << std::endl ;

	float cutoff = 0.8 ;
	int check = 1 ;  // 1 == overlap and 0 == no overlap

	std::cout << ">>> checking for overlap (xyz): cutoff = "<< cutoff << std::endl ;

	clock_t tStart = clock();

	check = mol.self_overlap(cutoff,frame) ;

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	if(check == 0) 
	{
		std::cout << ">>> NO OVERLAP FOUND <<<" << std::endl;
	}
	else
	{
		std::cout << ">>> OVERLAP FOUND <<<" << std::endl;
	}


	std::cout << "number_of_frames = " << mol.number_of_frames << std::endl ;

	mol.write_dcd(dcd_output_filename) ;	

	mol.write_pdb(output_filename,frame);

	const std::string second_dcd_output_filename="second_new_min3.dcd" ;
	
	FILE *outfile = sasio::open_write_dcd_file(second_dcd_output_filename, mol.natoms, mol.number_of_frames);

	std::vector<float> value = {0.0,0.0,0.0} ;
	float vcount = 0.0 ;
	std::cout << value[0] << " " << value[1] << " " << value[2] << std::endl ;
	
	tStart = clock();
	int step = 0 ;

	const std::string translated_pdb_filename = "t_min3.pdb" ;

	for(int j = 0 ; j < 1000 ; ++j)
	{
		mol.translate(frame,value) ;	

		mol.write_dcd_step(outfile,step,j);

		value[0] += 0.0001 ;
		value[1] += 0.0001 ;
		value[2] += 0.0001 ;
	
		vcount += value[0] ;
		if(j == 999) mol.write_pdb(translated_pdb_filename,frame) ;
	}

	int close_result = close_dcd_write(outfile) ;
	
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	std::cout << vcount << std::endl;

	std::cout << value[0] << " " << value[1] << " " << value[2] << std::endl ;

	std::cout << "TESTING DCD READ" << std::endl << std::endl ;

	std::string dcd_input_file_name = "new_min3.dcd" ;
	int input_natoms, input_nset, input_reverseEndian ;

	FILE *dcd_file_pointer ;

	dcd_file_pointer = sasio::open_dcd_read(dcd_input_file_name, input_natoms, input_nset, input_reverseEndian) ;

	std::cout << "BACK FROM OPEN DCD READ" << std::endl ;

	std::cout << "input_natoms = " << input_natoms << std::endl ;
	std::cout << "input_nset = " << input_nset << std::endl ;
	std::cout << "input_reverseEndian = " << input_reverseEndian << std::endl ;

	std::cout << "TESTING DCD READ STEP" << std::endl << std::endl ;

	sasmol::SasMol mol2 ; 
	mol2.read_pdb(pdb_filename);
	int input_frame = 0 ;

	mol2.read_dcd_step( dcd_file_pointer, input_frame, input_natoms, input_reverseEndian) ;	
	
	int close_read_result = close_dcd_read(dcd_file_pointer) ;

	std::string second_pdb_filename = "mol2_new_min3.pdb" ;

	mol2.write_pdb(second_pdb_filename,frame);

	std::cout << "BACK FROM OPEN DCD READ" << std::endl ;

	std::cout << "TESTING DCD READ (FULL)" << std::endl << std::endl ;

	std::string second_dcd_filename = "second_new_min3.dcd" ;
	
	mol.read_dcd(second_dcd_filename);
	
	std::cout << "BACK FROM DCD READ (FULL)" << std::endl ;

	std::cout << "number_of_frames = " << mol.number_of_frames << std::endl ;

	return 0 ;


}
