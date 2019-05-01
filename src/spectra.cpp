/* Written by Leonardo Silvestri 2007-2013 */

//#include <lattice.h>
//#include <N_Functions.h>
//#include <MultiPhonon.h>
//#include <PhononCloud.h>
//#include <OP_model.h>
//#include <OP_epsilon.h>
#include <Transfer_matrix.h>
#include <OPutils.h>
#include <iostream>
#include <vector>
#include <complex>

#ifndef Pi
#define Pi 3.141592654
#endif 
#define MYFLOAT double

namespace ublas = boost::numeric::ublas;
using namespace std;


int main (int argc, char* argv[]) {

//////////////////////////////////////////////////////////
// Reads input file supplied by user

	if (argc<2)
		nrerror("Please supply input file name!");
	
	string input_file=argv[1];
	
	ostringstream ss;
	vector<string> dipole_files;
	string epsilon_file;	
	vector<ublas::vector<complex<MYFLOAT> > > all_dipoles(0);
	vector<MYFLOAT> all_energies(0);
	
	MYFLOAT hw_min;
	MYFLOAT hw_max;
	MYFLOAT hw_step;
	MYFLOAT epsilon_inf; 
	MYFLOAT gamma_coeff;
	MYFLOAT d_slab;
	vector<MYFLOAT> inc_angle;
	vector<MYFLOAT> gamma_angle;
	MYFLOAT na;
	MYFLOAT nf;
	MYFLOAT cell_volume;
	
	read_spectrum_input_data_full(input_file, hw_min, hw_max, hw_step, epsilon_inf, gamma_coeff,
                             d_slab, inc_angle, gamma_angle, na, nf, cell_volume, dipole_files, epsilon_file);
	
	if (dipole_files.size()>0)
	{
		ss.str("");
		ss << input_file.substr(0,input_file.find('.')) << "_epsilon.txt";
		epsilon_file=ss.str();
		write_epsilon(dipole_files, epsilon_file, all_dipoles,
					  all_energies, hw_min, hw_max, hw_step,
					  epsilon_inf, cell_volume, gamma_coeff);
	}
	
	write_absorption_spectrum(hw_min, hw_max, hw_step, d_slab, inc_angle, gamma_angle, na, nf,
							  input_file.substr(0,input_file.find('.')), epsilon_file);
	
	return(0);
}

