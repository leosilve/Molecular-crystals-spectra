/* Written by Leonardo Silvestri 2007-2013 */

#include <lattice.h>
//#include <N_Functions.h>
//#include <MultiPhonon.h>
#include <PhononCloud.h>
//#include <OP_model.h>
//#include <OP_epsilon.h>
//#include <OP_absorption.h>
#include <OPutils.h>
#include <iostream>
#include <vector>
#include <complex>

#ifndef Pi
#define Pi 3.141592654
#endif 
#define MYFLOAT float

namespace ublas = boost::numeric::ublas;
using namespace std;

typedef ublas::matrix<MYFLOAT, ublas::column_major>	Matrix;
typedef ublas::vector<MYFLOAT>	Vector;

int main () {

//////////////////////////////////////////////////////////
// Reads input file OP_input.txt

	string OP_input_filename="OP_input.txt";
	ifstream infile(OP_input_filename.c_str());

	if (!infile.is_open()) nrerror("OP input file not found!");

	string lattice_file;
	string Bands_input_file;
	string Basis_file;
	string OPmodel_file;
	string emission_input_file;
	string absorption_input_file;
	string Problem_name;
	int FLAG_READ;

	get_next_string_from_file(infile, lattice_file);
	get_next_string_from_file(infile, Bands_input_file);
	get_next_string_from_file(infile, Basis_file);
	get_next_string_from_file(infile, OPmodel_file);
	get_next_string_from_file(infile, emission_input_file);
	get_next_string_from_file(infile, absorption_input_file);
	get_next_string_from_file(infile, Problem_name);
	get_next_from_file(infile, FLAG_READ);
	infile.close();

//////////////////////////////////////////////////////////
// Constructs lattice and computes excitonic bands

	lattice<MYFLOAT>	MyLattice(lattice_file);

//	if (Bands_input_file!="")
//		MyLattice.computes_lattice_bands(Bands_input_file);
	
	PhononCloud<MYFLOAT>	MyBasis(Basis_file, &MyLattice);
	ostringstream ss;
	ss.str("");
	ss << Problem_name << ".log";
	ofstream log_file(ss.str().c_str());
	MyLattice.print(log_file);
	MyBasis.print(log_file);
//	MyModel.print(log_file);
	log_file.close();
	
	return 0;
/*
	
//////////////////////////////////////////////////////////
// Constructs and solves Model problem

// Controllare la BD_BASIS e tutta la procedura di soluzione per la base MultiPhonon
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

	if (Basis_file!="" && OPmodel_file!="")
	{

	PhononCloud<MYFLOAT>	MyBasis(Basis_file, &MyLattice);
	OPmodel<MYFLOAT, PhononCloud<MYFLOAT> >		MyModel(OPmodel_file, &MyBasis);
//	MultiPhonon<MYFLOAT>	MyBasis(Problem_name, Basis_file, &MyLattice);
//	OPmodel<MYFLOAT, MultiPhonon<MYFLOAT> >		MyModel(OPmodel_file, &MyBasis);

	ostringstream ss;
	ss.str("");
	ss << Problem_name << ".log";
	ofstream log_file(ss.str().c_str());
	MyLattice.print(log_file);
	MyBasis.print(log_file);
	MyModel.print(log_file);
	log_file.close();

	string eigvec_filename, eigval_filename, dipoles_filename;
	ss.str("");
	ss << Problem_name << "_eigenvectors.txt";
	eigvec_filename=ss.str();
	ss.str("");
	ss << Problem_name << "_dipoles.txt";
	dipoles_filename=ss.str();
	ss.str("");
	ss << Problem_name << "_eigenvalues.txt";
	eigval_filename=ss.str();

	if (FLAG_READ!=1)
	{
		MyModel.solve();	
		MyModel.print_results(eigvec_filename, eigval_filename, dipoles_filename);
	}
	else
	{
		if (!MyModel.read_eigenstates(eigvec_filename))
			cout << "Error: could not read eigenstates!" << endl;
		MyModel.print_dipoles(dipoles_filename);
	}
//	return 0;
/////////////////////////////////////////////////////////
// Computes emission

	if (emission_input_file!="")
	{
		string emission_output_file;
		ss.str("");
		ss << Problem_name << "_emission.txt";
		emission_output_file=ss.str();
		MyModel.compute_DIMER_emission(emission_input_file, emission_output_file);
	}
//	return 0;

/////////////////////////////////////////////////////////
// Writes epsilon

	vector<string> dipole_files;
	dipole_files=MyModel.get_dipole_files();
	if (!dipole_files.size()>0)
	{
		dipole_files.resize(1);
		dipole_files[0]=dipoles_filename;
	}
	
	string epsilon_file;
	ss.str("");
	ss << Problem_name << "_epsilon.txt";
	epsilon_file=ss.str();

	vector<ublas::vector<complex<MYFLOAT> > > all_dipoles(0);
	vector<MYFLOAT> all_energies(0);
	MYFLOAT hw_step=MyModel.get_hw_step();
	MYFLOAT epsilon_inf=MyModel.get_eps_inf(); 
	MYFLOAT cell_volume=MyLattice.get_cell_volume();
	MYFLOAT	gamma_coeff=MyModel.get_gamma_coeff();
	write_epsilon(	dipole_files, epsilon_file, all_dipoles,
					all_energies, hw_step,
					epsilon_inf, cell_volume, gamma_coeff);

//////////////////////////////////////////////////////////
// Computes absorption spectrum

	string spectrum_output_file;
	ss.str("");
	ss << Problem_name << "_absorption_spectrum.txt";
	spectrum_output_file=ss.str();
	write_absorption_spectrum(absorption_input_file, spectrum_output_file, epsilon_file);
	return(0);
	
	}  // end if model and basis files exist
*/
}

