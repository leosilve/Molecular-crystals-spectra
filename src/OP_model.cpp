/* Written by Leonardo Silvestri 2007-2013 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <math.h>
#include <vector>
#include <map>
#include <OPutils.h>
#include <N_Functions.h> 
#include <physics.h>
#include <lattice.h>
#include <QuantumState.h>
#include <OP_model.h>
#include <PhononCloud.h>
#include <Transfer_matrix.h>
/*#include "arcomp.h" 
#include "arlnsmat.h" // ARluNonSymMatrix definition.
#include "arlscomp.h" // ARluCompStdEig definition.
#include "arlssym.h " // ARluSymStdEig definition.
#include "arlsmat.h " // ARluSymMatrix definition.*/
//#include <chrono>
#include <ctime>

#ifndef IMAG_TOL
#define IMAG_TOL 1e-5
#endif 
#ifndef Pi
#define Pi 3.141592654
#endif 
#ifndef ZERO_TOL
#define ZERO_TOL 1e-10
#endif 

using namespace std;

template <class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::initialize(string input_file)
{
	cout << "Reading model data from file: " << input_file << endl;
	
	model_name=input_file.substr(0,input_file.find('.'));
	string line;
	istringstream iss;
	ifstream infile(input_file.c_str());
	if (!infile.is_open()) nrerror("File not found!");
	
	std::string			lattice_file;
	std::string			Basis_file;
	get_next_string_from_file(infile, model_name);
	get_next_string_from_file(infile, calc);
	get_next_string_from_file(infile, lattice_file);
	get_next_string_from_file(infile, Basis_file);
	get_next_string_from_file(infile, emission_input_file);
	get_next_string_from_file(infile, absorption_input_file);
	get_next_from_file(infile, FLAG_REAL);
	get_next_from_file(infile, NEV);
	getline(infile, line, ':');
	getline(infile, line);
	iss.str( line );
	copy( istream_iterator<int>( iss ), istream_iterator<int>(), back_inserter( DCEV ) );
	infile.close();
	
	lattice<OPFLOAT>		MyLattice(lattice_file);
	basisPtr = new basis(Basis_file, lattice_file);
	
	ostringstream ss;
	ss.str("");
	ss << model_name << "_log.txt";
	ofstream log_file(ss.str().c_str());
	print(log_file);
	log_file.close();
	
	
	int i,j,n;
	if (norm_2(basisPtr->get_KVEC())>ZERO_TOL)
		basisPtr->change_KVEC(0.0,0.0,0.0);
	if (DCEV.size()>0 && calc[0]!='e' && calc[0]!='E')
	{
		basisPtr->compute_BD_BASIS();
		if (basisPtr->get_BD_BASIS_size()==0)
			nrerror("OPmodel<OPFLOAT,basis>::initialize: BD_BASIS not found!");
		all_eigenstates=0;
		DC_eigenstates.resize(DCEV.size());
		for (i=0; i<DCEV.size(); i++)
		{
			if (DCEV[i]>=basisPtr->get_BD_BASIS_size())
				nrerror("OPmodel<OPFLOAT,basis>::initialize: BD_BASIS[i] not found! Try with a different Davydov component (PC basis can be decomposed into only 2 subspaces!)");
			n=basisPtr->get_BD_BASIS_size(DCEV[i]);
			DC_eigenstates[i] = new std::vector<QuantumState<OPFLOAT> >(n, sizeof(QuantumState<OPFLOAT>(n)));
		}
	}
	else
	{
		n=(*basisPtr).get_BASIS_size();
		all_eigenstates = new std::vector<QuantumState<OPFLOAT> >(n, sizeof(QuantumState<OPFLOAT>(n)));
		if (DCEV.size()>0)
		{
			cout << "Cannot use Block Diagonal decomposition ... setting DCEV.size()=0!" << endl;
			DCEV.resize(0);
		}
		DC_eigenstates.resize(0);
	}	
	return;
}

template<class OPFLOAT, class basis> 
bool OPmodel<OPFLOAT,basis>::print_eigenvectors(string eigvec_filename, vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr)
{
	if ((*eigenstates_Ptr).size()<=0)	
	{ 
		cout << "OPmodel::print_eigenvectors: no eigenvectors " << endl;
		return(false);
	}
	
	ofstream vec_file(eigvec_filename.c_str());
	int i, j;
	int n = (*eigenstates_Ptr).size();
	int dim = (*eigenstates_Ptr)[0].size();
	vec_file << n << "\t" << dim << endl;
	for (i=0; i<(*eigenstates_Ptr).size(); i++)
	{
		if ((*eigenstates_Ptr)[i].size()!=dim)
			nrerror("OP_model::print_eigenvectors: not all eigenvectors have the same size!");
		if (FLAG_REAL==1)
			for (j=0; j<dim; j++)
					vec_file << /*std::fixed << std::setprecision(4) <<*/ (*eigenstates_Ptr)[i].get_coeff(j).real() << " ";
		else
			for (j=0; j<dim; j++)
				vec_file << /*std::fixed << std::setprecision(4) <<*/ (*eigenstates_Ptr)[i].get_coeff(j).real() << " "
				         << /*std::fixed << std::setprecision(4) <<*/ (*eigenstates_Ptr)[i].get_coeff(j).imag() << " ";
		vec_file << endl;
	}
	vec_file.close();
	return(true);
}
	
template<class OPFLOAT, class basis> 
bool OPmodel<OPFLOAT,basis>::print_eigenvalues(string eigval_filename, vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr)
{
	if ((*eigenstates_Ptr).size()<=0)	
	{
		cout << "OPmodel::print_eigenvalues: no eigenvectors " << endl;
		return(false);
	}
	ofstream val_file(eigval_filename.c_str());
	int i, j;
	int n = eigenstates_Ptr->size();
	val_file << n << endl;
	for (i=0; i<n; i++)
			val_file <<  /*") E: " << std::fixed << std::setprecision(5) <<*/ (*eigenstates_Ptr)[i].get_energy() << endl;
	val_file.close();
	return(true);
}

template<class OPFLOAT, class basis> 
bool OPmodel<OPFLOAT,basis>::print_dipoles(string dipoles_filename, vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr, int dc)
{
	if ((*eigenstates_Ptr).size()<=0)	
	{
		cout << "OPmodel::print_dipoles: no eigenvectors " << endl;
		return(false);
	}
	int i,j;
	int n = (*eigenstates_Ptr).size();
	int dim = (*eigenstates_Ptr)[0].size();
	ofstream dipole_file(dipoles_filename.c_str());
	ComplexVector temp_dipole(3);	
	OPFLOAT total_dipole=0.0;
	for (i=0; i<n; i++)
	{
		if ((*eigenstates_Ptr)[i].size()!=dim)
			nrerror("OP_model::print_eigenvectors: not all eigenvectors have the same size!");
		temp_dipole[0]=0.0;temp_dipole[1]=0.0;temp_dipole[2]=0.0;	
		if (dc<0)
		{
			for (j=0; j<dim; j++)
				temp_dipole+=(*eigenstates_Ptr)[i].get_coeff(j)*(*basisPtr).get_BASIS_dipole(j);
		}
		else
		{
			for (j=0; j<dim; j++)
				temp_dipole+=(*eigenstates_Ptr)[i].get_coeff(j)*(*basisPtr).get_BD_BASIS_dipole(dc,j);
		}	
		if (norm_2(temp_dipole)>ZERO_TOL )
			dipole_file << endl << std::fixed << std::setprecision(5) << (*eigenstates_Ptr)[i].get_energy()	<< "\t" 
									<< std::fixed << std::setprecision(4) 
				            << temp_dipole[0].real() << "\t" << temp_dipole[0].imag() << "\t" 
							<< temp_dipole[1].real() << "\t" << temp_dipole[1].imag() << "\t" 
							<< temp_dipole[2].real() << "\t" << temp_dipole[2].imag(); 										
		total_dipole+=pow(norm_2(temp_dipole),2);
	}
	dipole_file.close();
	cout << "Basis set total dipole: " << total_dipole << endl;
	total_dipole=0.0;
	for (i=0; i<(*(*basisPtr).get_latticePtr()).get_N_MOL(); i++)
		total_dipole+=pow(norm_2((*(*basisPtr).get_latticePtr()).get_MOL_DIPOLE(i)),2);
	cout << "Expected total dipole: " << total_dipole << endl;
	return(true);
}


template<class OPFLOAT, class basis> 
bool OPmodel<OPFLOAT,basis>::print_eigenvectors()
{
	vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr;
	ostringstream ss;
	int dc;
		
	if (DCEV.size()==0)
	{
		if (all_eigenstates!=0)
		{
			eigenstates_Ptr=all_eigenstates;
			ss.str("");
			ss << model_name << "_eigvec_all.txt";
			print_eigenvectors(ss.str(), eigenstates_Ptr);
			return(true);
		}
		else
			return(false);
	}

	for (dc=0; dc<DCEV.size(); dc++)
	{
			if (DC_eigenstates[dc]!=0)
			{
				eigenstates_Ptr=DC_eigenstates[dc];
				ss.str("");
				ss << model_name << "_eigvec_" << DCEV[dc] << ".txt";
				print_eigenvectors(ss.str(), eigenstates_Ptr);
			}
			else
				return(false);
	}
	return(true);
}

template<class OPFLOAT, class basis> 
bool OPmodel<OPFLOAT,basis>::print_eigenvalues()
{
	vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr;
	ostringstream ss;
	int dc;
		
	if (DCEV.size()==0)
	{
		if (all_eigenstates!=0)
		{
			eigenstates_Ptr=all_eigenstates;
			ss.str("");
			ss << model_name << "_eigval_all.txt";
			print_eigenvalues(ss.str(), eigenstates_Ptr);
			return(true);
		}
		else
			return(false);
	}

	for (dc=0; dc<DCEV.size(); dc++)
	{
			if (DC_eigenstates[dc]!=0)
			{
				eigenstates_Ptr=DC_eigenstates[dc];
				ss.str("");
				ss << model_name << "_eigval_" << DCEV[dc] << ".txt";
				print_eigenvalues(ss.str(), eigenstates_Ptr);
			}
			else
				return(false);
				
	}
	return(true);	
}

template<class OPFLOAT, class basis> 
vector<string> OPmodel<OPFLOAT,basis>::print_dipoles()
{
	vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr;
	ostringstream ss;
	int dc;
	vector<string> dipole_files(0);
		
	if (DCEV.size()==0)
	{
		if (all_eigenstates!=0)
		{
			eigenstates_Ptr=all_eigenstates;
			ss.str("");
			ss << model_name << "_dipoles_all.txt";
			print_dipoles(ss.str(), eigenstates_Ptr,-1);
			dipole_files.push_back(ss.str());
		}
		return(dipole_files);
	}

	for (dc=0; dc<DCEV.size(); dc++)
	{
		if (DC_eigenstates[dc]!=0)
		{
				eigenstates_Ptr=DC_eigenstates[dc];
				ss.str("");
				ss << model_name << "_dipoles_" << DCEV[dc] << ".txt";
				print_dipoles(ss.str(), eigenstates_Ptr, dc);
				dipole_files.push_back(ss.str());
		}
	}
	return(dipole_files);
}
	

template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::solve_complex_full(int dc, int flag_macro)
{
	int n;	// Dimension of the problem.
	int i,j;
	vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr;
	if (dc<0)
	{
		n = (*basisPtr).get_BASIS_size();
		eigenstates_Ptr=all_eigenstates;
	}
	else
	{
		if (dc<DC_eigenstates.size())
		{
			n = basisPtr->get_BD_BASIS_size(DCEV[dc]);
			eigenstates_Ptr=DC_eigenstates[dc];
		}
		else
			nrerror("OP_model::solve_BASIS: Davydov component not found!");
	}

// Solve with heev from LAPACK
		cout << "Solving the full BASIS set with LAPACK heev:  dc = " << dc << ", n = " << n << " ..." << endl; 
		char jobz = 'V';	// = 'N':  Compute eigenvalues only;
							// = 'V':  Compute eigenvalues and eigenvectors.
		char uplo = 'U';	// = 'U':  Upper triangle of A is stored;
							// = 'L':  Lower triangle of A is stored.
		vector<OPFLOAT>					eigenvalues(n);
		vector<Complex>					eigenvectors(n*n,0.0); 
		for (i=0; i<n; i++)
			for (j=i; j<n; j++)
			{
//					cout << i << " " << j << " " <<  (*basisPtr).Hint(i,j,dc,flag_macro) << endl;
					eigenvectors[j*n+i]=(*basisPtr).Hint(i,j,dc,flag_macro);
//					if (i==0 && j==0) cout << i << "\t" << j << "\t" << eigenvectors[j][i] << endl;
			}
		heev_cpp(jobz, uplo, n, eigenvectors, eigenvalues, 'O');
		for (i=0; i<n; i++)
		{
			(*eigenstates_Ptr)[i].set_energy(eigenvalues[i]);
			(*eigenstates_Ptr)[i].resize(n);
			for (j=0; j<n; j++)
				(*eigenstates_Ptr)[i].set_coeff(j,eigenvectors[i*n+j]);
			(*eigenstates_Ptr)[i].normalize();
		}
		cout << " ... done" << endl;
	return;
}


template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::solve_real_full(int dc, int flag_macro)
{
	int n;	// Dimension of the problem.
	int i,j;
	vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr;
	if (dc<0)
	{
		n = (*basisPtr).get_BASIS_size();	// Dimension of the problem.
		eigenstates_Ptr=all_eigenstates;
	}
	else
	{
		if (dc<DC_eigenstates.size())
		{
			n = (*basisPtr).get_BD_BASIS_size(dc);	// Dimension of the problem.
			eigenstates_Ptr=DC_eigenstates[dc];
		}
		else
			nrerror("OP_model::solve_BASIS: Davydov component not found!");
	}

// Solve with syev from LAPACK
	cout << "Finding all the eigenvalues of the real H with LAPACK syev: dc=" << dc << ", n=" << n << endl;
		char jobz = 'V';	// = 'N':  Compute eigenvalues only;
							// = 'V':  Compute eigenvalues and eigenvectors.
		char uplo = 'U';	// = 'U':  Upper triangle of A is stored;
							// = 'L':  Lower triangle of A is stored.
		vector<OPFLOAT>		eigenvalues(n);
		vector<OPFLOAT>		eigenvectors(n*n,0);
		Complex comp_val;

		for (i=0; i<n; i++)
			for (j=i; j<n; j++)
			{
				comp_val=(*basisPtr).Hint(i,j,dc,flag_macro);
				if (abs(comp_val.imag()/comp_val.real())>IMAG_TOL)
				{	
					cout << "comp_val = " << comp_val << endl; 
					nrerror("solve_real: H matrix element is not real!");
				}
				eigenvectors[j*n+i]=comp_val.real();
			}

		syev_cpp(jobz, uplo, n, eigenvectors, eigenvalues,'O');
		(*eigenstates_Ptr).resize(n);
		for (i=0; i<n; i++)
		{
			(*eigenstates_Ptr)[i].set_energy(eigenvalues[i]);
			(*eigenstates_Ptr)[i].resize(n);
			for (j=0; j<n; j++)
				(*eigenstates_Ptr)[i].set_coeff(j,eigenvectors[i*n+j]);
			(*eigenstates_Ptr)[i].normalize();
		}
	return;
}

template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::solve(int flag_macro)
{
	if (DCEV.size()>0)
	{
//		cout << "Check BD decomposition: " << (*basisPtr).check_BD_decomposition(0,1,1e-5) << endl;
		if (FLAG_REAL==1)
		{
			for (int i=0; i<DCEV.size(); i++)
			{
				if (NEV>0)
					solve_real_nev(i, flag_macro);
				else
					solve_real_full(i, flag_macro);
			}
		}
		else
		{ 
			for (int i=0; i<DCEV.size(); i++)
			{
				if (NEV>0)
					solve_complex_nev(i, flag_macro);
				else
					solve_complex_full(i, flag_macro);
			}	
		}
	}
	else
	{
		if (FLAG_REAL==1)
		{
			if (NEV>0)
				solve_real_nev(-1, flag_macro);
			else
				solve_real_full(-1, flag_macro);
		}
		else
		{
			if (NEV>0)
				solve_complex_nev(-1, flag_macro);
			else
				solve_complex_full(-1, flag_macro);
		}
	}
	return;
}

template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::solve()
{
	if (calc[0]=='e' || calc[0]=='E')
		compute_emission();
	else
	    compute_absorption();
	return;
}

template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::compute_absorption()
{
	int flag_macro = 0;
	solve(flag_macro);	
	if (!print_eigenvalues())
		nrerror("OPmodel<OPFLOAT,basis>::compute_absorption: Cannot print eigenvalues!");
	if (!print_eigenvectors())
		nrerror("OPmodel<OPFLOAT,basis>::compute_absorption: Cannot print eigenvectors!");
	
	ostringstream ss;
	vector<string> dipole_files = print_dipoles();
	
	if (dipole_files.size()==0)
		nrerror("OPmodel<OPFLOAT,basis>::compute_absorption: Cannot compute absorption: no dipole files!");
	
	string epsilon_file;
	ss.str("");
	ss << model_name << "_epsilon.txt";
	epsilon_file=ss.str();
	
	vector<ublas::vector<complex<OPFLOAT> > > all_dipoles(0);
	vector<OPFLOAT> all_energies(0);
	
	OPFLOAT hw_min;
	OPFLOAT hw_max;
	OPFLOAT hw_step;
	OPFLOAT epsilon_inf; 
	OPFLOAT gamma_coeff;
	OPFLOAT d_slab;
	vector<OPFLOAT> inc_angle;
	vector<OPFLOAT> gamma_angle;
	OPFLOAT na;
	OPFLOAT nf;
	OPFLOAT cell_volume=basisPtr->get_latticePtr()->get_VCELL();
	
	read_spectrum_input_data(absorption_input_file, hw_min, hw_max, hw_step, epsilon_inf, gamma_coeff,
                             d_slab, inc_angle, gamma_angle, na, nf);
	
	write_epsilon(	dipole_files, epsilon_file, all_dipoles,
					  all_energies, hw_step,
					  epsilon_inf, cell_volume, gamma_coeff);
	
	write_absorption_spectrum(hw_min, hw_max, hw_step, d_slab, inc_angle, gamma_angle, na, nf,
							  model_name, epsilon_file);
	return;
}

template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::print(ofstream& log)
{
	time_t now = time(NULL);
	
	// Writes log file
	cout << "Printing model data ... " << endl;
    log << "Started computation at " << asctime(localtime(&now));
	log <<  endl;
	log << "==================================================================================================" << endl;
	log << "Model name: " << model_name << endl;
	log << "==================================================================================================" << endl << endl;
	log << "Type of calculation : " << calc << endl;
	log << "Real H matrix (1=yes): " << FLAG_REAL << endl;
	log << "N. of eigenvalues seeked (0=all): " << NEV << endl;
	log << "Davydov components to evaluate: ";
	for (int i=0; i<DCEV.size(); i++) log << DCEV[i] << "/t";
	
	string line;
	ifstream copy_file;
	if (calc[0]=='e' || calc[0]=='E')
	{
		log << "Emission input file:" << endl;
		copy_file.open(emission_input_file.c_str());
	}
	else
	{
		log << "Absorption input file:" << endl;
		copy_file.open(absorption_input_file.c_str());
	}
	
	while (!copy_file.eof())
	{
		getline(copy_file,line);
		log << "> " << line << endl;
	}
	
	basisPtr->print(log);
	basisPtr->get_latticePtr()->print(log);
	return;
}

/*
template<class OPFLOAT, class basis> 
bool OPmodel<OPFLOAT,basis>::print_results()
{
	vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr;
	ostringstream ss;
	int dc;
		
	if (DCEV.size()==0)
	{
		if (all_eigenstates!=0)
		{
		eigenstates_Ptr=all_eigenstates;
		ss.str("");
		ss << model_name << "_eigvec_all.txt";
		print_eigenvectors(ss.str(), eigenstates_Ptr);
		ss.str("");
		ss << model_name << "_eigval_all.txt";
		print_eigenvalues(ss.str(), eigenstates_Ptr);
		ss.str("");
		ss << model_name << "_dipoles_all.txt";
		print_dipoles(ss.str(), eigenstates_Ptr,-1);
		}
		else
			return(false);
	}
	else
	{
		for (dc=0; dc<DCEV.size(); dc++)
		{
			if (DC_eigenstates[dc]!=0)
			{
				eigenstates_Ptr=DC_eigenstates[dc];
				ss.str("");
				ss << model_name << "_eigvec_" << DCEV[dc] << ".txt";
				print_eigenvectors(ss.str(), eigenstates_Ptr);
				ss.str("");
				ss << model_name << "_eigval_" << DCEV[dc] << ".txt";
				print_eigenvalues(ss.str(), eigenstates_Ptr);
				ss.str("");
				ss << model_name << "_dipoles_" << DCEV[dc] << ".txt";
				print_dipoles(ss.str(), eigenstates_Ptr, dc);
			}
			else
				return(false);
		}
	}
	return(true);
}*/

template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::read_emission_input_data(string file_name, OPFLOAT& hw_min, OPFLOAT& hw_max, OPFLOAT& hw_step, 
							  OPFLOAT& width, vector<OPFLOAT>& T_vec, vector<Vector>& em_dir_vec, int& n0, int& n1, int& n2)
{	 
	string line;
	istringstream iss, iss2;
	ifstream infile(file_name.c_str());
	if (!infile.is_open()) nrerror("read_emission_input_data: File not found!");
	cout << "Reading emission input ..." << endl;
	T_vec.resize(0);
	em_dir_vec.resize(0);
	vector<OPFLOAT> dir_temp(0);
	
	get_next_from_file(infile, hw_min);
	get_next_from_file(infile, hw_max);
	get_next_from_file(infile, hw_step);
	get_next_from_file(infile, width);	
	
//	cout << hw_min << endl;
//	cout << hw_max << endl;
//	cout << hw_step << endl;
	getline(infile, line, ':');
	getline(infile, line);
	iss.str( line );
	copy( istream_iterator<OPFLOAT>( iss ), istream_iterator<OPFLOAT>(), back_inserter( T_vec ) );
	iss.clear();
//	cout << T_vec[3] << endl;
	getline(infile, line, ':');
	getline(infile, line);
	iss.str( line );
	copy( istream_iterator<OPFLOAT>( iss ), istream_iterator<OPFLOAT>(), back_inserter( dir_temp ) );
	iss.clear();
//	cout << dir_temp << endl;
	if (dir_temp.size()%3 !=0)
		nrerror("Wrong input dir");
	em_dir_vec.resize(dir_temp.size()/3);
	for (int d=0; d<em_dir_vec.size(); d++)
	{
		em_dir_vec[d].resize(3);
		for (int i=0; i<3; i++)
			em_dir_vec[d][i]=dir_temp[d*3+i];
	}
	get_next_from_file(infile, n0);
	infile >> n1;
	infile >> n2;

	infile.close();
	return;
};

/*	vector<int>		FinVibsMax(2);
	FinVibsMax[0]=3;	// Max n. of vibrations for dw mode
	FinVibsMax[1]=2;	// Max n. of vibrations for high energy mode*/

template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::compute_emission()
{
	if (DCEV.size()>0)
		nrerror("OPmodel<OPFLOAT,basis>::compute_emission: emission can only be computed for all states (No Davydov components)!");
	//  Reads parameters
	OPFLOAT hw_min, hw_max, hw_step;	
	OPFLOAT width;
	vector<OPFLOAT> T_vec;
	vector<Vector> em_dir_vec;
	int n0, n1,n2;
	read_emission_input_data(emission_input_file, hw_min, hw_max, hw_step, width, T_vec, em_dir_vec, n0, n1, n2);
	
	int i0,i1,i2,i,j,m, dc, d, n;
	int dwtot, hitot;
	Vector ka, kb, kc, k;
	int flag_macro=1;
	OPFLOAT hwLBE=1e10;
	OPFLOAT hwLBE_temp;
	vector<vector<OPFLOAT> > state_emission;
//	vector<int> lowest_state(4);
	Vector em_dir;
	OPFLOAT Tmax=0;
	for (i=0; i<T_vec.size(); i++)
		if (T_vec[i]>Tmax)
			Tmax=T_vec[i];
	

	//////////// Checks if emission data file exists
	ifstream file_in;
	ofstream file_out;
	ostringstream ss;
	int dwmax, himax;
	
	dwmax=basisPtr->get_MAX_VIB(0);
	if (basisPtr->get_NMODES()>1)
		himax=basisPtr->get_MAX_VIB(1);
	else
		himax=0;
	
	int n0_file, n1_file, n2_file;
	int dwmax_file, himax_file;
	bool emission_data_exist=false;
	
	ss.str("");
	ss << model_name << "_emission_data.txt";
	file_in.open(ss.str().c_str());
	if (file_in.is_open())
	{
		file_in >> n0_file;
		file_in >> n1_file;
		file_in >> n2_file;
		file_in >> dwmax_file;
		file_in >> himax_file;
		if (n0_file==n0 && n1_file==n1 && n2_file==n2 && dwmax_file==dwmax && himax_file==himax)
		{
			emission_data_exist=true;
			cout << "Emission data file exists" << endl;
		}
		file_in.close();
		
	}
	
	//////////// Writing emission intensities to file 
	
	if (!emission_data_exist)
	{
	cout << "Writing emission data file ... " << endl;
	
	file_out.open(ss.str().c_str());
	
	file_out << n0 << " " << n1 << " " << n2 << "  "  << dwmax << " " << himax << endl;
	// Controllare come mettere le dimensioni !!!!!!!!!!!
	// Ho cambiato e messo himax e dwmax veri (cioe' non ho aggiunto 1). Significa che i loop devono andare fino a <=

	for (i0=0; i0<n0; i0++)	
	{
		ka=basisPtr->get_latticePtr()->get_KCELL(0)*(i0-(n0-2+n0%2)/2)/n0;
		for (i1=0; i1<n1; i1++)
		{
			kb=basisPtr->get_latticePtr()->get_KCELL(1)*(i1-(n1-2+n1%2)/2)/n1;
			for (i2=0; i2<n2; i2++)
			{
				kc=basisPtr->get_latticePtr()->get_KCELL(2)*(i2-(n2-2+n2%2)/2)/n2;
				
				k=ka+kb+kc;
				
				if (norm_2(k)<ZERO_TOL)	
					k=(*((*basisPtr).get_latticePtr())).get_KLIGHT();
				
				(*basisPtr).change_KVEC(k);				
				solve(flag_macro);	
				

			
					std::sort(all_eigenstates->begin(),all_eigenstates->end());
					hwLBE_temp=(*all_eigenstates)[0].get_energy();
					if (hwLBE_temp<hwLBE)	
					{
						hwLBE=hwLBE_temp;
//						lowest_state[0]=i0;
//						lowest_state[1]=i1;
//						lowest_state[2]=i2;
//						lowest_state[3]=-1;
					}
					for (i=0; i<(*all_eigenstates).size(); i++)
					{
						if (Boltzmann_factor((*all_eigenstates)[i].get_energy()-hwLBE_temp, Tmax)>1e-3)
						{
							for (d=0; d<em_dir_vec.size(); d++)
							{
								em_dir=em_dir_vec[d];
								state_emission=basisPtr->compute_emission((*all_eigenstates)[i], em_dir);
								file_out << i0 << " " << i1 << " " << i2 << " " << d << " " 
								         << (*all_eigenstates)[i].get_energy() << " ";
								for (dwtot=0; dwtot<state_emission.size(); dwtot++)	
									for (hitot=0; hitot<state_emission[dwtot].size(); hitot++)		
									{
										if (i0==(n0-2+n0%2)/2 && i1==(n1-2+n1%2)/2 && i2==(n2-2+n2%2)/2 && dwtot+hitot==0) 					
											file_out << n0*n1*n2*state_emission[dwtot][hitot] << " ";
										else
											file_out << state_emission[dwtot][hitot] << " ";
									}											
								file_out << endl;
							}
						}
					}


			}
		}
	}
	file_out.close();
	}
	 
    ///////////////// Computing emission spectrum	
	ifstream infile;
	vector<vector<OPFLOAT> > Z(T_vec.size(), vector<OPFLOAT> (em_dir_vec.size(), 0.0) );
	OPFLOAT T;

	int t,nhw;
	OPFLOAT energy, hw, Si,em_energy;
	int ntothw=(hw_max-hw_min)/hw_step+1;
	vector<vector<vector<OPFLOAT> > > intensity(T_vec.size(), vector<vector<OPFLOAT> > (em_dir_vec.size(), vector<OPFLOAT> (ntothw,0.0)));
	vector<vector<vector<OPFLOAT> > > emission_data(0);
	vector<OPFLOAT> energy_data(0);
	vector<int> em_dir_data(0);
	
	cout << "Reading emission data ..." << endl;
	
//	ss.str("");
//	ss << model_name << "_emission_data.txt";
	infile.open(ss.str().c_str());
	infile >> n0;
	infile >> n1;
	infile >> n2;
	infile >> dwmax;
	infile >> himax;
	state_emission.resize(dwmax+1);
	for (dwtot=0; dwtot<=dwmax; dwtot++)
		state_emission[dwtot].resize(himax+1);
	
	hwLBE=1e10;
	while (!infile.eof())
	{
		infile >> i0;
		infile >> i1;
		infile >> i2; 
		infile >> d;
		infile >> energy;
		for (dwtot=0; dwtot<=dwmax; dwtot++)
			for (hitot=0; hitot<=himax; hitot++)
				infile >> state_emission[dwtot][hitot];		
		emission_data.push_back(state_emission);	
		energy_data.push_back(energy);
		em_dir_data.push_back(d);
		if (energy<hwLBE)
			hwLBE=energy;
	}
	infile.close();
	
	cout << "Computing emission spectra with hwLBE = " << hwLBE << " ..." << endl;
	for (i=0; i<emission_data.size(); i++)
	{
		for (t=0; t<T_vec.size(); t++) 
		{	
			T=T_vec[t];
			d=em_dir_data[i];
			em_dir=em_dir_vec[d];
			Z[t][d]+=Boltzmann_factor(energy_data[i]-hwLBE, T);
//			cout << energy_data[i]-hwLBE << " " << T << " " << Boltzmann_factor(energy_data[i]-hwLBE, T) << endl;
			for (nhw=0; nhw<ntothw; nhw++)
			{
					hw=hw_min+nhw*hw_step;
					Si=0.0;
					for (dwtot=0; dwtot<=dwmax; dwtot++)	
					for (hitot=0; hitot<=himax; hitot++)
					{
						em_energy=energy_data[i]-dwtot*(*basisPtr).get_vib_En_0()-hitot*(*basisPtr).get_vib_En_1();
						Si+=emission_data[i][dwtot][hitot]*pow(em_energy/hwLBE,3)*Gaussian_dist(hw-em_energy, width)/Gaussian_dist(OPFLOAT(0.0), width);
					}
					intensity[t][d][nhw]+=Boltzmann_factor(energy_data[i]-hwLBE, T)*Si;					
			}
		}	
	}

//////////////////// Writes output file
	cout << "Writing emission spectra ..." << endl;
	for (t=0; t<T_vec.size(); t++) 
	{
		T=T_vec[t];
		for (d=0; d<em_dir_vec.size(); d++) 
		{
			ss.str("");
			ss << model_name << ".emission." << T << "." << d << ".txt";
			file_out.open(ss.str().c_str());
			for (nhw=0; nhw<ntothw; nhw++)
				file_out << hw_min+nhw*hw_step << "\t" << intensity[t][d][nhw]/Z[t][d] << endl;
			file_out.close();
		} 
	}
	return;
}

/*
template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::compute_emission(string input_filename, int dc, OPFLOAT hwLBE)
{
	vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr;
	if (dc<0)
	{
		eigenstates_Ptr=all_eigenstates;
	}
	else
	{
		if (dc<DC_eigenstates.size())
		{
			eigenstates_Ptr=DC_eigenstates[dc];
		}
		else
			nrerror("OP_model::read_eigenstates: Davydov component not found!");
	}
	cout << "Computing emission ... " << endl;


// Computes emission intensities
	int i,j,t,d;
	OPFLOAT Z=0.0;
	OPFLOAT dummy;
	Vector em_dir;
	OPFLOAT T;
	vector<vector<vector<OPFLOAT> > > states_emission(0);
	vector<OPFLOAT> states_energy(0);

	std::sort((*eigenstates_Ptr).begin(),(*eigenstates_Ptr).end());

	for (t=0; t<T_vec.size(); t++) for (d=0; d<em_dir_vec.size(); d++)
	{
		T=T_vec[t];
		em_dir=em_dir_vec[d];
		cout << T << em_dir << endl;
		if (T<ZERO_TOL)
		{
			Z=1.0;
//			hwLBE=(*eigenstates_Ptr)[0].get_energy();
			states_emission.push_back((*basisPtr).compute_emission((*eigenstates_Ptr)[0], em_dir));
			states_energy.push_back((*eigenstates_Ptr)[0].get_energy());
		}
		else
		{
//			hwLBE=(*eigenstates_Ptr)[0].get_energy();
			for (i=0; i<(*eigenstates_Ptr).size(); i++)
			{
				dummy=Boltzmann_factor((*eigenstates_Ptr)[i].get_energy()-hwLBE, T);
				if (dummy>1e-3)
				{
					states_emission.push_back((*basisPtr).compute_emission((*eigenstates_Ptr)[i], em_dir));
					states_energy.push_back((*eigenstates_Ptr)[i].get_energy());
					Z+=Boltzmann_factor((*eigenstates_Ptr)[i].get_energy()-hwLBE, T);
				}
			}
		}	

//		for (i=0; i<eigenstates[0].get_emission().size(); i++) for (j=0; j<eigenstates[0].get_emission()[i].size(); j++)
//			cout << i << "\t" << j << "\t" << eigenstates[0].get_emission()[i][j] << endl;

// Computes emission spectrum
		cout << "Computing emission spectrum ..." << endl;

		OPFLOAT hw;
		OPFLOAT Si,S;
		OPFLOAT em_energy=0.0;
		int dwtot, hitot;	
		ostringstream ss;
		ss.str("");
		if (dc<0)
			ss << model_name << "emission_" << T << "_dir" << d << ".txt";
		else
			ss << model_name << "emission_" << dc << "_" << T << "_dir" << d << ".txt";
		cout << ss.str() << endl;
		ofstream file_out;
		file_out.open(ss.str().c_str());
		for (hw=hw_min; hw<=hw_max; hw+=hw_step)
		{
			S=0.0;
			for (i=0; i<states_emission.size(); i++)
			{
					Si=0.0;
					for (dwtot=0; dwtot<states_emission[i].size(); dwtot++)	
					for (hitot=0; hitot<states_emission[i][dwtot].size(); hitot++)
					{
						em_energy=states_energy[i]-dwtot*(*basisPtr).get_vib_En_0()-hitot*(*basisPtr).get_vib_En_1();
						Si+=states_emission[i][dwtot][hitot]*pow(em_energy/hwLBE,3)*Gaussian_dist(hw-em_energy, width)/Gaussian_dist(OPFLOAT(0.0), width);
					}
					S+=1.0/Z*Boltzmann_factor(states_energy[i]-hwLBE, T)*Si;
			}
			file_out << hw << "\t" << S << endl;
		}
	}
	return;
}*/


///////////////////////////////////////////////////////////
////   ARPACK ROUTINES

template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::solve_real_nev(int dc, int flag_macro)
{
	return;
}
template<class OPFLOAT, class basis> 
void OPmodel<OPFLOAT,basis>::solve_complex_nev(int dc, int flag_macro)
{
	return;
}

/*
 template<class OPFLOAT, class basis> 
 void OPmodel<OPFLOAT,basis>::solve_complex_nev(int dc, int flag_macro)
 {
 cout << "Solving problem for dc=" << dc << endl;
 int n;	// Dimension of the problem.
 int i,j;
 vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr;
 if (dc<0)
 {
 n = (*basisPtr).get_BASIS_size();
 eigenstates_Ptr=all_eigenstates;
 }
 else
 {
 if (dc<DC_eigenstates.size())
 {
 n = (*basisPtr).get_BD_BASIS_size(dc);
 eigenstates_Ptr=DC_eigenstates[dc];
 }
 else
 nrerror("OP_model::solve_BASIS: Davydov component not found!");
 }
 cout << "n=" << n << endl;
 
 if (NEV<=0 || NEV>=n)
 nrerror("OP_model::solve_complex_nev: Wrong nev!");
 
 // Solve with ARPACK
 cout << "Finding the first " << NEV << " eigenvalues of the complex H with ARPACK ..." << endl;
 int				nnz;	// Number of nonzero elements in A.
 vector<int>		irow;	// pointer to an array that stores the row
 // indices of the nonzeros in A.
 vector<int>		pcol;	// pointer to an array of pointers to the
 // beginning of each column of A in valA.
 vector<Complex> valH;	// pointer to an array that stores the
 // nonzero elements of A.
 // Creating a sparse complex matrix.
 nnz=CreateCSCNonSymComplexMatrix(n, valH, irow, pcol, ZERO_TOL, dc);
 ARluNonSymMatrix<Complex, OPFLOAT > H(n, nnz, &valH[0], &irow[0], &pcol[0]);
 
 // Defining the eigenvalue problem. 
 
 ARluCompStdEig<OPFLOAT> prob(NEV, H); 
 prob.ChangeWhich("SM");
 //			prob.ChangeNcv(n-1);
 //			prob.ChangeMaxit(1000*nev);
 // Declaring output variables.
 vector<Complex>		ARPACK_eigenvalues(NEV);		// Eigenvalues.
 vector<Complex>		ARPACK_eigenvectors(NEV*n);		// Eigenvectors.
 
 // Finding eigenvalues and eigenvectors.
 
 bool ischur=false;
 int conv_eig;
 Complex*   ValPtr= &ARPACK_eigenvalues[0];
 Complex*   VecPtr= &ARPACK_eigenvectors[0];	
 
 conv_eig  = prob.EigenValVectors(VecPtr, ValPtr, ischur);
 cout << "NEV: " << NEV << " conv_eig: " << conv_eig << endl;
 
 (*eigenstates_Ptr).resize(conv_eig);
 for (i=0; i<conv_eig; i++)
 {
 if (abs(ARPACK_eigenvalues[i].imag()/ARPACK_eigenvalues[i].real())>IMAG_TOL)
 nrerror("OPmodel::solve(): ARPACK returned imaginary eigenvalues!");
 (*eigenstates_Ptr)[i].set_energy(ARPACK_eigenvalues[i].real());
 (*eigenstates_Ptr)[i].resize(n);
 for (j=0; j<n; j++)
 (*eigenstates_Ptr)[i].set_coeff(j,ARPACK_eigenvectors[i*n+j]);
 }
 
 cout << "Eigenvectors stored" << endl;
 return;
 }*/


/*
 int CreateCSCNonSymComplexMatrix(int n, vector<Complex> &A, vector<int> &irow, vector<int> &pcol, 
 OPFLOAT zero_tol, int dc)
 {
 int   row, col, first_done;
 Complex  val=0.0;
 A.resize(0);
 irow.resize(0);
 pcol.resize(0);
 
 int nnz=0;
 pcol.push_back(0);
 first_done=1;
 for (col=0; col<n; col++)
 {	  
 for (row=0; row<n; row++)
 {
 val=(*basisPtr).Hint(row,col,dc);
 if (abs(val)>zero_tol)
 {
 A.push_back(val);
 irow.push_back(row);
 if (first_done==0)
 {
 pcol.push_back(A.size()-1);
 first_done=1;
 }
 nnz++;
 }
 }
 if (first_done==0)
 pcol.push_back(A.size()-1);
 first_done=0;
 }
 pcol.push_back(nnz);
 return(nnz);
 } // CSCNonSymmetricComplexMatrix A with function to compute A[i][j].
 
 Complex GetCSCNonSymComplexMatrix(vector<Complex> &A, vector<int> &irow, vector<int> &pcol, int i, int j)
 {
 int pos;
 for (pos=pcol[j]; pos<pcol[j+1]; pos++)
 if (irow[pos]==i)
 return(A[pos]);
 return(0.0);
 }
 */



/*
 int CreateCSCSymRealMatrix(int n, vector<OPFLOAT> &A, vector<int> &irow, 
 vector<int> &pcol, OPFLOAT zero_tol, int dc)
 {
 int   row, col, first_done;
 OPFLOAT  val=0.0;
 Complex comp_val;
 A.resize(0);
 irow.resize(0);
 pcol.resize(0);
 
 int nnz=0;
 pcol.push_back(0);
 first_done=1;
 for (col=0; col<n; col++)
 {	  
 for (row=0; row<=col; row++)
 {
 comp_val=(*basisPtr).Hint(row,col,dc);
 if (abs(comp_val.imag()/comp_val.real())>IMAG_TOL)
 nrerror("CreateCSCSymRealMatrix: H matrix element is not real!");
 val=comp_val.real();
 if (abs(val)>zero_tol)
 {
 A.push_back(val);
 irow.push_back(row);
 if (first_done==0)
 {
 pcol.push_back(A.size()-1);
 first_done=1;
 }
 nnz++;
 }
 }
 if (first_done==0)
 pcol.push_back(A.size()-1);
 first_done=0;
 }
 pcol.push_back(nnz);
 return(nnz);
 } */


/*
 template<class OPFLOAT, class basis> 
 void OPmodel<OPFLOAT,basis>::solve_real_nev(int dc, int flag_macro)
 {
 cout << "Solving problem for dc=" << dc << endl;
 int n;	// Dimension of the problem.
 int i,j;
 vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr;
 if (dc<0)
 {
 n = (*basisPtr).get_BASIS_size();	// Dimension of the problem.
 eigenstates_Ptr=all_eigenstates;
 }
 else
 {
 if (dc<DC_eigenstates.size())
 {
 n = (*basisPtr).get_BD_BASIS_size(dc);	// Dimension of the problem.
 eigenstates_Ptr=DC_eigenstates[dc];
 }
 else
 nrerror("OP_model::solve_BASIS: Davydov component not found!");
 }
 cout << "n=" << n << endl;
 
 //	int n = (*basisPtr).get_BASIS_size();	// Dimension of the problem.
 //	int i,j;
 
 
 if (NEV<=0 || NEV>=n)
	nrerror("OP_model::solve_real_nev: Wrong nev!");
 
 // Solve with ARPACK
 cout << "Finding the first " << NEV << " eigenvalues of the real H with ARPACK ..." << endl;
 int				nnz;	// Number of nonzero elements in A.
 vector<int>		irow;	// pointer to an array that stores the row
 // indices of the nonzeros in A.
 vector<int>		pcol;	// pointer to an array of pointers to the
 // beginning of each column of A in valA.
 vector<OPFLOAT> valH;	// pointer to an array that stores the
 // nonzero elements of A.
 char   uplo='U';		// Variable that indicates whether the upper
 // (uplo='U') or the lower (uplo='L') part of
 // H will be stored in A, irow and pcol.
 
 // Creating a sparse symmetric real matrix.
 cout << "Computing H matrix ..." << endl;
 nnz = CreateCSCSymRealMatrix(n, valH, irow, pcol, ZERO_TOL, dc);
 cout << nnz << " " << float(nnz/(n*n)*100) << "%" << endl;
 cout << "Creating ARPACK matrix ..." << endl;
 ARluSymMatrix<OPFLOAT> H(n, nnz, &valH[0], &irow[0], &pcol[0], uplo);
 
 // Defining the eigenvalue problem. 
 
 cout << "Creating ARPACK problem ..." << endl;
 ARluSymStdEig<OPFLOAT> prob(NEV, H);
 prob.ChangeWhich("SM");
 //		prob.ChangeNcv(n-1);
 //		prob.ChangeMaxit(1000*nev);
 // Declaring output variables.
 vector<OPFLOAT>		ARPACK_eigenvalues(NEV);		// Eigenvalues.
 vector<OPFLOAT>		ARPACK_eigenvectors(NEV*n);		// Eigenvectors.
 
 // Finding eigenvalues and eigenvectors.
 
 cout << "Ncv=" << prob.GetNcv() << endl;
 cout << "Maxit=" << prob.GetMaxit() << endl;
 bool ischur=false;
 int conv_eig;
 OPFLOAT*   ValPtr= &ARPACK_eigenvalues[0];
 OPFLOAT*   VecPtr= &ARPACK_eigenvectors[0];
 
 conv_eig  = prob.EigenValVectors(VecPtr, ValPtr, ischur);
 cout << "NEV: " << NEV << " conv_eig: " << conv_eig << endl;
 
 (*eigenstates_Ptr).resize(conv_eig);
 for (i=0; i<conv_eig; i++)
 {
 (*eigenstates_Ptr)[i].set_energy(ARPACK_eigenvalues[i]);
 (*eigenstates_Ptr)[i].resize(n);
 for (j=0; j<n; j++)
 (*eigenstates_Ptr)[i].set_coeff(j,ARPACK_eigenvectors[i*n+j]);
 }
 return;
 }
 */

/*
 // Controllare questa routine
 OPFLOAT GetCSCSymRealMatrix(vector<Complex> &A, vector<int> &irow, vector<int> &pcol, int i, int j)
 {
 int n;
 if (j>=i)
 {
 n=irow[i];
 while (n<irow[i+1])
 {		
 if (pcol[n]>j)
 return(0.0);
 if (pcol[n]==j)
 return(A[n]);
 n++;
 }
 return(0.0);
 }
 else 
 return(GetCSCSymRealMatrix(j,i));
 }*/

// Explicit instantiation
template class OPmodel<float, PhononCloud<float> >;
template class OPmodel<double, PhononCloud<double> >;

