/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _OP_MODEL_H
#define _OP_MODEL_H

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
//#include <N_Functions.h> 
//#include <physics.h>
#include <lattice.h>
#include <QuantumState.h>
/*#include "arcomp.h" 
#include "arlnsmat.h" // ARluNonSymMatrix definition.
#include "arlscomp.h" // ARluCompStdEig definition.
#include "arlssym.h " // ARluSymStdEig definition.
#include "arlsmat.h " // ARluSymMatrix definition.*/

template<class OPFLOAT, class basis>
class OPmodel 
{
	typedef std::complex<OPFLOAT>							Complex;
	typedef ublas::vector<OPFLOAT>							Vector;
	typedef ublas::vector<Complex>							ComplexVector;
	typedef ublas::matrix<OPFLOAT,	ublas::column_major>	Matrix;
	typedef ublas::matrix<Complex,	ublas::column_major>	ComplexMatrix;
	
	private:
		basis*				basisPtr;
		std::string			model_name;
		std::string			calc;
		std::string			emission_input_file;
		std::string			absorption_input_file;
		int					NEV;
		std::vector<int>	DCEV;
		int					FLAG_REAL;
		std::vector<QuantumState<OPFLOAT> >*				all_eigenstates;
		std::vector<std::vector<QuantumState<OPFLOAT> >* >	DC_eigenstates;
	
	void initialize(std::string input_file);
	bool print_eigenvectors(std::string eigvec_filename, std::vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr);
	bool print_eigenvalues(std::string eigval_filename, std::vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr);
	bool print_dipoles(std::string dipoles_filename, std::vector<QuantumState<OPFLOAT> >*  eigenstates_Ptr, int dc);
	void solve(int flag_macro);
	void print(std::ofstream& log);
	void compute_emission();
	void compute_absorption();

//  CSCNonSymmetricComplexMatrix A with function to compute A[i][j].
//	int CreateCSCNonSymComplexMatrix(int n, std::vector<Complex> &A, std::vector<int> &irow, std::vector<int> &pcol, OPFLOAT zero_tol, int dc); 
//	Complex GetCSCNonSymComplexMatrix(std::vector<Complex> &A, std::vector<int> &irow, std::vector<int> &pcol, int i, int j);
//	int CreateCSCSymRealMatrix(int n, std::vector<OPFLOAT> &A, std::vector<int> &irow, std::vector<int> &pcol, OPFLOAT zero_tol, int dc);

	void solve_complex_full(int dc, int flag_macro);
	void solve_real_full(int dc, int flag_macro);
	void solve_complex_nev(int dc, int flag_macro);
	void solve_real_nev(int dc, int flag_macro);
	void read_emission_input_data(std::string file_name, OPFLOAT& hw_min, OPFLOAT& hw_max, OPFLOAT& hw_step, 
								  OPFLOAT& width, std::vector<OPFLOAT>& T_vec, std::vector<Vector>& em_dir_vec, 
								  int& n0, int& n1, int& n2);

	public:
	// Long Constructor
	OPmodel(std::string input_file) {initialize(input_file);}
	virtual ~OPmodel( void ) { 
		delete all_eigenstates;
		for (int i=0; i<DC_eigenstates.size(); i++)
			 delete DC_eigenstates[i];
			 delete basisPtr;
			 return;
			 };

	void solve();
	bool print_eigenvectors();
	bool print_eigenvalues();
	std::vector<std::string> print_dipoles();
//	int print_dipoles(std::string& dipoles_filename);
//	bool read_eigenstates(std::string eigvec_filename);
//	int read_eigenstates(std::string eigvec_filename, int dc);


//	void compute_emission(std::string input_filename, int dc, OPFLOAT hwLBE);

};  // end of class OPmodel

#endif
