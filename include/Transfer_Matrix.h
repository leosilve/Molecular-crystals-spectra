/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _TRANSFER_MATRIX_H
#define _TRANSFER_MATRIX_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <math.h>
#include <vector>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <OPutils.h>
#include <N_Functions.h> 

namespace ublas = boost::numeric::ublas;

typedef std::complex<double>							dComplex;
typedef ublas::vector<double>							dVector;
typedef ublas::vector<dComplex>							dComplexVector;
typedef ublas::matrix<double,	ublas::column_major>	dMatrix;
typedef ublas::matrix<dComplex,	ublas::column_major>	dComplexMatrix;

// Note: if there is an extra (empty) line in the data file this routine will put 1 wrong extra term
// DataFile must have 13 columns: energy (in eV), Re[eps_11], Im[eps_11], Re[eps_12], Im[eps_12],
// Re[eps_13], Im[eps_13], Re[eps_22], Im[eps_22], Re[eps_23], Im[eps_23], Re[eps_33], Im[eps_33]
// Epsilon is then assumed to be symmetric
void read_epsilon_data(std::string DataFile,
					   std::vector<double>& energy,
					   std::vector<std::vector<dComplex> >& epsilon_Data	);

// Given the energy, this routine finds the correct epsilon by interpolating data 
void set_epsilon( double hw, dComplexMatrix& epsilon, std::vector<double>& energy, std::vector<std::vector<dComplex> >&  epsilon_Data);
// A routine to manually set a diagonal epsilon
void set_diag_epsilon( dComplex e11, dComplex e22, dComplex e33, dComplexMatrix&  eps);
// A routine to manually set an isotropic epsilon
void set_iso_epsilon( dComplex e11, dComplexMatrix& eps);

// Checks if the used epsilons are symmetric
int check_symmetry(dComplexMatrix&  eps);

// Checks if epsilon corresponds to an isotropic material
int check_isotropy(dComplexMatrix&  eps);

// This routine computes the Tp matrix (see Schubert) 
// It also stores all the properties of the 4 waves propagating in the material
// kz (kz[i] contains the value for wave i=0,..,3), 
// polarization (pol[i] is a vector containing {Ex, Ey, Hx, Hy}, for each wave i=0,..,3)
// Poynting vector direction (the vectors s[i], for each wave i=0,..,3) 
void set_Tp(  
			dComplexMatrix&  Tp,
			dComplexMatrix&  eps,
			double kx,
			double hw,
			double d,
			std::vector<dComplex>& kz,
			std::vector<dComplexVector>& pol,
			std::vector<dComplexVector>& s,
			int info_flag       // if different from 0 data about the waves are printed on a log file
			);

// Computes matrix T (see Schubert)
void set_T(					   std::vector<dComplexMatrix >&  epsilon,
							   dComplexMatrix&  T,
							   dComplexMatrix& La,
							   dComplexMatrix& Lf,
							   double kx,
							   double hw,
							   std::vector<double>& d,
							   int info_flag
		   );

// Rotates dMatrix M1 and puts the rotated matrix in M2
// alpha, beta , gamma are Euler angles in degrees
void rotate_matrix(	 dComplexMatrix&  M1, 
					 dComplexMatrix&  M2, 
					 double alpha,
					 double beta,
					 double gamma);

///////////////////////////////////////////////////////////////////////////////////////////
// Useful function to access coefficients and Poynting vectors

dComplex Bs(dComplexMatrix& T, dComplex As, dComplex Ap);
dComplex Bp(dComplexMatrix& T, dComplex As, dComplex Ap);
dComplex Cs(dComplexMatrix& T, dComplex As, dComplex Ap);
dComplex Cp(dComplexMatrix& T, dComplex As, dComplex Ap);

double  SzI( dComplexMatrix& T, 
			dComplexMatrix& La,
			dComplexMatrix& Lf,
			dComplex As, 
			dComplex Ap);

double  SzRs( dComplexMatrix& T, 
			 dComplexMatrix& La,
			 dComplexMatrix& Lf,
			 dComplex As, 
			 dComplex Ap);

double  SzRp( dComplexMatrix& T, 
			 dComplexMatrix& La,
			 dComplexMatrix& Lf,
			 dComplex As, 
			 dComplex Ap);

double  SzTs( dComplexMatrix& T, 
			 dComplexMatrix& La,
			 dComplexMatrix& Lf,
			 dComplex As, 
			 dComplex Ap);

double  SzTp( dComplexMatrix& T, 
			 dComplexMatrix& La,
			 dComplexMatrix& Lf,
			 dComplex As, 
			 dComplex Ap);
					 
double  Rs_coeff(dComplexMatrix& T, 
				 dComplexMatrix& La,
				 dComplexMatrix& Lf,
				 dComplex As, 
				 dComplex Ap);

double  Rp_coeff(	  dComplexMatrix& T, 
				 dComplexMatrix& La,
				 dComplexMatrix& Lf,
				 dComplex As, 
				 dComplex Ap);

double  Ts_coeff(	  dComplexMatrix& T, 
				 dComplexMatrix& La,
				 dComplexMatrix& Lf,
				 dComplex As, 
				 dComplex Ap);

double  Tp_coeff(	  dComplexMatrix& T, 
				 dComplexMatrix& La,
				 dComplexMatrix& Lf,
				 dComplex As, 
				 dComplex Ap);

///////////////////////////////////////////////////////////////////////////////////////////
// Example of a routine that computes absorption for a slab 
// placed between two isotropic media, as a function of energy. 
void spectrum_vs_energy(std::string DataFile,			// name of the file where epsilon data are
						std::string OutputFile,			// name of the file where the spectrum is written
						double na,                  // ambient n
						double nf,					// substrate n
						double d_slab,				// thickness of the first slab
						double inc_angle,			// angle of incidence (°)
						double gamma				// angle between the incident plane and xz plane 
						);

template<typename T>
void read_spectrum_input_data(std::string file_name, T& hw_min, T& hw_max, T& hw_step, T& epsilon_inf, T& gamma_coeff, T& d_slab, 
							  std::vector<T>& inc_angle, std::vector<T>& gamma_angle, T& na, T& nf);

template<typename T>
void read_spectrum_input_data_full(std::string file_name, T& hw_min, T& hw_max, T& hw_step, T& epsilon_inf, T& gamma_coeff, T& d_slab, 
								   std::vector<T>& inc_angle, std::vector<T>& gamma_angle, T& na, T& nf, 
								   T& cell_volume, std::vector<std::string>& dipole_files, std::string& epsilon_file);

template <typename T>
void write_absorption_spectrum(T hw_min, T hw_max,  T hw_step, T d_slab,
							   std::vector<T> inc_angle,
							   std::vector<T> gamma_angle,
							   T na,
							   T nf,
							   std::string output_filename, std::string epsilon_filename);

template<typename T>
void read_dipoles( std::vector<std::string> input_files, std::vector<ublas::vector<std::complex<T> > > &all_dipoles,
				  std::vector<T> &all_energies );

template<typename T>
bool write_epsilon(	std::vector<std::string> input_files, std::string output_file, std::vector<ublas::vector<std::complex<T> > > &all_dipoles,
				   std::vector<T> &all_energies, T hw_step,
				   T epsilon_inf, T cell_volume, T gamma_coeff);

template<typename T>
bool write_epsilon(	std::vector<std::string> input_files, std::string output_file, std::vector<ublas::vector<std::complex<T> > > &all_dipoles,
				   std::vector<T> &all_energies, T hw_min, T hw_max, T hw_step,
				   T epsilon_inf, T cell_volume, T gamma_coeff);


#endif
	  
