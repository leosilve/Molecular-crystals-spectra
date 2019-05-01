/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _PHONONCLOUD_H
#define _PHONONCLOUD_H

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
#include <algorithm>
#include <OPutils.h>
#include <pos.h>
#include <N_Functions.h> 
#include <ME.h>
#include <physics.h>
#include <lattice.h>
#include <QuantumState.h>
#include <mpstate.h>

namespace ublas = boost::numeric::ublas;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////    Phonon Cloud class (two modes: mode 0 is the DW and mode 1 is the high energy)     /////////  
/////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class FLOAT>
class PhononCloud
{

	typedef std::complex<FLOAT>							    Complex;
	typedef ublas::vector<FLOAT>							Vector;
	typedef ublas::vector<Complex>							ComplexVector;
	typedef ublas::matrix<FLOAT,	ublas::column_major>	Matrix;
	typedef ublas::matrix<Complex,	ublas::column_major>	ComplexMatrix;

private:

	lattice<FLOAT>*	latticePtr;	// pointer to the crystal lattice class

	FLOAT		lambda_0;		// Huang-Rhys factor (for each species)
	FLOAT		lambda_1;		// Huang-Rhys factor high energy mode (for each species)
	FLOAT		elec_En;		// Elec transition energy  (for each species)
	FLOAT		vib_En_0;		// Vibration energy (for each species)
	FLOAT		vib_En_1;		// Vibration high energy mode (for each species)

	int											N_SYM_CLOUDS;	// Num of symmetric clouds for each Davydov component
	std::vector<std::vector<int>	>			MAX_VIBS;		// Max n. of vibrations at each NN for multiple modes
	std::vector<ME<FLOAT> >						BASIS;			// Basis states
	std::vector<std::vector<MPSTATE<FLOAT> > >	BD_BASIS;		// Symmetrized basis set
	Vector						KVEC;			// wave vector
	std::string					prob_name;		// Name of the PhononCloud problem
	FLOAT						EXP_FAC;		// factor appearing in the calc of MAX_VIBS
	std::vector<int>			MAX_VIB;		// total maximum of vibrations per exciton for each mode
	int							NP;				// Maximum number of particles for the problem
	int							DC;				// Davydov component solved for
	int							FLAG_DW;		// If ==1 considers a symmetric FLOAT well vibrational potential for the ground state molecule
	int							FLAG_2D;
	int							NMODES;
	int							MAX_CLOUD;
	FLOAT						CLOUD_RADIUS;
	std::vector<std::vector<pos<FLOAT> > >	CLOUD_NN;			// CLOUD_NN contains the nearest neighbours positions (a 3d vector)

// Symmetric double well parameters
	double						DW_hw0;			// in eV
	double						DW_alpha;		// in eV
	double						DW_A;			// in eV 
	double						DW_Vd;			// in eV (energy of the FLOAT well bottom)
	double						DW_hwe;			// in eV 
	double						DW_lambdae;		// displacement in excited state natural units
	unsigned					DW_basis_dim;
	unsigned					DW_nug_max;
	unsigned					DW_nue_max;
	std::vector<std::vector<double> >		DW_FCF_table;	// Table of FCF
	std::vector<double>						DW_levels;		// List of DW energy levels

// Private functions
	
	void		initialize(std::string input_file, std::string lattice_file);
	void		set_n_vib(int cloud_pos, int v0_max, int v1_max, ME<FLOAT>& exc); 
	void		set_n_vib_one_mode(int cloud_pos, int v0_max, ME<FLOAT>& exc);
	ME<FLOAT>	compute_symmetric_state(int sym, ME<FLOAT>& exc);
	void		compute_BASIS();
	FLOAT		FC_overlap(int mode, int nue, int nug);
	FLOAT		Vibrational_energies(int mode, int flag_exc, int n);
	
public:
// Long Constructor
	PhononCloud(std::string input_file, std::string lattice_file) { initialize(input_file, lattice_file); };
	virtual ~PhononCloud( void ) {
		delete latticePtr;
		return;
	}

// Functions to access private members
	FLOAT		get_lambda_0() const { return(lambda_0); };
	FLOAT 		get_lambda_1() const { return(lambda_1); };
	FLOAT 		get_elec_En()  const { return(elec_En); };
	FLOAT 		get_vib_En_0() const { return(vib_En_0); };
	FLOAT 		get_vib_En_1() const { return(vib_En_1); };
	
	FLOAT							get_EXP_FAC()					const {return(EXP_FAC); };
	std::vector<int>				get_MAX_VIB()					const {return(MAX_VIB); };
	int								get_MAX_VIB(int i)				const {return(MAX_VIB[i]); };
	std::vector<std::vector<int> >	get_MAX_VIBS()					const {return(MAX_VIBS); };
	int								get_MAX_VIBS(int mode, int i)	const {return(MAX_VIBS[mode][i] ); };
	int								get_BASIS_size()				const {return(BASIS.size());}
	ME<FLOAT>						get_BASIS(int i)				const {return(BASIS[i]); };
	MPSTATE<FLOAT>					get_BD_BASIS(int i, int j)		const {return(BD_BASIS[i][j]); };
	int								get_BD_BASIS_size()				const {return(BD_BASIS.size()); };
	int								get_BD_BASIS_size(int i)		const {return(BD_BASIS[i].size()); };
	int								get_NMODES()					const {return(NMODES); };
	int								get_NP()						const {return(NP); };
	int								get_MAX_CLOUD()					const {return(MAX_CLOUD); };
	Vector							get_KVEC()						const {return(KVEC); };
	int								get_FLAG_DW()					const {return(FLAG_DW); };
	int								get_CLOUD_RADIUS()				const {return(CLOUD_RADIUS); };
	
	lattice<FLOAT>*					get_latticePtr()				const {return(latticePtr);};
	ComplexVector					get_BASIS_dipole(int i)				  {return(get_dipole(BASIS[i]));}
	ComplexVector					get_BD_BASIS_dipole(int dc, int i)	  {return(get_dipole(BD_BASIS[dc][i]));}

	
	ComplexVector					get_dipole(ME<FLOAT>& exc);
	ComplexVector					get_dipole(MPSTATE<FLOAT>& mps);
	
	int								get_CLOUD_NN_size()				const { return(CLOUD_NN.size()); };
	int								get_CLOUD_NN_size(int i)		const { return(CLOUD_NN[i].size()); };
	
	pos<FLOAT> get_CLOUD_NN(int i, int j) const
	{ 
		if (i<0 || i>=(*latticePtr).get_N_MOL() || j<0 || j>=get_CLOUD_NN_size(i))
			nrerror("get_CLOUD_NN: index out of bounds");
		return(CLOUD_NN[i][j]);
	};
	
	int get_CLOUD_NN_alpha(int i, int j) const
	{ 
		if (i<0 || i>=(*latticePtr).get_N_MOL() || j<0 || j>=get_CLOUD_NN_size(i))
			nrerror("get_CLOUD_NN: index out of bounds");
		return(CLOUD_NN[i][j].alpha);
	};
	
	Vector get_CLOUD_NN_n(int i, int j) const
	{ 
		if (i<0 || i>=(*latticePtr).get_N_MOL() || j<0 || j>=get_CLOUD_NN_size(i))
			nrerror("get_CLOUD_NN: index out of bounds");
		return(CLOUD_NN[i][j].n);
	};
	
	std::vector<std::vector<FLOAT> > compute_emission(QuantumState<FLOAT>&	emitting_state, Vector& em_dir);
//	FLOAT compute_emission_I0(QuantumState<FLOAT>&	emitting_state, Vector& em_dir);
//	FLOAT compute_emission_I1(QuantumState<FLOAT>&	emitting_state, Vector& em_dir, int vtot);
//	FLOAT compute_emission_I2(QuantumState<FLOAT>&	emitting_state, Vector& em_dir, int v1, int v2);
	Complex Hint(ME<FLOAT> &exc1, ME<FLOAT> &exc2, int flag_macro);
	Complex Hint(MPSTATE<FLOAT>  &mps1, MPSTATE<FLOAT>  &mps2, int flag_macro);
	Complex Hint(int i, int j, int DC, int flag_macro);
	bool check_BD_decomposition(int dc1, int dc2, int flag_macro, FLOAT tol);
	
	void print(std::ofstream& log);
	bool change_KVEC(Vector new_KVEC);
	bool change_KVEC(FLOAT k0, FLOAT k1, FLOAT k2);
	void compute_BD_BASIS();
	
	

};	// End of PhononCloud class


#endif
