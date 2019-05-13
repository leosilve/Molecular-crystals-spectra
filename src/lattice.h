/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _LATTICE_H
#define _LATTICE_H

#include <boost/math/special_functions/erf.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <Complex>
#include <math.h>
#include <vector>
#include <OPutils.h>
#include <pos.h>
#include <sym_op.h>

namespace ublas = boost::numeric::ublas;

//////////////////////////////////////////
///     Molecular crystal's lattice.
///
///		It defines all the properties of the molecular crystal's lattice
///
//////////////////////////////////////////

template <typename LATFLOAT>
class lattice {
	
	typedef std::complex<LATFLOAT>							Complex;
	typedef ublas::vector<LATFLOAT>							Vector;
	typedef ublas::vector<Complex>							ComplexVector;
	typedef ublas::matrix<LATFLOAT,	ublas::column_major>	Matrix;
	typedef ublas::matrix<Complex,	ublas::column_major>	ComplexMatrix;
	typedef std::vector<std::vector<Complex> >				ComplexDoubleVector;

protected:

 // Lattice parameters.
	int				N_MOL;			///< Number of molecules per cell
	int				N_SP;			///< Number of crystallographic species
	int				N_SYM_OP;		///< Number of symmetry operations
	std::string		lattice_name;	///< Name used for the lattice .log file
	std::vector<sym_op<LATFLOAT> >		SYM;			///< Symmetry operations
	Matrix								epsilon_0;		///< Static epsilon (used for screening)
	std::vector<Vector>					RHO;			///< Molecular positions inside the unit cell in orthogonal coordinates
	std::vector<Vector>					MOL_DIPOLE;		///< Molecular dipole moments
	std::vector<std::vector<Vector> >	ATOM;			///< Atomic positions with respect to the CoM
	std::vector<std::vector<LATFLOAT> >	CHARGE;			///< Atomic "transition" charges
 	std::vector<Vector>					CELL;			///< Crystal axes
	std::vector<Vector>					KCELL;			///< Reciprocal lattice vectors
	LATFLOAT							VCELL;		///< Cell volume
	std::string							MODE; 		///< Method used to calculate the resonance-interaction matrix	
											///< Case 'F' -> Finite_sums; case 'S' -> Single_layer_sums; 
											///< 'M'-> Mixed finite charge dist + Ewald; default -> Ewald_sums
	std::string							MOL_MOL_INT;	///< Type of interaction "N"= nearest neighbors; 
														///< "T"=screend distributed transition charge; default=screened_dipole;
	LATFLOAT							NON_SCR_RADIUS;	///< Radius beyond which distributed transition charge interactions and dipole interactions are considered screened
											///< (set it to 0 to always screen interactions and very large to disregard screening)
	LATFLOAT							RMAX;		///< Radius of the sphere (for MODE='F') or circle (for MODE='S') used for finite sum calculations
	Vector								KLIGHT;         ///< k vector of incident light, used to compute bands at the BZ center 
											///< in the 3D grid mode, corresponds to hw = |KLIGHT|*1970
	std::vector<std::vector<LATFLOAT> > NN_int;			///< Nearest neighbours interactions

// Internal variables

	std::vector<std::vector<pos<LATFLOAT> > >	NN;			///< NN contains the nearest neighbours positions (a 3d vector) 
										///< for each type of molecule (1st index) and each neighbour (2nd index)
	Matrix						OtoF;		///< Matrix for transformation from Orthogonal to Fractional coordinates

////////////////////////////////////////////////////////////////////////////////////
// Utilities

	Vector FracToOrth(Vector frac_vec)
	{
		Vector orth_vec;
		if (frac_vec.size()==3)
			orth_vec=get_CELL(0)*frac_vec(0)+get_CELL(1)*frac_vec(1)+get_CELL(2)*frac_vec(2);
		else 
			nrerror("lattice::FracToOrth: frac_vec has the wrong size!");
		return(orth_vec);
	}

	Vector OrthToFrac(Vector orth_vec)
	{
		Vector frac_vec;
		if (orth_vec.size()==3)
			frac_vec=prod(OtoF,orth_vec);
		else 
			nrerror("lattice::OrthToFrac: orth_vec has the wrong size!");
		return(frac_vec);
	}

	Vector move_into_primitive_cell_orth(Vector orth_vec)
	{
		Vector frac_vec, out_orth_vec;
		frac_vec=OrthToFrac(orth_vec);
		for (int i=0; i<3; i++)
				frac_vec(i)-=floor(frac_vec(i));
		out_orth_vec=FracToOrth(frac_vec);
		return(out_orth_vec);
	}
	
	///< Initialises the lattice from input file
	void initialize(std::string input_file);
	
	// CONTROLLARE SCREENING (anche le nuove cose dello screening parziale) !!!!!!!!!!!!!!!!!!!!
	///< Returns the 
	///< Returns energy in eV provided distances are in A and dipoles in Debye
	///< The long-wavelength term is included only if macro_flag!=0
	Complex	Ewald_sums(Vector k, int alpha, int beta, int flag_macro);

	Complex Finite_sums(int flag_2D, Vector k, int alpha, int beta);

    	///< This routine DOES NOT exploit inversion symmetry
	Complex Correction_to_Ewald_sums(Vector k, int alpha, int beta);
	
	///< Computes free exciton bands along given directions
	void compute_FE_bands(std::vector<Vector> dir,  int npoints, int n0,  int n1,  int n2, int flag_macro); 
	
	///< Interaction between transition charge distributions 
    	///< rba is the distance between the center of mass of the two molecules
	LATFLOAT distributed_charge_interaction(int alpha, int beta, Vector rba); 
	
	///< Screened interaction between transition charge distributions 
    	///< rba is the distance between the center of mass of the two molecules
	///< e00, e11, e22 are the three diagonal components of the screening dielectric tensor 
	///< (the non-diagonal components are ignored)
	LATFLOAT screened_distributed_charge_interaction(int alpha, int beta, Vector rba, 
													 LATFLOAT e00, LATFLOAT e11, LATFLOAT e22);

public:

////////////////////////////////////////////////////////////////////////////////////
// Constructors and destructor.

	/// Reads lattice data from file
	lattice(std::string input_file) {initialize(input_file); }
	virtual ~lattice(void) {return;};

////////////////////////////////////////////////////////////////////////////////////
// Functions to read private members (they are all const functions)

	int				get_N_MOL()			const { return(N_MOL); };
	LATFLOAT		get_VCELL()			const { return(VCELL ); };
	int				get_N_SP()			const { return(N_SP); };
	int				get_N_SYM_OP()		const { return(N_SYM_OP); };
	Matrix			get_epsilon_0()		const { return(epsilon_0); };
	std::string 	get_lattice_name()	const { return(lattice_name); };
	std::vector<Vector>	get_KCELL()		const { return(KCELL); };
	std::vector<Vector>	get_CELL()		const { return(CELL); };
	Matrix			get_OtoF()			const { return(OtoF); };
	int				get_NN_size()		const { return(NN.size()); };
	int				get_NN_size(int i)	const { return(NN[i].size()); };
	Vector			get_KLIGHT()		const { return(KLIGHT); };

	sym_op<LATFLOAT> get_SYM(int i)	const
	{ 
		if (i<0 || i>=get_N_SYM_OP())
			nrerror("get_SYM: index out of bounds");
		return(SYM[i]);
	};

	Vector	get_RHO(int i) const
	{ 
		if (i<0 || i>=get_N_MOL())
			nrerror("get_RHO: index out of bounds");
		return(RHO[i]);
	};

	Vector	get_MOL_DIPOLE(int i) const
	{ 
		if (i<0 || i>=get_N_MOL())
			nrerror("get_MOL_DIPOLE: index out of bounds");
		return(MOL_DIPOLE[i]);
	};

	int	get_ATOM_size() const {return(ATOM.size()); };
	int	get_ATOM_size(int i) const {return(ATOM[i].size()); };
	Vector	get_ATOM(int alpha, int i) const
	{ 
		if (alpha<0 || alpha>=get_N_MOL() )
			nrerror("get_ATOM: molecule index out of bounds");
		if (i<0 || i>=get_ATOM_size(alpha) ) 
			nrerror("get_ATOM: atomic index out of bounds");
		return(ATOM[alpha][i]);
	};

	LATFLOAT	get_CHARGE(int alpha, int i) const
	{ 
		if (alpha<0 || alpha>=get_N_MOL() )
			nrerror("get_CHARGE: molecule index out of bounds");
		if (i<0 || i>=get_ATOM_size(alpha) ) 
			nrerror("get_CHARGE: atomic index out of bounds");
		return(CHARGE[alpha][i]);
	};

	Vector	get_CELL(int i) const
	{ 
		if (i<0 || i>=3)
			nrerror("get_CELL: index out of bounds");
		return(CELL[i]);
	};

	Vector	get_KCELL(int i) const
	{ 
		if (i<0 || i>=3)
			nrerror("get_KCELL: index out of bounds");
		return(KCELL[i]);
	};

	std::vector< pos<LATFLOAT> > get_NN(int i) const
	{ 
		if (i<0 || i>=get_N_MOL())
			nrerror("get_NN: index out of bounds");
		return(NN[i]);
	};
	
	pos<LATFLOAT> get_NN(int i, int j) const
	{ 
		if (i<0 || i>=get_N_MOL() || j<0 || j>=get_NN_size(i))
			nrerror("get_NN: index out of bounds");
		return(NN[i][j]);
	};
	
	int get_NN_alpha(int i, int j) const
	{ 
		if (i<0 || i>=get_N_MOL() || j<0 || j>=get_NN_size(i))
			nrerror("get_NN: index out of bounds");
		return(NN[i][j].alpha);
	};
	
	Vector get_NN_n(int i, int j) const
	{ 
		if (i<0 || i>=get_N_MOL() || j<0 || j>=get_NN_size(i))
			nrerror("get_NN: index out of bounds");
		return(NN[i][j].n);
	};

	LATFLOAT	get_epsilon_0(int i, int j) const
	{ 
		if (i<0 || i>=3 || j<0 || j>=3)
			nrerror("get_epsilon_0: index out of bounds");
		return(epsilon_0(i,j));
	};

////////////////////////////////////////////////////////////////////////////////////
// Public functions

	void print(std::ofstream& log_file); ///< Prints information about the lattice to file
	std::vector<std::vector<pos<LATFLOAT> > > compute_NN(LATFLOAT rmax, int flag_2D); ///< Computes the nearest neighbours for each inequivalent molecule
											  ///< Check if it exploits all the symmetry operations 
											  ///< If flag_2D==1 it considers only NN on the same xy plane
											  ///< This function accesses NN without get_ functions
	//	void set_NN(std::vector<std::vector<pos<LATFLOAT> > > new_NN);	
	//	void set_lattice_name(std::string _lattice_name) { if (_lattice_name.length() != 0) lattice_name=_lattice_name; };
	LATFLOAT interaction( int alpha, int beta, Vector dr);	///< Interaction between molecules of type alpha and beta at a distance dr (3D vector)
	LATFLOAT interaction( int alpha, int nn);		///< Interaction between a molecule of type alpha and its nn-th nearest neighboour 

	///< flag_macro is passed to the routine Ewald_sums (if used)
	Complex compute_Ltilde(Vector kvec, int i, int j, int flag_macro);

	void computes_lattice_bands(std::string& band_input_file);

};  // end of lattice class definition

#endif
