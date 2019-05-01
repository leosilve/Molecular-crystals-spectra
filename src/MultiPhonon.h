#ifndef _MULTIPHONON_H
#define _MULTIPHONON_H

#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <math.h>
#include <vector>
#include <map>
#include "../My_include/My_error.h"
#include "../My_boost_numeric/My_boost_N_Functions.h"   
#include "../OP_include/physics.h"
#include "../OP_include/lattice.h"
#include "../OP_include/vector_per.h"   
#include "../OP_include/mpstate.h"

namespace ublas = boost::numeric::ublas;
namespace math = boost::math;
using namespace std;

#ifndef Pi
#define Pi 3.141592654
#endif 
#ifndef ZERO_TOL
#define ZERO_TOL 1e-10
#endif
#ifndef IMAG_TOL
#define IMAG_TOL 1e-5
#endif 

/////////////////////////////////////////////////////////
/// Phonon structure

struct phonon {
	vector_per q;
	int beta;
};

inline bool operator<(phonon& lhs, phonon& rhs)
{
	if (lhs.q < rhs.q)
		return true;
	else
	{
		if (lhs.q==rhs.q && lhs.beta<rhs.beta)
			return true;
		else
			return false;
	}
}

inline bool operator==(phonon& lhs, phonon& rhs)
{
	if (lhs.q==rhs.q && lhs.beta==rhs.beta)
		return true;
	else
		return false;
}

bool ph_comp (phonon& lhs, phonon& rhs)
{
	return (lhs<rhs);
}
///////////////////////////////////////////////////////////////////////////
/////////          Multi Phonon Basis Functions

class MP {
private:
	vector_per k;
	int mol;
	int mu;
	int nph;
	vector<phonon> ph;
public:	
	MP() {
		vector_per k; 
		mol=0;
		mu=0;
		nph=0;
		ph.resize(nph);
	};

	MP(int _nph) {
		vector_per k; 
		mol=0;
		mu=0;
		nph=_nph;
		ph.resize(_nph);
		for (int i=0; i<_nph; i++)
		{
			ph[i].q.set(0,0,0);
			ph[i].beta=0;
		}
	}

	void set_k_per(vector_per param) {k=param;};
	void set_k_per(int _kx, int _ky, int _kz) {k.set(_kx,_ky,_kz);};
	void set_mol(int _mol) {mol=_mol;};
	void set_mu(int _mu) {mu=_mu;};
	void set_nph(int _nph) {
		nph=_nph;
		ph.resize(_nph);
		for (int i=0; i<_nph; i++)
		{
			ph[i].q.set(0,0,0);
			ph[i].beta=0;
		}
	};
	void set_ph(int _i, vector_per _q, int _beta)  {
			ph[_i].q.set_3d(_q);
			ph[_i].beta=_beta;
	}
	
//	vector<phonon> get_raw_ph() {return(ph);};
	int get_mol() {return(mol);};
	int get_mu() {return(mu);};
	int get_nph() {return(nph);}
	phonon get_ph(int i) { return(ph[i]);}
	vector_per get_k() {return(k);};
	void ph_sort() { sort(ph.rbegin(),ph.rend(), ph_comp );}
	int operator== (MP);
	string print_string();
};

ostream &operator<<(ostream& stream, MP mp)   // Screen output of a molecular exciton 
{
	stream << "|" << mp.get_k() << "," << mp.get_mol() << "," << mp.get_mu() << ";" ;
	for (int i=0; i<mp.get_nph(); i++)
			stream    << mp.get_ph(i).q << "," << mp.get_ph(i).beta  << ";";
	stream << "> (nph=" << mp.get_nph() << ") ";
	return stream;
};


//////////////////////////////////////////////////////////////////////////////////////////
//         MP Member functions

string MP::print_string()
{
	ostringstream ss;
	ss << k << mol << mu;
	for (int i=0; i<nph; i++)
			ss << ph[i].q << ph[i].beta;
	return(ss.str());
}

int MP::operator ==(MP mp)
{
	int i;
	if (k==mp.get_k() && mol==mp.get_mol() && mu==mp.get_mu() && nph==mp.get_nph())
	{
		for (i=0; i<nph; i++) 	
			if (ph[i].q!=mp.get_ph(i).q || ph[i].beta!=mp.get_ph(i).beta )
				return(false);
		return(true);
	}
	return(false);
}

///////////////////////////////////////////////////////////////////////////
/////////          Multi Phonon Basis Class

template<class MPFLOAT>
class MultiPhonon 
{

	typedef std::complex<MPFLOAT>							complex;
	typedef ublas::vector<MPFLOAT>							Vector;
	typedef ublas::vector<complex>							ComplexVector;
	typedef ublas::matrix<MPFLOAT,	ublas::column_major>	Matrix;
	typedef ublas::matrix<complex,	ublas::column_major>	ComplexMatrix;

private:

// Problem parameters

	lattice<MPFLOAT>*	latticePtr;	// pointer to the crystal lattice class

	MPFLOAT		lambda_0;		// Huang-Rhys factor (for each species)
	MPFLOAT		lambda_1;		// Huang-Rhys factor high energy mode (for each species)
	MPFLOAT		elec_En;		// Elec transition energy  (for each species)
	MPFLOAT		vib_En_0;		// Vibration energy (for each species)
	MPFLOAT		vib_En_1;		// Vibration high energy mode (for each species)

	int					MAX_PH;			// Max number of phonons per state 
	int					MAX_HI_VIB;		// Max number of high energy vibrations
//	int					N_POL_LEV;		// Number of polarons included
	vector<int>			POL_LEVS;		// List of polaronic levels
	vector<int>			N_STATES;		// Total number of basis states for each assigned number of phonons
	int					N1;				// Used in vector_per
	int					N2;				// Used in vector_per
	int					N3;				// Used in vector_per
	vector_per			kprob;			// Wave vector of the subspace studied
	string				prob_name;		// Name of the MultiPhonon problem
	string				H_print;		// If == "N" or "n" H matrix is not printed
	vector<MPFLOAT>		EN_MAX;			// Maximum energy of the BD states included in the spectrum for each phonon n subspace
	vector<MPFLOAT>		Sn0_1;			// Franck Condon factors stored for quick usage
	int					N_MOL;			// n. of molecules in the crystal unit cell
	string				MODE;
	string				mol_mol_int;
	MPFLOAT				NON_SCR_RADIUS;
	MPFLOAT				RMAX;
	Vector				KLIGHT;
	vector<vector<MPFLOAT> > NN_int;

// Problem flags 

	int					EXACT_Q0;		// If ==1 then the q=0 phonon is treated exactly
	int					FLAT_BAND;		// If ==1 then a flat band is assumed
/*	int					SUPER_BASIS;	// If ==1 and BD FLAG==1 then uses the approximate basis set of my model
	int					HBASIS;			// If ==1 computes Hbasis and uses it
	int					CHECK_BD;		// If ==1 checks BD_BASIS set
	int					NEW_METHOD;		// If ==1 uses the new method for polaronic states*/

// Stored variables

	vector<vector<vector<ComplexMatrix > > >	Ltilde_all;			// Ltilde matrix elements for each k
	vector<vector<vector<ComplexMatrix > > >	udc_all;			// Ltilde matrix eigenvectors for each k
	ComplexMatrix								udc_kprob;				// Ltilde matrix eigenvectors at kprob
	vector<MP>									BASIS;				// Basis states
	vector<vector<MPSTATE<MPFLOAT>> >			BD_BASIS;			// Block diagonal basis set

public:

	int get_BASIS_size() {return(BASIS.size());}
	int get_BD_BASIS_size(int i) 
	{
		if (i>=0 && i<BD_BASIS.size()) 
			return(BD_BASIS[i].size());
		else
			nrerror("MultiPhonon::get_BD_BASIS_size: index out of bounds!");
	}
	lattice<MPFLOAT>*	get_latticePtr() {return(latticePtr);}

//	Constructors 

	// Short Constructor
	MultiPhonon() {};
	// Long Constructor
	MultiPhonon(string Problem_name, string input_file, lattice<MPFLOAT>* _latticePtr) { latticePtr=_latticePtr; initialize(Problem_name, input_file); }

	MPFLOAT			get_lambda_0() { return(lambda_0); };
	MPFLOAT 		get_lambda_1() { return(lambda_1); };
	MPFLOAT 		get_elec_En() { return(elec_En); };
	MPFLOAT 		get_vib_En_0() { return(vib_En_0); };
	MPFLOAT 		get_vib_En_1() { return(vib_En_1); };
//	MPFLOAT 		get_epsilon_inf() { return(epsilon_inf); };
//	MPFLOAT 		get_gamma_coeff() { return(gamma_coeff); };
	void			set_lambda_1(MPFLOAT _lambda_1) { lambda_1=_lambda_1; };

// Reads input data
	void initialize(string Problem_name, string input_file)
	{
//		prob_name=input_file.substr(0,input_file.find('.'));
		prob_name=Problem_name;// + "_MultiPhonon";

		cout << "Reading MP problem data from file: " << input_file << endl;

		int i,j,k,n;
		string line;
		istringstream  iss;
		ifstream infile(input_file.c_str());
		if (!infile.is_open()) nrerror("File not found!");

		MPFLOAT		d_temp;
		int			i_temp;

		getline(infile, line, ':');
		infile >> lambda_0;
		getline(infile, line, ':');
		infile >> lambda_1;
		getline(infile, line, ':');
		infile >> elec_En;
		getline(infile, line, ':');
		infile >> vib_En_0;
		getline(infile, line, ':');
		infile >> vib_En_1;
		getline(infile, line, ':');
		getline(infile, line, '"');
		getline(infile, MODE, '"');	
		getline(infile, line, ':');
		getline(infile, line, '"');
		getline(infile, mol_mol_int, '"');	
		getline(infile, line, ':');
		infile >> NON_SCR_RADIUS;
		getline(infile, line, ':');
		infile >> RMAX;
		getline(infile, line, ':');
		infile >> i_temp;
		if (i_temp!=(*latticePtr).get_N_MOL())
		{
			cout << "No valid NN_interactions, setting NN_int size to 0!" << endl;
			NN_int.resize(0);
			getline(infile, line, ':');
			getline(infile, line, ':');
		}
		else
		{
			NN_int.resize((*latticePtr).get_N_MOL());
			getline(infile, line, ':');
			infile >> i_temp;
			for (i=0; i<(*latticePtr).get_N_MOL(); i++)
				NN_int[i].resize(i_temp);
			getline(infile, line, ':');
			for (i=0; i<(*latticePtr).get_N_MOL(); i++)
				for (j=0; j<NN_int[i].size(); j++)
					infile >> NN_int[i][j];
		}	
		getline(infile, line, ':');
		KLIGHT.resize(3);
		for (i=0;i<3;i++) 
			infile >> KLIGHT[i];
		getline(infile, line, ':');
		infile >> N1;	infile >> N2;	infile >> N3;
		vector_per::N1=N1;
		vector_per::N2=N2;
		vector_per::N3=N3;
		getline(infile, line, ':');
		infile >> i_temp; kprob.set_n1(i_temp);	
		infile >> i_temp; kprob.set_n2(i_temp);	
		infile >> i_temp; kprob.set_n3(i_temp);
		getline(infile, line, ':');
		infile >> MAX_PH;
		if (MAX_PH<0)
		{
			MAX_PH=0;
			lambda_0=0.0;
		}
		getline(infile, line, ':');
		infile >> MAX_HI_VIB;
		if (MAX_HI_VIB<0)
		{
			MAX_HI_VIB=0;
			set_lambda_1(0.0);
		}	
		
		getline(infile, line, ':');
		getline(infile, line);
		iss.str( line );
		copy( istream_iterator<int>( iss ), istream_iterator<int>(), back_inserter( POL_LEVS ) );
		if (POL_LEVS.size()>MAX_HI_VIB+1)
			POL_LEVS.resize(MAX_HI_VIB+1);
		getline(infile, line, ':');
		EN_MAX.resize(4);
		for (i=0; i<4; i++)
			infile >> EN_MAX[i];
		getline(infile, line, ':');
		EXACT_Q0=0;
		FLAT_BAND=0;
/*		infile >> EXACT_Q0;
		if (EXACT_Q0==1)
			cout << "WARNING: check computation of state dipoles for EXACT_Q0==1 !!!!!!!" << endl;
		else
			EXACT_Q0=0;*/
/*		getline(infile, line, ':');
		infile >> FLAT_BAND;*/
/*		getline(infile, line, ':');
		infile >> SUPER_BASIS;
		getline(infile, line, ':');
		infile >> HBASIS;
		getline(infile, line, ':');
		infile >> CHECK_BD;
		getline(infile, line, ':');
		infile >> NEW_METHOD;*/
		infile.close();

		H_print="N";
/*		ostringstream ss;
		if (MAX_HI_VIB>=0)
		{
			ss << prob_name << "_N" << N1 << "_P" << MAX_PH  
			//<< "_k" << k.n1() << k.n2() << k.n3() 
			<< "_V" << MAX_HI_VIB << "_L";
			for (i=0; i<POL_LEVS.size(); i++)
				ss << POL_LEVS[i];
		}
		else
		{
			ss << prob_name << "_N" << N1 << "_P" << MAX_PH  
			//<< "_k" << k.n1() << k.n2() << k.n3() 
			<< "_NV" ;
		}	
		prob_name = ss.str();*/
		Sn0_1.resize(MAX_HI_VIB+1);
		for (i=0; i<MAX_HI_VIB+1; i++)
			Sn0_1[i]=S_factor(i, 0, lambda_1);
		N_MOL=(*latticePtr).get_N_MOL();

// This routine copies run_data input file
/*		ifstream infile("Run_data_MP.txt", std::ios_base::binary);
		string run_filename = prob_name + "_run_data.txt";
		ofstream outfile(run_filename.c_str(), std::ios_base::binary);
		outfile << infile.rdbuf();
		outfile.close();*/

		int n1, n2, n3;
		ComplexMatrix zero_mat(N_MOL,N_MOL);
		for (i=0; i<N_MOL; i++)	
			for (j=0; j<N_MOL; j++)
				zero_mat(i,j)=0.0;

		Ltilde_all.resize(N1);
		for (n1=0; n1<N1; n1++) 
			Ltilde_all[n1].resize(N2);
		for (n1=0; n1<N1; n1++)
			for (n2=0; n2<N2; n2++)
				Ltilde_all[n1][n2].resize(N3);
		for (n1=0; n1<N1; n1++)
			for (n2=0; n2<N2; n2++)
				for (n3=0; n3<N3; n3++)
				{
					Ltilde_all[n1][n2][n3].resize(N_MOL,N_MOL);
					Ltilde_all[n1][n2][n3]=zero_mat;
				}
	
		udc_all.resize(N1);
		for (n1=0; n1<N1; n1++) 
			udc_all[n1].resize(N2);
		for (n1=0; n1<N1; n1++)
			for (n2=0; n2<N2; n2++)
				udc_all[n1][n2].resize(N3);
		for (n1=0; n1<N1; n1++)
			for (n2=0; n2<N2; n2++)
				for (n3=0; n3<N3; n3++)
				{
					udc_all[n1][n2][n3].resize(N_MOL,N_MOL);
					udc_all[n1][n2][n3]=zero_mat;
				}

//		eigenvalues.resize(N_MOL);
//		eigenvectors.resize(N_MOL);
//		dipoles.resize(N_MOL);

/*		string log_filename = prob_name + ".log";
		ofstream outfile(log_filename.c_str());
		print(outfile);
		outfile.close();*/
	
		compute_free_exciton_dispersion();
/*	
		/*		pol_eigenvec.resize(N1*N2*N3*N_MOL);
		for (i=0; i<N1*N2*N3*N_MOL; i++) 
			pol_eigenvec[i].resize((MAX_HI_VIB+1)*(MAX_HI_VIB+1));
	
		pol_eigenval.resize(N1*N2*N3*N_MOL);
		for (i=0; i<N1*N2*N3*N_MOL; i++) 
			pol_eigenval[i].resize(MAX_HI_VIB+1);

		if (N_POL_LEV>0)
		{
			vector_per prob_k=k;
			for (n1=0; n1<N1; n1++)	for (n2=0; n2<N2; n2++)	for (n3=0; n3<N3; n3++)
			{	
				k.set(n1,n2,n3);			
				compute_polaronic_states(k);
			}
			k=prob_k;
			udc_k = get_udc(k).copy();
		}
*/
		compute_BASIS(kprob); 
		return;
	}

void print(ofstream& log)
{
		cout << "Printing MultiPhonon Basis data ... " << endl;
		int N=N1*N2*N3;
//		string log_filename = prob_name + ".log";
//		ofstream log(log_filename.c_str());
		log << endl;
		log << "==================================================================================================" << endl;
		log << "MultiPhonon problem: " << prob_name << endl;
		log << "==================================================================================================" << endl << endl;
//		log << "Lattice : " << (*latticePtr).get_lattice_name() << endl;
		log << "lambda_0: " << get_lambda_0() << endl;
		log << "lambda_1: " << get_lambda_1() << endl;
		log << "elec_En: " << get_elec_En() << endl;
		log << "vib_En_0: " << get_vib_En_0() << endl;
		log << "vib_En_1: " << get_vib_En_1() << endl;
		log << "MODE: " <<	MODE << endl;
		log << "mol_mol_int: " <<	mol_mol_int << endl;
		log << "NON_SCR_RADIUS: " <<	NON_SCR_RADIUS << endl;
		log << "RMAX: " <<	RMAX << endl;
		log << "KLIGHT: " <<	KLIGHT << endl;
		if (NN_int.size()>0)
			log << "NN_int: " <<	NN_int << endl;
		log << "N1 : " << N1 << endl;
		log << "N2 : " << N2 << endl;
		log << "N3 : " << N3 << endl;
		log << "MAX_PH : " << MAX_PH << endl;
		log << "EN_MAX : " << EN_MAX << endl;
		log << "MAX_HI_VIB : " << MAX_HI_VIB << endl;
		if (POL_LEVS.size()>0)
			log << "POL_LEVS: " << endl << POL_LEVS << endl;
/*		log << "EXACT_Q0 : " << EXACT_Q0 << endl;
		log << "FLAT_BAND : " << FLAT_BAND << endl;
		log << "HBASIS : " << HBASIS << endl;
		log << "CHECK_BD : " << CHECK_BD << endl;
		log << "NEW_METHOD : " << NEW_METHOD << endl;*/
		log << "k : " << kprob << endl;
		log << "N. of states with 0 phonons: " << N_MOL*(MAX_HI_VIB+1) << endl;
		if (MAX_PH>0)
			log << "N. of states with 1 phonons: " << N*(MAX_HI_VIB+1)*N_MOL*N_MOL << endl;
		if (MAX_PH>1)
			log << "N. of states with 2 phonons: " << N_MOL*(N*leo::CRnk(N_MOL,2)+math::binomial_coefficient<MPFLOAT>(N,2)*N_MOL*N_MOL)*(MAX_HI_VIB+1) << endl;
		if (MAX_PH>2)
			log << "N. of states with 3 phonons: " << N_MOL*(N*leo::CRnk(N_MOL,3)+leo::Dnk(N,2)*leo::CRnk(N_MOL,2)*N_MOL+math::binomial_coefficient<MPFLOAT>(N,3)*N_MOL*N_MOL*N_MOL)*(MAX_HI_VIB+1) << endl;
		log << endl;
		log << endl << "------------------------------------------------------"  << endl;
		log << "MP dispersion " << endl << endl;
		int n1,n2,n3;
		vector_per		kper;
		vector<MPFLOAT> Values;
		ComplexMatrix udc;
		for (n1=0; n1<N1; n1++)		
		for (n2=0; n2<N2; n2++)
		for (n3=0; n3<N3; n3++)
		{
			kper.set(n1,n2,n3);
//			cout << Ltilde_all[n1][n2][n3] << endl;
			udc=Ltilde_all[n1][n2][n3];
			heev_cpp('V', 'U', udc, Values);
			log << "k = " << kper  << endl << "Ltilde = " << Ltilde_all[n1][n2][n3] << endl 
				<< "eigenvectors = " << udc_all[n1][n2][n3] << endl 
				<< "eigenvalues = " << endl << Values << endl;
//			cout << Ltilde_all[n1][n2][n3] << endl;
//			getchar();
		}

	log << endl << "---------------------------------------" << endl;
	log << " BASIS" << endl << endl;
	log << "Number of states : " << BASIS.size() << endl;
	log << N_STATES << endl << endl;
	for (int i=0; i<BASIS.size(); i++)
		log << "State n." << i << " : " << BASIS[i] << " En: " << Hint(BASIS[i],BASIS[i]) << endl;
	return;
}

ComplexVector get_BASIS_dipole(int i)
{
	return(get_dipole(BASIS[i]));
}
ComplexVector get_BD_BASIS_dipole(int dc, int i)
{
	return(get_dipole(BD_BASIS[dc][i]));
}
ComplexVector get_dipole(MP& mp)
{
	complex I(0.0,1.0);
	ComplexVector  dipole(3);
	dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
/*	if ((*latticePtr).get_MODE()=="Spano")
// This function is to compare with Spano's papers
// Check what happens for EXACT_Q0==1 !!!!!!!!!!!!!!
	{
		vector_3d<double> bpol(0.0,0.025,0.0);
		vector_3d<double> acpol(0.0433,0.0,0.9987);

		if (EXACT_Q0!=1)
		{
			if (mp.get_k().n1()==0 && mp.get_k().n2()==0 && mp.get_k().n3()==0 && mp.get_nph()==0)
				return( bpol * S_factor(mp.get_mu(), 0, lambda_1 ) );

			if (mp.get_k().n1()==N1/2 && mp.get_k().n2()==N2/2 && mp.get_k().n3()==N3/2 && mp.get_nph()==0
				&& N1%2==0 && N2%2==0 )
				return( acpol * S_factor(mp.get_mu(), 0, lambda_1 ) );
		}
		else
		{
			if (mp.get_k().n1()==0 && mp.get_k().n2()==0 && mp.get_k().n3()==0 )
			{
				int allq0_flag=1;
				for (int i=0; i<mp.get_nph(); i++)
					if (mp.get_ph(i).q.n1()!=0 || mp.get_ph(i).q.n2()!=0 || mp.get_ph(i).q.n3()!=0) 
					{
						allq0_flag=0;
						break;
					}
				if (allq0_flag==1)
				{
					return( bpol * S_factor(mp.get_mu(), 0, lambda_1 ) 
							     * S_factor(mp.get_nph(), 0, lambda_0/sqrt(MPFLOAT(N1*N2*N3)) ));
				}
			}

			if (mp.get_k().n1()==N1/2 && mp.get_k().n2()==N2/2 && mp.get_k().n3()==N3/2 
				&& N1%2==0 && N2%2==0 )
			{
				int allq0_flag=1;
				for (int i=0; i<mp.get_nph(); i++)
					if (mp.get_ph(i).q.n1()!=0 || mp.get_ph(i).q.n2()!=0 || mp.get_ph(i).q.n3()!=0) 
					{
						allq0_flag=0;
						break;
					}
				if (allq0_flag==1)
				return( acpol * S_factor(mp.get_mu(), 0, lambda_1 ) 
							  * S_factor(mp.get_nph(), 0, lambda_0/sqrt(MPFLOAT(N1*N2*N3)) ));
			}
		}

		return( (*latticePtr).get_MOL_DIPOLE(mp.get_mol()) * 0.0);
	}
	else
// This function is the correct one for my model
// Check what happens for EXACT_Q0==1 !!!!!!!!!!!!!!
	{
		if (EXACT_Q0!=1)
		{*/

			if (mp.get_k().n1()==0 && mp.get_k().n2()==0 && mp.get_k().n3()==0 && mp.get_nph()==0)
//				dipole=(*latticePtr).get_MOL_DIPOLE(mp.get_mol())*exp(-I*inner_prod(get_k_3d(kprob),(*latticePtr).get_RHO(mp.get_mol())))* S_factor(mp.get_mu(), 0, lambda_1 );
				dipole=(*latticePtr).get_MOL_DIPOLE(mp.get_mol())*S_factor(mp.get_mu(), 0, lambda_1 );
			return(dipole);
/*		}
		else
		{
			if (N_MOL==1)
			{
			if (mp.get_k().n1()==0 && mp.get_k().n2()==0 && mp.get_k().n3()==0)
			{
				int i;
				int allq0_flag=1;
				for (i=0; i<mp.get_nph(); i++)
					if (mp.get_ph(i).q.n1()!=0 || mp.get_ph(i).q.n2()!=0 || mp.get_ph(i).q.n3()!=0 
						|| mp.get_ph(i).beta!=mp.get_mol() ) 
					{
						allq0_flag=0;
						break;
					}
				if (allq0_flag==1)
					return( (*latticePtr).get_MOL_DIPOLE(mp.get_mol()) * S_factor(mp.get_mu(), 0, lambda_1 ) 
															  * S_factor(mp.get_nph(), 0, lambda_0/sqrt(MPFLOAT(N1*N2*N3)) ));
			}
			return( (*latticePtr).get_MOL_DIPOLE(mp.get_mol()) * 0.0);
			}
			else
			{
				nrerror("get_dipole: Routine works only with 1 molecule per cell!");
			}
		}
	}*/
}

ComplexVector get_dipole(MPSTATE<MPFLOAT>& mps)
{
	int i;
	ComplexVector  dipole(3);
	dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
	for (i=0; i<mps.size(); i++)
		dipole += ( get_dipole(BASIS[mps.get_state(i)]) * mps.get_coeff(i) );
	return(dipole);
}

Vector get_k_3d(vector_per &kper)
{
	return((*latticePtr).get_KCELL(0)*(kper.n1()*1.0/N1)+(*latticePtr).get_KCELL(1)*(kper.n2()*1.0/N2)
		   +(*latticePtr).get_KCELL(2)*(kper.n3()*1.0/N3));
}

void compute_free_exciton_dispersion()
{
	cout << "Computing free exciton dispersion" << endl;
	int i,j;
	complex			I(0.0,1.0);
	ComplexMatrix	Ltilde_temp(N_MOL,N_MOL);
	ComplexMatrix	udc_temp(N_MOL,N_MOL);
	vector<MPFLOAT>	values(N_MOL);
	vector_per		kper, ksym;
	Vector			k3d;
	int				flag_macro=0;
//	string			log_filename = prob_name + ".log";
//	ofstream		outfile(log_filename.c_str(), ios::app);
//	string			band_filename = prob_name + "_dispersion.txt";
//	ofstream		bandfile(band_filename.c_str());
//	outfile << endl << "------------------------------------------------------"  << endl;
//	outfile << "MP dispersion " << endl << endl;

	int dc, n1, n2, n3;
	if (FLAT_BAND==1)
	{
		nrerror("MultiPhonon:: compute_free_exciton_dispersion: Flat Band dispersion not available!");
	}
	else
	{
		for (n1=0; n1<N1; n1++)		
		for (n2=0; n2<N2; n2++)
		for (n3=0; n3<N3; n3++)
		{	
			kper.set(n1,n2,n3);	
			k3d=get_k_3d(kper);
//			if (norm_2(k3d)<ZERO_TOL)	
//				k3d+=KLIGHT;
			(*latticePtr).set_Ltilde_and_u(k3d, MODE, mol_mol_int, NON_SCR_RADIUS, NN_int, RMAX, 
			                  flag_macro, Ltilde_temp, udc_temp, values);
//			outfile << kper  << endl << Ltilde_temp << endl << udc_temp << endl << values << endl;
/*			bandfile << kper.n1()  << "\t" << kper.n2()  << "\t" << kper.n3();
			for (dc=0; dc<N_MOL; dc++)
				bandfile << "\t" << values[dc];
			bandfile << endl;*/
			Ltilde_all[n1][n2][n3]=Ltilde_temp;
			udc_all[n1][n2][n3]=udc_temp;
			if (kper==kprob) 
				udc_kprob = udc_temp;	
//	This routine was probably added to have exactly the same eigenvalues for equivalent k points
/*			if (N_POL_LEV>0)
			{
				ksym=get_sym_k(1, kper);
				TNT::Array2D<complex<double> u_sym(N_MOL,N_MOL);
				u_sym[0][0]=u_temp[1][0];
				u_sym[1][0]=u_temp[0][0];
				u_sym[0][1]=u_temp[1][1];
				u_sym[1][1]=u_temp[0][1];
				udc[vector_per_mod::posmod(ksym.n1(),N1)][vector_per_mod::posmod(ksym.n2(),N2)][vector_per_mod::posmod(ksym.n3(),N3)]=u_sym.copy();
				u_sym[0][0]=e_temp[1][1];
				u_sym[1][0]=e_temp[0][1];
				u_sym[0][1]=e_temp[1][0];
				u_sym[1][1]=e_temp[0][0];
				Lmat[vector_per_mod::posmod(ksym.n1(),N1)][vector_per_mod::posmod(ksym.n2(),N2)][vector_per_mod::posmod(ksym.n3(),N3)]=u_sym.copy();
				if (k==kper) 
					udc_k = udc[vector_per_mod::posmod(k.n1(),N1)][vector_per_mod::posmod(k.n2(),N2)][vector_per_mod::posmod(k.n3(),N3)].copy();
			}*/
		}	
	}
//    outfile << endl << "Davydov components matrix at kprob (each column is a Davydov state) " << endl;
//	outfile << udc_kprob << endl;
//	outfile.close();
//	bandfile.close();
	return;
}

int get_sym_mol(int sym, int i)
{
	return ((*latticePtr).get_SYM(sym).get_map(i));
}

// This routine returns the wave vector obtained by operating on k with a symmetry operation of the factor group
// funziona solo se la matrice di rotazione e' diagonale !!!!!!!!!!
vector_per get_sym_k(int sym, vector_per k)
{
	vector_per temp;
	temp.set_n1(k.n1()*(*latticePtr).get_SYM(sym).get_rot()(0,0));
	temp.set_n2(k.n2()*(*latticePtr).get_SYM(sym).get_rot()(1,1));
	temp.set_n3(k.n3()*(*latticePtr).get_SYM(sym).get_rot()(2,2));
	return(temp);
}

// Check if it is ok for k!=0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void compute_BASIS(vector_per k)
{
/*	if ((*latticePtr).get_N_SP()>1)
		nrerror("compute_BASIS: cannot compute BASIS for N_SP>1!");
	if ((k.n1()!=0 || k.n2()!=0 || k.n3()!=0) && N_MOL>1)
		nrerror("compute_BASIS: cannot compute BASIS for k!=0 and N_MOL>1!"); // Perche'?*/
	cout << "Computing BASIS" << endl;

// -----------------------------------Computes the total number of BASIS states
	N_STATES.resize(MAX_PH+1);
	int N=N1*N2*N3;
	int nstates=N_MOL*(MAX_HI_VIB+1);		// n. di stati con 0 fononi
	N_STATES[0]=N_MOL*(MAX_HI_VIB+1);

	if (MAX_PH>0)
	{
		nstates+=N*(MAX_HI_VIB+1)*N_MOL*N_MOL;  // n. di stati con 1 fonone
		N_STATES[1]=N*N_MOL*N_MOL*(MAX_HI_VIB+1);
	}

	if (MAX_PH>1)
	{
		nstates+=N_MOL*(N*leo::CRnk(N_MOL,2)+math::binomial_coefficient<MPFLOAT>(N,2)*N_MOL*N_MOL)*(MAX_HI_VIB+1);  // n. di stati con 2 fononi
		N_STATES[2]=N_MOL*(N*leo::CRnk(N_MOL,2)+math::binomial_coefficient<MPFLOAT>(N,2)*N_MOL*N_MOL)*(MAX_HI_VIB+1);
	}
	if (MAX_PH>2)
	{
		nstates+=N_MOL*(N*leo::CRnk(N_MOL,3)+leo::Dnk(N,2)*leo::CRnk(N_MOL,2)*N_MOL+math::binomial_coefficient<MPFLOAT>(N,3)*N_MOL*N_MOL*N_MOL)*(MAX_HI_VIB+1);  // n. di stati con 3 fononi
		N_STATES[3]=N_MOL*(N*leo::CRnk(N_MOL,3)+leo::Dnk(N,2)*leo::CRnk(N_MOL,2)*N_MOL+math::binomial_coefficient<MPFLOAT>(N,3)*N_MOL*N_MOL*N_MOL)*(MAX_HI_VIB+1);
	}

	cout << N_STATES << endl;
	int nmol,nsp;
	int n1k, n2k, n3k;
	int beta, beta1, beta2, beta3;
	int	n1q1, n2q1, n3q1, n1q2, n2q2, n3q2, n1q3, n2q3, n3q3;
	int sym, tempcount;
	int np;
	vector_per q1,q2,q3;
	MP mp_sym;

// ------------------------------  Finds the BASIS states
	BASIS.resize(nstates);
	int count=0;
	int nm=N_MOL/(*latticePtr).get_N_SP();

	for (nsp=0; nsp<(*latticePtr).get_N_SP(); nsp++)
	{
		nmol=nsp*nm;

// 0 phonon states (they have all the same k_per=k)

		BASIS[count].set_k_per(k);
		BASIS[count].set_mol(nmol);
		BASIS[count].set_mu(0);
		BASIS[count].set_nph(0);
		count++;
		tempcount=count-1;
		for (sym=1; sym<nm; sym++)
		{
			mp_sym.set_k_per(k);
			mp_sym.set_mol(get_sym_mol(sym, BASIS[tempcount].get_mol()));
			mp_sym.set_mu(0);
			mp_sym.set_nph(0);
			BASIS[count]=mp_sym;
			count++;
		}

// 1 phonon states
		if (MAX_PH>0)
		{
			for (beta=0; beta<N_MOL; beta++)
			for (n1k=0; n1k<N1; n1k++)		
				for (n2k=0; n2k<N2; n2k++)
					for (n3k=0; n3k<N3; n3k++)
			{
				BASIS[count].set_k_per(n1k,n2k,n3k);
				BASIS[count].set_mol(nmol);
				BASIS[count].set_mu(0);
				BASIS[count].set_nph(1);
				BASIS[count].set_ph(0, k-BASIS[count].get_k(), beta);
				count++;
				tempcount=count-1;
				for (sym=1; sym<nm; sym++)
				{
					compute_symmetric_state(BASIS[tempcount], mp_sym, sym);
					BASIS[count]=mp_sym;
					count++;
				}
			}
		}

// 2 phonons states
		if (MAX_PH>1)
		{
// These are the states with q1=q2 
			for (beta1=0; beta1<N_MOL; beta1++)
				for (beta2=beta1; beta2<N_MOL; beta2++)
			for (n1q1=0; n1q1<N1; n1q1++)		
				for (n2q1=0; n2q1<N2; n2q1++)
					for (n3q1=0; n3q1<N3; n3q1++)
			{
				q1.set(n1q1,n2q1,n3q1);
				q2.set_3d(q1);
				BASIS[count].set_k_per(k-q1-q2);
				BASIS[count].set_mol(nmol);
				BASIS[count].set_mu(0);
				BASIS[count].set_nph(2);
				BASIS[count].set_ph(0, q1, beta1);
				BASIS[count].set_ph(1, q2, beta2);
				count++;
				tempcount=count-1;
				for (sym=1; sym<nm; sym++)
				{
					compute_symmetric_state(BASIS[tempcount], mp_sym, sym);
					BASIS[count]=mp_sym;
					count++;
				}
			}

// These are the states with q1>q2 
			for (beta1=0; beta1<N_MOL; beta1++)
						for (beta2=0; beta2<N_MOL; beta2++)
			for (n1q1=0; n1q1<N1; n1q1++)		
				for (n2q1=0; n2q1<N2; n2q1++)
					for (n3q1=0; n3q1<N3; n3q1++)
			{
				q1.set(n1q1,n2q1,n3q1);
				for (n1q2=0; n1q2<N1; n1q2++)		
					for (n2q2=0; n2q2<N2; n2q2++)
						for (n3q2=0; n3q2<N3; n3q2++)
				{
					q2.set(n1q2,n2q2,n3q2);
					if (q1>q2)
					{
						BASIS[count].set_k_per(k-q1-q2);
						BASIS[count].set_mol(nmol);
						BASIS[count].set_mu(0);
						BASIS[count].set_nph(2);
						BASIS[count].set_ph(0, q1, beta1);
						BASIS[count].set_ph(1, q2, beta2);
						count++;	
						tempcount=count-1;
						for (sym=1; sym<nm; sym++)
						{
							compute_symmetric_state(BASIS[tempcount], mp_sym, sym);
							BASIS[count]=mp_sym;
							count++;
						}
					}
				}
			}
		}   // end if MAX_PH>1

// 3 phonons states
		if (MAX_PH>2)
		{
// These are the states with two equal vectors q1=q2 or q2=q3 
			for (beta1=0; beta1<N_MOL; beta1++)
				for (beta2=beta1; beta2<N_MOL; beta2++)
			for (n1q1=0; n1q1<N1; n1q1++)		
				for (n2q1=0; n2q1<N2; n2q1++)
					for (n3q1=0; n3q1<N3; n3q1++)
			{
				q1.set(n1q1,n2q1,n3q1);
				q2.set_3d(q1);
				for (beta3=0; beta3<N_MOL; beta3++)
				for (n1q3=0; n1q3<N1; n1q3++)		
					for (n2q3=0; n2q3<N2; n2q3++)
						for (n3q3=0; n3q3<N3; n3q3++)
						{
							q3.set(n1q3,n2q3,n3q3);
							if (q1>q3)
							{
								BASIS[count].set_k_per(k-q1-q2-q3);
								BASIS[count].set_mol(nmol);
								BASIS[count].set_mu(0);
								BASIS[count].set_nph(3);
								BASIS[count].set_ph(0, q1, beta1);
								BASIS[count].set_ph(1, q2, beta2);
								BASIS[count].set_ph(2, q3, beta3);
								count++;
								tempcount=count-1;
								for (sym=1; sym<nm; sym++)
								{
									compute_symmetric_state(BASIS[tempcount], mp_sym, sym);
									BASIS[count]=mp_sym;
									count++;
								}
							}
							if (q1<q3)
							{
								BASIS[count].set_k_per(k-q1-q2-q3);
								BASIS[count].set_mol(nmol);
								BASIS[count].set_mu(0);
								BASIS[count].set_nph(3);
								BASIS[count].set_ph(0, q3, beta3);
								BASIS[count].set_ph(1, q1, beta1);
								BASIS[count].set_ph(2, q2, beta2);
								count++;
								tempcount=count-1;
								for (sym=1; sym<nm; sym++)
								{
									compute_symmetric_state(BASIS[tempcount], mp_sym, sym);
									BASIS[count]=mp_sym;
									count++;
								}
							}
						}
			}

// These are the states with three equal vectors q1=q2=q3 

			for (n1q1=0; n1q1<N1; n1q1++)		
				for (n2q1=0; n2q1<N2; n2q1++)
					for (n3q1=0; n3q1<N3; n3q1++)
			{
				q1.set(n1q1,n2q1,n3q1);
				q2.set_3d(q1);
				q3.set_3d(q1);
				for (beta1=0; beta1<N_MOL; beta1++)
					for (beta2=beta1; beta2<N_MOL; beta2++)
						for (beta3=beta2; beta3<N_MOL; beta3++)
						{
								BASIS[count].set_k_per(k-q1-q2-q3);
								BASIS[count].set_mol(nmol);
								BASIS[count].set_mu(0);
								BASIS[count].set_nph(3);
								BASIS[count].set_ph(0, q1, beta1);
								BASIS[count].set_ph(1, q2, beta2);
								BASIS[count].set_ph(2, q3, beta3);
								count++;
								tempcount=count-1;
								for (sym=1; sym<nm; sym++)
								{
									compute_symmetric_state(BASIS[tempcount], mp_sym, sym);
									BASIS[count]=mp_sym;
									count++;
								}
						}
			}

// These are the states with q1>q2>q3 
			for (beta1=0; beta1<N_MOL; beta1++)
				for (beta2=0; beta2<N_MOL; beta2++)
					for (beta3=0; beta3<N_MOL; beta3++)
			for (n1q1=0; n1q1<N1; n1q1++)		
				for (n2q1=0; n2q1<N2; n2q1++)
					for (n3q1=0; n3q1<N3; n3q1++)
			{
				q1.set(n1q1,n2q1,n3q1);
				for (n1q2=0; n1q2<N1; n1q2++)		
					for (n2q2=0; n2q2<N2; n2q2++)
						for (n3q2=0; n3q2<N3; n3q2++)
				{
					q2.set(n1q2,n2q2,n3q2);
					for (n1q3=0; n1q3<N1; n1q3++)		
						for (n2q3=0; n2q3<N2; n2q3++)
							for (n3q3=0; n3q3<N3; n3q3++)
					{
						q3.set(n1q3,n2q3,n3q3);
						if (q1>q2 && q2>q3)
						{
							BASIS[count].set_k_per(k-q1-q2-q3);
							BASIS[count].set_mol(nmol);
							BASIS[count].set_mu(0);
							BASIS[count].set_nph(3);
							BASIS[count].set_ph(0, q1, beta1);
							BASIS[count].set_ph(1, q2, beta2);
							BASIS[count].set_ph(2, q3, beta3);
							count++;
							tempcount=count-1;
							for (sym=1; sym<nm; sym++)
							{
								compute_symmetric_state(BASIS[tempcount], mp_sym, sym);
								BASIS[count]=mp_sym;
								count++;
							}
						}
					}
				}
			}
		}   // end if MAX_PH>2

	} // end for (nsp=0; nsp<(*latticePtr).get_N_SP(); nsp++)

// Adding states with mu>0

	tempcount=count;
	int n, n_vib;
	int num_ph;
	for (n_vib=1; n_vib<=MAX_HI_VIB; n_vib++)
		for (n=0; n<tempcount; n++)
		{	
			BASIS[count].set_k_per(BASIS[n].get_k());
			BASIS[count].set_mol(BASIS[n].get_mol());
			BASIS[count].set_mu(n_vib);
			num_ph=BASIS[n].get_nph();
			BASIS[count].set_nph(num_ph);
			for (np=0; np<num_ph; np++)
			{			
				BASIS[count].set_ph(np, BASIS[n].get_ph(np).q, BASIS[n].get_ph(np).beta);
			}
			count++;
		}
	tempcount=count;

	int n_comp_states=count;
	if (n_comp_states!=nstates)
		nrerror("compute_BASIS: number of computed states is different from prediction!");

//  Sorting all BASIS states phonon part

	cout << "Sorting BASIS states phonon part ..." << endl;
	int i;
	for (i=0; i<n_comp_states; i++)
		BASIS[i].ph_sort();

//  Writing BASIS to logfile
	
/*	cout << "Writing to logfile ..." << endl;
	string log_filename = prob_name + ".log";
	ofstream outfile(log_filename.c_str(), ios::app);
	outfile << "---------------------------------------" << endl;
	outfile << " BASIS" << endl << endl;
	outfile << "Number of states : " << count << endl;
	outfile << N_STATES << endl << endl;
	for (int i=0; i<BASIS.size(); i++)
		outfile << "State n." << i << " : " << BASIS[i] << " En: " << Hint(BASIS[i],BASIS[i]) << endl;
	outfile.close(); */

// Computing BD_BASIS states 

/*	BD_BASIS.resize(0);
	if (DCEV.size()>0)
	{
//	if (N_POL_LEV>0)
//		compute_polaron_BD_BASIS();
//	else
		compute_BD_BASIS();
	}

//  Writing BD_BASIS to logfile

	if (BD_BASIS.size()>0)
	{
		cout << "Writing BD_BASIS to logfile ..." << endl;
		log_filename = prob_name + ".log";
		outfile.open(log_filename.c_str(), ios::app);
		outfile << "---------------------------------------" << endl;
		outfile << " BD BASIS" << endl << endl;
		for (i=0; i<N_MOL; i++)
		outfile << BD_BASIS[i]  << endl;
		outfile.close();
	}*/
	return;
} 

void compute_SYM_BASIS()
{
		cout << "Computing BD_BASIS ..." << endl;
		complex I(0.0,1.0);
		complex energy;
		MPSTATE<MPFLOAT> mps;
		vector_per k_temp;
		Vector kexc, ksym, rhoexc, rhosym;
		int i, j, dc, z;

		BD_BASIS.resize(N_MOL);
		for (i=0; i< N_MOL; i++)
			BD_BASIS[i].resize(0);

		for (i=0; i< BASIS.size(); i++)
		{	
//			cout << i << endl;
			for (dc=0; dc<N_MOL; dc++)
			{
				mps.resize(0);
				for (z=0; z<N_MOL; z++)
					mps.push_back(i+z, udc_kprob(z,dc));
				BD_BASIS[dc].push_back(mps);
				energy=Hint(mps, mps);
				if (abs(energy.imag()/energy.real())>IMAG_TOL)
					nrerror("compute_SYM_BASIS: diagonal energy is not real!");
				BD_BASIS[dc].back().set_energy(energy.real());
				if (BD_BASIS[dc].back().get_energy()>EN_MAX[BASIS[i].get_nph()])
					BD_BASIS[dc].erase(BD_BASIS[dc].end()-1);
			}
			i+=N_MOL-1;
		}	// end for (i=0; i< BASIS.size(); i++)
	return;
}

void compute_symmetric_state(MP& mp_ini, MP& mp_fin, int nsym)
{
	if (nsym < (*latticePtr).get_N_SYM_OP())
	{
		mp_fin.set_k_per(get_sym_k(nsym, mp_ini.get_k()));
		mp_fin.set_mol(get_sym_mol(nsym, mp_ini.get_mol()));
		mp_fin.set_mu(mp_ini.get_mu());
		int num_ph=mp_ini.get_nph();
		mp_fin.set_nph(num_ph);
		for (int np=0; np<num_ph; np++)
			mp_fin.set_ph(np, get_sym_k(nsym,mp_ini.get_ph(np).q), get_sym_mol(nsym, mp_ini.get_ph(np).beta));
		mp_fin.ph_sort();
		return;
	}
	else
		nrerror("compute_symmetric_state: called with an invalid nsym!");
}

complex get_J(vector_per k, int alpha, int beta)
{
	return(Ltilde_all[vector_per_mod::posmod(k.n1(),N1)][vector_per_mod::posmod(k.n2(),N2)][vector_per_mod::posmod(k.n3(),N3)](alpha,beta));
}

ComplexMatrix get_udc(vector_per k)
{
		return(udc_all[vector_per_mod::posmod(k.n1(),N1)][vector_per_mod::posmod(k.n2(),N2)][vector_per_mod::posmod(k.n3(),N3)]);
}

complex Hint(int i, int j)
{
	return(Hint(BASIS[i],BASIS[j]));
}

// This routine only works for 1 vibration for each phonon wave vector and with phonon wave vectors
// ordered in some way !!!!!!!!!!!!!!!!!!!!!!!! 
// Explicit energies have been written only for the case mp1.get_nph()<=mp2.get_nph()
complex Hint(MP &mp1, MP &mp2)
{
	if (mp1.get_nph()>mp2.get_nph())
	{
		return( conj(Hint(mp2, mp1)) );
	}

//	complex I(0.0,1.0);
	int np = mp2.get_nph();
	int npp = mp1.get_nph();
	if (np-npp>1) return(0.0);

//	int check=0;

	int a=mp2.get_mol();
	int ap=mp1.get_mol();
	int nmu = mp2.get_mu();
	int nmup = mp1.get_mu();
	vector_per k0=mp2.get_k();
	vector_per kp=mp1.get_k();
	complex energy=0.0;

	if (npp==0)
	{
		if (np==0)  // this is for 00 interaction
		{
			if (kp == k0) 
			{
				energy+=Sn0_1[nmup]*Sn0_1[nmu]*get_J(kp,ap,a);
				if (ap==a && nmup == nmu )
						energy+=elec_En+nmup*vib_En_1+lambda_0*lambda_0*(1.0-MPFLOAT(EXACT_Q0)/(N1*N2*N3))*vib_En_0;	
			}
//			check+=1;
		}

		if (np==1) // this is for 01 interaction
		{
			if (kp == k0+mp2.get_ph(0).q && ap==a 
				&& ap==mp2.get_ph(0).beta && nmup == nmu
				&& (kp != k0 || EXACT_Q0!=1) ) 
			{
				energy+=-lambda_0*vib_En_0/sqrt(MPFLOAT(N1*N2*N3));//*exp(I*inner_prod(get_k_3d(k0-kp),(*latticePtr).get_RHO(a)));
			}
//			check+=1;
		}

	}

	if (npp==1)
	{
		if (np==1)  // this is for 11 interaction
		{
			if (mp1.get_ph(0).beta==mp2.get_ph(0).beta)
			{
				if (kp == k0 && mp1.get_ph(0).q == mp2.get_ph(0).q  ) 
				{
					energy+=Sn0_1[nmup]*Sn0_1[nmu]*get_J(kp,ap,a);
					if (ap==a && nmup == nmu)
							energy+=elec_En+nmup*vib_En_1+(1+lambda_0*lambda_0*(1.0-MPFLOAT(EXACT_Q0)/(N1*N2*N3)))*vib_En_0;
				}
			}
//			check+=1;
		}
		
		if (np==2)  // this is for 12 interaction
		{
			if (ap==a && nmup == nmu && (kp != k0 || EXACT_Q0!=1) )
			{
				if (
					(kp == k0+mp2.get_ph(1).q && mp1.get_ph(0).q == mp2.get_ph(0).q && 
					ap == mp2.get_ph(1).beta && mp1.get_ph(0).beta ==mp2.get_ph(0).beta)	
					)
					energy+=-lambda_0*vib_En_0/sqrt(MPFLOAT(N1*N2*N3));

				if (
					(kp == k0+mp2.get_ph(0).q && mp1.get_ph(0).q == mp2.get_ph(1).q && 
					ap == mp2.get_ph(0).beta && mp1.get_ph(0).beta ==mp2.get_ph(1).beta)
					)
					energy+=-lambda_0*vib_En_0/sqrt(MPFLOAT(N1*N2*N3));

				if (
					mp2.get_ph(0).q == mp2.get_ph(1).q && mp2.get_ph(0).beta == mp2.get_ph(1).beta
					)
					energy/=sqrt(2.0);
			}
//			check+=1;
		}
	}

	if (npp==2)  
	{
		if (np==2)  // this is for 22 interaction
					// it works only for sorted states !!!!!!!!!!!!!!!!
		{
			if (kp == k0 && mp1.get_ph(0).q == mp2.get_ph(0).q && mp1.get_ph(1).q == mp2.get_ph(1).q &&
				mp1.get_ph(0).beta == mp2.get_ph(0).beta && mp1.get_ph(1).beta == mp2.get_ph(1).beta )
			{
					energy+=Sn0_1[nmup]*Sn0_1[nmu]*get_J(kp,ap,a);
					if (ap == a && nmup == nmu)
						energy+=elec_En+nmup*vib_En_1+(2+lambda_0*lambda_0*(1.0-MPFLOAT(EXACT_Q0)/(N1*N2*N3)))*vib_En_0;
			}
/*			if (kp == k0 && mp1.get_ph(0).q == mp2.get_ph(1).q && mp1.get_ph(1).q == mp2.get_ph(0).q &&
				mp1.get_ph(0).beta == mp2.get_ph(1).beta && mp1.get_ph(1).beta == mp2.get_ph(0).beta )
			{
					energy+=Sn0_1[nmup]*Sn0_1[nmu]*get_J(kp,ap,a);
					if (ap == a && nmup == nmu)
						energy+=elec_En+nmup*vib_En_1+(2+lambda_0*lambda_0*(1.0-MPFLOAT(EXACT_Q0)/(N1*N2*N3)))*vib_En_0;
			}
			if (mp1.get_ph(0).q == mp2.get_ph(0).q && mp1.get_ph(1).q == mp2.get_ph(1).q &&
				mp1.get_ph(0).beta == mp2.get_ph(0).beta && mp1.get_ph(1).beta == mp2.get_ph(1).beta &&
				mp1.get_ph(0).q == mp1.get_ph(1).q && mp1.get_ph(0).beta == mp1.get_ph(1).beta)
			{
					energy-=Sn0_1[nmup]*Sn0_1[nmu]*get_J(kp,ap,a);
					if (ap == a && nmup == nmu)
						energy-=elec_En+nmup*vib_En_1+(2+lambda_0*lambda_0*(1.0-MPFLOAT(EXACT_Q0)/(N1*N2*N3)))*vib_En_0;
			}*/  // This part is only needed if states are unsorted

//			check+=1;
		}
		
		if (np==3)  // this is for 23 interaction
		{
			int b1=mp2.get_ph(0).beta;
			int b2=mp2.get_ph(1).beta;
			int b3=mp2.get_ph(2).beta;

			int b1p=mp1.get_ph(0).beta;
			int b2p=mp1.get_ph(1).beta;

			vector_per q1=mp2.get_ph(0).q;
			vector_per q2=mp2.get_ph(1).q;
			vector_per q3=mp2.get_ph(2).q;

			vector_per q1p=mp1.get_ph(0).q;
			vector_per q2p=mp1.get_ph(1).q;

			if ( ap == a && nmup == nmu && (k0!=kp || EXACT_Q0!=1))
			{

				if ( k0+q1==kp && a==b1 && q2==q1p && b2==b1p && q3==q2p && b3==b2p )
					energy+= 1.0;
				if ( k0+q1==kp && a==b1 && q3==q1p && b3==b1p && q2==q2p && b2==b2p )
					energy+= 1.0;
				if ( k0+q2==kp && a==b2 && q1==q1p && b1==b1p && q3==q2p && b3==b2p )
					energy+= 1.0;
				if ( k0+q2==kp && a==b2 && q3==q1p && b3==b1p && q1==q2p && b1==b2p )
					energy+= 1.0;
				if ( k0+q3==kp && a==b3 && q1==q1p && b1==b1p && q2==q2p && b2==b2p )
					energy+= 1.0;
				if ( k0+q3==kp && a==b3 && q2==q1p && b2==b1p && q1==q2p && b1==b2p )
					energy+= 1.0;

				energy*=-lambda_0*vib_En_0/sqrt(MPFLOAT(N1*N2*N3));

				MPFLOAT f=1.0;
				MPFLOAT fp=1.0;

				if (q1==q2 && b1==b2)	f+= 1.0/sqrt(2.0)-1.0;
				if (q1==q3 && b1==b3)	f+= 1.0/sqrt(2.0)-1.0;
				if (q2==q3 && b2==b3)	f+= 1.0/sqrt(2.0)-1.0;
				if (q1==q2 && b1==b2 && q1==q3 && b1==b3) f+=2.0-3.0/sqrt(2.0)+1.0/sqrt(6.0);

				if (q1p==q2p && b1p==b2p)	fp+= 1.0/sqrt(2.0)-1.0;

				energy*=f*fp;

			} // end if (ap== a && nmup == nmu && k0!=kp)
//			check+=1;
		}
	}

	if (npp==3)
	{
		if (np==3)	// this is for 33 interaction 
					// it works only for sorted states !!!!!!!!!!!!!!!!
		{
			if (kp == k0 && mp1.get_ph(0).q == mp2.get_ph(0).q && mp1.get_ph(1).q == mp2.get_ph(1).q && 
				mp1.get_ph(2).q == mp2.get_ph(2).q && mp1.get_ph(0).beta == mp2.get_ph(0).beta && 
				mp1.get_ph(1).beta == mp2.get_ph(1).beta && mp1.get_ph(2).beta == mp2.get_ph(2).beta )
			{
					energy+=Sn0_1[nmup]*Sn0_1[nmu]*get_J(kp,ap,a);
					if (ap == a && nmup == nmu)
						energy+=elec_En+nmup*vib_En_1+(3+lambda_0*lambda_0*(1.0-MPFLOAT(EXACT_Q0)/(N1*N2*N3)))*vib_En_0;
			}
		}
//		check+=1;
	}

/*	if (check!=1 && abs(mp1.get_nph()-mp2.get_nph())<=1)
	{
		cout << mp1 << endl << mp2 << endl << check << endl;
		nrerror("MultiPhonon::Hint: check!=1!");
	}*/
	return(energy);
}



void  write_spectrum(string spectrum_mode, string input_filename)
{
	if (spectrum_mode == "SL")
		write_spectrum_Single_layer(input_filename, prob_name, dipoles.size(), (*latticePtr).get_gamma_coeff());
	else
	{
		if (spectrum_mode == "N")
			write_Im_n(input_filename, prob_name, dipoles.size(), (*latticePtr).get_epsilon_inf(), (*latticePtr).get_cell_volume(), 
								(*latticePtr).get_gamma_coeff());
//		if (spectrum_mode == "TE")
//			write_spectrum_transverse_epsilon(input_filename, prob_name, dipoles.size(), (*latticePtr).get_epsilon_inf(), (*latticePtr).get_cell_volume(), 
//								(*latticePtr).get_gamma_coeff());
		else
			write_spectrum_epsilon(input_filename, prob_name, dipoles.size(), (*latticePtr).get_epsilon_inf(), (*latticePtr).get_cell_volume(), 
								(*latticePtr).get_gamma_coeff());
	}
	return;
}

complex Hint(MPSTATE<MPFLOAT> &mps1, MPSTATE<MPFLOAT> &mps2)
{
	complex energy = 0.0;
	for (int i = 0; i< mps1.size(); i++)
		for (int j = 0; j< mps2.size(); j++)
			energy += conj(mps1.get_coeff(i))*mps2.get_coeff(j)*Hint(BASIS[mps1.get_state(i)],BASIS[mps2.get_state(j)]);
	return(energy);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Queste due funzioni sono da controllare!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex Hint(int i, int j, int DC)
{
	if (DC<0)
	{
		return(Hint(BASIS[i], BASIS[j]));
	}
	else
	{
		return(Hint(BD_BASIS[DC][i], BD_BASIS[DC][j]));
	}
}
int check_BD_decomposition(int dc1, int dc2, MPFLOAT tol)
{
	int i,j;
	cout << "Checking BD decomposition: " << endl;
	for (i=0; i< BD_BASIS[dc1].size(); i++)
		for (j=0; j< BD_BASIS[dc2].size(); j++)
		{
	//		cout << i << "\t" << j << endl;
			if (abs(Hint(BD_BASIS[dc1][i], BD_BASIS[dc2][j])) > tol)
			{
				cout << i << "\t" << j << "\t" << Hint(BD_BASIS[dc1][i], BD_BASIS[dc2][j]) << endl;
				return(0);
			}
		}
	return(1);
}

}; // End of class MultiPhonon declaration

/*


complex Hint_Hbasis(MPSTATE<MPFLOAT> &mps1, MPSTATE<MPFLOAT> &mps2)
{
	complex energy = 0.0;
	for (int i = 0; i< mps1.size(); i++)
		for (int j = 0; j< mps2.size(); j++)
			energy += conj(mps1.get_coeff(i))*mps2.get_coeff(j)*Hbasis(mps1.get_state(i), mps2.get_state(j));
	return(energy);
}*/

/*
// Checks BD_BASIS
void check_BD_BASIS()
{
	cout << "Checking BD_BASIS set" << endl;
	int i,j, dc1,dc2, max_dc1=-1,max_dc2=-1,max_i=-1,max_j=-1;
	MPFLOAT max_H=-1.0;
	complex temp;
	ofstream HBD_file("HBD_file.dat");

	for (dc1=0; dc1<N_MOL; dc1++)
		for (dc2=dc1; dc2<N_MOL; dc2++)
		{
			max_H=-1.0;
			for (i=0; i<BD_BASIS[dc1].size(); i++)
				for (j=0; j<BD_BASIS[dc2].size(); j++)
				{
					if (HBASIS==1)
						temp=Hint_Hbasis(BD_BASIS[dc1][i], BD_BASIS[dc2][j]);
					else
						temp=Hint_std(BD_BASIS[dc1][i], BD_BASIS[dc2][j]);
//					if (i==0)
//						cout << dc1 << dc2 << " " << " i=" << i << " j=" << j 
//							   << " H=" << temp << endl;
					if ((H_print!="N" || H_print!="n") && abs(temp)>ZERO_TOL ) 
						HBD_file << i << "\t" << j << "\t" << temp << endl;
					if ( abs(temp) > max_H )
					{
						max_H=abs(temp);
						max_dc1=dc1;
						max_dc2=dc2;
						max_i=i;
						max_j=j;
					}
				}
			cout << "dc1: " << max_dc1 << " dc2: " << max_dc2 << " Largest off diagonal term: i=" 
				 << max_i << " j=" << max_j << " H=" << max_H << endl;
		}
	HBD_file.close();
//	nrerror("check_BD_BASIS(): found H matrix element between different dc!");					
	return;
}
*/

/*
complex Hij_bd(int i, int j)
{
//	complex val;
	if (HBASIS==1)
		return(Hint_Hbasis(BD_BASIS[DC][i], BD_BASIS[DC][j]));
	else
		return(Hint_std(BD_BASIS[DC][i], BD_BASIS[DC][j]));
//	cout << "i j val: " << i << " " << j << " " << val << endl;
//	return(val);
}*/
/*
complex Hbasis(int i, int j)
{
	int n;
	if (j>=i)
	{
		n=HB_rowInd[i];
		while (n<HB_rowInd[i+1])
		{		
			if (HB_col[n]>j)
				return(0.0);
			if (HB_col[n]==j)
				return(HB_val[n]);
			n++;
		}
		return(0.0);
	}
	else 
		return(conj(Hbasis(j,i)));
}

void compute_Hbasis()
{
	cout << "Computing Hbasis ... " << endl;
	HB_col.resize(0);
	HB_rowInd.resize(0);
	HB_val.resize(0);

	int i,j;
	complex val;
	for (i=0; i<BASIS.size(); i++)
	{
		HB_rowInd.push_back(HB_val.size());
		for (j=i; j<BASIS.size(); j++)
		{
			val=Hint(BASIS[i],BASIS[j]);
			if (abs(val)>ZERO_TOL)
			{
				HB_val.push_back(val);
				HB_col.push_back(j);
			}
		}
	}
	HB_rowInd.push_back(HB_val.size());
//	cout << "Checking Hbasis" << endl;
//	for (i=0; i<BASIS.size(); i++)
//		for (j=0; j<BASIS.size(); j++)
//			if (abs(Hbasis(i,j)-Hint(BASIS[i],BASIS[j]))>ZERO_TOL)
//				cout << i << " " << j << " error" << endl;
	if (H_print!="N" || H_print!="n") 
		print_Hbasis();
	return;
}

void set_H_print(string _H_print) { H_print=_H_print; };

void print_Hbasis()
{
	cout << "Printing Hbasis ..." << endl;
	ofstream HB_file("Hbasis.dat"); 
	int i,j;
	complex val;
	for (i=0; i<BASIS.size(); i++)
		for (j=i; j<BASIS.size(); j++)
		{
			val=Hbasis(i,j);
			if (abs(val)>ZERO_TOL)
				HB_file << i << "\t" << j << "\t" << val << endl;
		}
	HB_file.close();
	return;
}*/

/*
void solve()
{
	cout << "Solving MultiPhonon problem ...." << endl;
	if (HBASIS==1)
	{
//		if (NEW_METHOD==1 && N_POL_LEV>0)
//			compute_H1basis();
//		else
			compute_Hbasis();
	}

//	if (BD_BASIS.size()>0 && SUPER_BASIS==1)
//		compute_super_BASIS();   // This line computes the approximated states as described in my report

//	if (DCEV.size()==0 || BD_BASIS.size()==0 )
		solve_BASIS();
//	else
//		solve_BD_BASIS();
	return;
}*/

/*

void compute_super_BASIS()
{
	int i,j,dc;
	vector<MPSTATE> super1(N_MOL), super2(N_MOL), super3(N_MOL);
	vector<int> count(N_MOL);
//	ofstream testfile("testsuper.txt");
	if (BD_BASIS.size()>0)
	{
		cout << "Computing super BASIS ... " << endl;
		if (MAX_PH>0)
		{
			for (dc=0; dc<N_MOL; dc++)
			{
				count[dc]=1;
				while (BASIS[BD_BASIS[dc][count[dc]].get_state(0)].get_nph()==1 && count[dc]<BD_BASIS[dc].size())
				{
					super1[dc]=super1[dc]+BD_BASIS[dc][count[dc]]*Hint(BD_BASIS[dc][0],BD_BASIS[dc][count[dc]]);
					count[dc]++;
					if ( count[dc]>= BD_BASIS[dc].size() )
						break;
				}
				super1[dc].normalize();
				super1[dc].set_energy(Hint(super1[dc],super1[dc]));
			}
		}

		if (MAX_PH>1)
		{
			for (dc=0; dc<N_MOL; dc++)
			{
				while (BASIS[BD_BASIS[dc][count[dc]].get_state(0)].get_nph()==2 )
				{
					super2[dc]=super2[dc]+BD_BASIS[dc][count[dc]]*Hint(super1[dc],BD_BASIS[dc][count[dc]]);
					count[dc]++;
					if ( count[dc]>= BD_BASIS[dc].size() )
						break;
				}
				super2[dc].normalize();
				super2[dc].set_energy(Hint(super2[dc],super2[dc]));
			}
		}

		if (MAX_PH>2)
		{
			for (dc=0; dc<N_MOL; dc++)
			{
				while (BASIS[BD_BASIS[dc][count[dc]].get_state(0)].get_nph()==3 )
				{
					super3[dc]=super3[dc]+BD_BASIS[dc][count[dc]]*Hint(super2[dc],BD_BASIS[dc][count[dc]]);
//					testfile << dc << "\t" << count[dc] << "\t" << Hint(super2[dc],BD_BASIS[dc][count[dc]]) << endl;
					count[dc]++;
					if ( count[dc]>= BD_BASIS[dc].size() )
						break;
				}
				super3[dc].normalize();
				super3[dc].set_energy(Hint(super3[dc],super3[dc]));
			}
		}

//		testfile.close();

		for (dc=0; dc<N_MOL; dc++)
		{
			BD_BASIS[dc].erase(BD_BASIS[dc].begin()+1, BD_BASIS[dc].end());
//			BD_BASIS[dc].erase(BD_BASIS[dc].begin()+1, BD_BASIS[dc].begin()+count[dc]);
			if (MAX_PH>0)	
				BD_BASIS[dc].push_back(super1[dc]);
			//	BD_BASIS[dc].insert(BD_BASIS[dc].begin()+1,super1[dc]);
				
			if (MAX_PH>1)
				BD_BASIS[dc].push_back(super2[dc]);
			//	BD_BASIS[dc].insert(BD_BASIS[dc].begin()+2,super2[dc]);

			if (MAX_PH>2)
				BD_BASIS[dc].push_back(super3[dc]);
		}

		ofstream outfile("Super_BASIS.txt");
		outfile << "---------------------------------------" << endl;
		outfile << " Super BASIS" << endl << endl;
		for (i=0; i<N_MOL; i++)
			outfile << BD_BASIS[i]  << endl;
		outfile.close();
	}
	return;
}*/
////////////////////////////////////////////////////////////
// H1 functions

/*
// This routine only works for 1 vibration for each phonon wave vector and with phonon wave vectors
// ordered in some way !!!!!!!!!!!!!!!!!!!!!!!! 
// Explicit energies have been written only for the case nph < mp.get_nph()

complex H1int(MP &mp1, MP &mp2)
{

	int np = mp2.get_nph();
	int npp = mp1.get_nph();
	if (abs(np-npp)!=1)
		return(0.0);
	if (np<npp)
		return(H1int(mp2,mp1));
//	if (np-npp!=1) 
//		nrerror("H1int: unexpected number of phonons!");

//	int check=0;

	int a=mp2.get_mol();
	int ap=mp1.get_mol();
	int nmu = mp2.get_mu();
	int nmup = mp1.get_mu();
	vector_per k0=mp2.get_k();
	vector_per kp=mp1.get_k();
	MPFLOAT energy=0.0;

	if (npp==0) // this is for 01 interaction
	{

			if (kp == k0+mp2.get_ph(0).q && ap==a 
				&& ap==mp2.get_ph(0).beta && nmup == nmu
				&& (kp != k0 || EXACT_Q0!=1) ) 
			{
				energy+=-lambda_0*vib_En_0/sqrt(MPFLOAT(N1*N2*N3));
			}
//			check+=1;
	}

	if (npp==1)  // this is for 12 interaction
	{

			if (ap==a && nmup == nmu && (kp != k0 || EXACT_Q0!=1) )
			{
				if (
					(kp == k0+mp2.get_ph(1).q && mp1.get_ph(0).q == mp2.get_ph(0).q && 
					ap == mp2.get_ph(1).beta && mp1.get_ph(0).beta ==mp2.get_ph(0).beta)	
					)
					energy+=-lambda_0*vib_En_0/sqrt(MPFLOAT(N1*N2*N3));

				if (
					(kp == k0+mp2.get_ph(0).q && mp1.get_ph(0).q == mp2.get_ph(1).q && 
					ap == mp2.get_ph(0).beta && mp1.get_ph(0).beta ==mp2.get_ph(1).beta)
					)
					energy+=-lambda_0*vib_En_0/sqrt(MPFLOAT(N1*N2*N3));

				if (
					mp2.get_ph(0).q == mp2.get_ph(1).q && mp2.get_ph(0).beta == mp2.get_ph(1).beta
					)
					energy/=sqrt(2.0);
			}
//			check+=1;
	}

	if (npp==2)  // this is for 23 interaction
	{

			int b1=mp2.get_ph(0).beta;
			int b2=mp2.get_ph(1).beta;
			int b3=mp2.get_ph(2).beta;

			int b1p=mp1.get_ph(0).beta;
			int b2p=mp1.get_ph(1).beta;

			vector_per q1=mp2.get_ph(0).q;
			vector_per q2=mp2.get_ph(1).q;
			vector_per q3=mp2.get_ph(2).q;

			vector_per q1p=mp1.get_ph(0).q;
			vector_per q2p=mp1.get_ph(1).q;

			if ( ap == a && nmup == nmu && (k0!=kp || EXACT_Q0!=1))
			{

				if ( k0+q1==kp && a==b1 && q2==q1p && b2==b1p && q3==q2p && b3==b2p )
					energy+= 1.0;
				if ( k0+q1==kp && a==b1 && q3==q1p && b3==b1p && q2==q2p && b2==b2p )
					energy+= 1.0;
				if ( k0+q2==kp && a==b2 && q1==q1p && b1==b1p && q3==q2p && b3==b2p )
					energy+= 1.0;
				if ( k0+q2==kp && a==b2 && q3==q1p && b3==b1p && q1==q2p && b1==b2p )
					energy+= 1.0;
				if ( k0+q3==kp && a==b3 && q1==q1p && b1==b1p && q2==q2p && b2==b2p )
					energy+= 1.0;
				if ( k0+q3==kp && a==b3 && q2==q1p && b2==b1p && q1==q2p && b1==b2p )
					energy+= 1.0;

				energy*=-lambda_0*vib_En_0/sqrt(N1*N2*N3);

				MPFLOAT f=1.0;
				MPFLOAT fp=1.0;

				if (q1==q2 && b1==b2)	f+= 1.0/sqrt(2.0)-1.0;
				if (q1==q3 && b1==b3)	f+= 1.0/sqrt(2.0)-1.0;
				if (q2==q3 && b2==b3)	f+= 1.0/sqrt(2.0)-1.0;
				if (q1==q2 && b1==b2 && q1==q3 && b1==b3) f+=2.0-3.0/sqrt(2.0)+1.0/sqrt(6.0);

				if (q1p==q2p && b1p==b2p)	fp+= 1.0/sqrt(2.0)-1.0;

				energy*=f*fp;

			} // end if (ap== a && nmup == nmu && k0!=kp)
//			check+=1;
	}

	return(energy);
}*/
/*

void compute_H1basis()
{
	cout << "Computing H1basis ... " << endl;
	HB_col.resize(0);
	HB_rowInd.resize(0);
	HB_val.resize(0);

	int i, j, mu, nph;
	int mu_block=BASIS.size()/(MAX_HI_VIB+1);
	vector<int> istart(MAX_PH+2);
	istart[0]=0;
	for (i=1; i<=MAX_PH; i++)
		istart[i]=istart[i-1]+N_STATES[i-1]/(MAX_HI_VIB+1);
	istart[MAX_PH+1]=mu_block;

//	cout << "istart :" << istart << endl;
//	cout << "mu_block :" << mu_block << endl;
	MPFLOAT val;
	for (mu=0; mu<=max(0,MAX_HI_VIB); mu++) 
	{
		for (nph=0; nph<MAX_PH; nph++)
		{
			for (i=mu*mu_block+istart[nph]; i<mu*mu_block+istart[nph+1]; i++)
			{
//				cout << "i: " << i << endl;
				HB_rowInd.push_back(HB_val.size());
				for (j=mu*mu_block+istart[nph+1]; j<mu*mu_block+istart[nph+2]; j++)
				{
	//				cout << "j: " << j << endl;
					val=H1int(BASIS[i],BASIS[j]);
					if (abs(val)>ZERO_TOL)
					{
						HB_val.push_back(val);
						HB_col.push_back(j);
					}
				}
			}
		} // end for (nph=0; nph<MAX_PH; nph++)
		
		for (i=mu*mu_block+istart[MAX_PH]; i<(mu+1)*mu_block; i++)
		{
//			cout << "i: " << i << endl;
			HB_rowInd.push_back(HB_val.size());
		}

	} // end for (mu=0; mu<=max(0,MAX_HI_VIB); mu++) 
	HB_rowInd.push_back(HB_val.size());

//	cout << "Checking H1basis" << endl;
//	for (i=0; i<BASIS.size(); i++)
//		for (j=0; j<BASIS.size(); j++)
//			if (abs(Hbasis(i,j)-Hint(BASIS[i],BASIS[j]))>ZERO_TOL)
//				cout << i << " " << j << " error" << endl; 
	return;
} */
/*

complex H1int_slow(MPSTATE<MPFLOAT> &mps1, MPSTATE<MPFLOAT> &mps2)
{
	complex energy = 0.0;
	for (int i = 0; i< mps1.size(); i++)
		for (int j = 0; j< mps2.size(); j++)
			energy += mps1.get_coeff(i)*mps2.get_coeff(j)*H1int(BASIS[mps1.get_state(i)],BASIS[mps2.get_state(j)]);
	return(energy);
}
*/
/*

complex H1ij_bd(int i, int j)
{
	double val;
	if (HBASIS==1)
	{
		if (i==j)
			return(BD_BASIS[DC][i].get_energy());
		else
		{
			if (abs(BASIS[BD_BASIS[DC][i].get_state(0)].get_nph()
				-BASIS[BD_BASIS[DC][j].get_state(0)].get_nph())!=1)
				return(0.0);
			else
				return(Hint(BD_BASIS[DC][i], BD_BASIS[DC][j]));
		}
// This is the old working and checked routine  
//		if (abs(BASIS[BD_BASIS[DC][i].get_state(0)].get_nph()
//			   -BASIS[BD_BASIS[DC][j].get_state(0)].get_nph())>1)
//			return(0.0);
//		else
//		{
//			if (i==j)
//				return(BD_BASIS[DC][i].get_energy());
//			else
//				return(Hint(BD_BASIS[DC][i], BD_BASIS[DC][j]));
//		}
	}
	else
		return(H1int_slow(BD_BASIS[DC][i], BD_BASIS[DC][j]));
//	cout << "i j val: " << i << " " << j << " " << val << endl;
//	return(val);
}*/
/////////////////////////////////////////////////////////////////////////////////
// Polaron functions

/*

double get_pol_eigenval(vector_per k, int dc, int lev)
{
	return(pol_eigenval[vector_per_mod::posmod(k.n1(),N1)*(N2*N3*N_MOL)+vector_per_mod::posmod(k.n2(),N2)*(N3*N_MOL)
				        +vector_per_mod::posmod(k.n3(),N3)*N_MOL+dc][lev]);
}


complex<double> get_pol_eigenvec(vector_per k, int dc, int lev, int n_coeff)
{
	return(pol_eigenvec[vector_per_mod::posmod(k.n1(),N1)*(N2*N3*N_MOL)+vector_per_mod::posmod(k.n2(),N2)*(N3*N_MOL)
			+vector_per_mod::posmod(k.n3(),N3)*N_MOL + dc][lev*(MAX_HI_VIB+1)+n_coeff]);
}*/
/*

void compute_polaronic_states(vector_per k)
{
	cout << "Computing polaronic states" << endl;

	udc_k = get_udc(k).copy();

	int nm=N_MOL;
	int nstates=nm*(MAX_HI_VIB+1);		// n. di stati con 0 fononi
	int count=0;
	int nmol, mu;

	BASIS.resize(nstates);

// 0 phonon states

	for (mu=0; mu<=MAX_HI_VIB; mu++)
		for (nmol=0; nmol<nm; nmol++)
		{
			BASIS[count].set_k_per(k);
			BASIS[count].set_mol(nmol);
			BASIS[count].set_mu(mu);
			BASIS[count].set_nph(0);
			count++;
		}

// Computing BD_BASIS states with respect to symmetry

	int i, j, dc, z;
	if ((*latticePtr).get_N_SYM_OP()>1)
	{
		cout << "Computing 0 phonon BD_BASIS ..." << endl;
		MPSTATE mps;

		BD_BASIS.resize(N_MOL);
		for (i=0; i< N_MOL; i++)
			BD_BASIS[i].resize(0);

		for (i=0; i< BASIS.size(); i++)
		{
			for (dc=0; dc<N_MOL; dc++)
			{
					mps.resize(0);
					for (z=0; z<nm; z++)
						mps.push_back(i+z, udc_k[z][dc]);
					BD_BASIS[dc].push_back(mps);
			}
			i+=nm-1;
		}	// end for (i=0; i< BASIS.size(); i++)
	}	// end if ((*latticePtr).get_N_SYM_OP()>1)

// Solve problem

	if (HBASIS==1)
		compute_Hbasis();
	if (CHECK_BD==1 && BD_BASIS.size()>0)
		check_BD_BASIS();
	for (dc=0; dc<nm; dc++)
			solve_bd(dc, BASIS.size(), 'P');

	return;
}*/
/*

void compute_polaron_BD_BASIS()
{

// Defining a map for BASIS
	
	int i;
	map <string, int> bmap;
	for (i=0; i<BASIS.size(); i++)
		bmap[BASIS[i].print_string()]=i;

// Computing BD_BASIS states for polaronic states (it works? only for 4T planar)
// in this routine mu corresponds to the polaronic level and alpha to the polaronic dc

	vector<int> flag(BASIS.size());			// Flags tells if states have already been included in BD_BASIS
	for (i=0; i<BASIS.size(); i++)	flag[i]=0;

	if ((*latticePtr).get_N_SYM_OP()>1)
	{
		cout << "Computing polaron BD_BASIS ..." << endl;
		MPSTATE mps1, mps2, mps;
		vector_per k_temp;
		vector<int> mystates(2);	// 2 polaronic states to be canceled from list by number
		vector<MP> mympp(2);		// 2 polaronic states to be canceled from list by description
		int j, dc, mu;
		int pdc, lev, lev_ok;
		int shift=BASIS.size()/(MAX_HI_VIB+1);
		MP mymp, mp;
		int nmp, alpha;

		BD_BASIS.resize(N_MOL);
		for (i=0; i< N_MOL; i++)
			BD_BASIS[i].resize(0);

		for (i=0; i< BASIS.size(); i++)
		{
			lev=BASIS[i].get_mu();
			pdc=BASIS[i].get_mol();
			lev_ok=0;
			for (j=0; j<N_POL_LEV; j++)
				if (POL_LEVS[j]==lev)
					lev_ok=1;
		if (flag[i]==0 && lev_ok==1)
		{
	//		cout << "Stato polaronico: " << i << " LEV: " << lev << endl;
// Here we are talking about polarons
			mystates[0]=i;
			mympp[0]=BASIS[i];
			compute_symmetric_state(mympp[0], mympp[1], 1);
			mympp[1].set_mol(mympp[0].get_mol());
			for (j=1; j<2; j++)
				mystates[j]=bmap[mympp[j].print_string()];
			for (j=0; j<2; j++)
				flag[mystates[j]]=1;
	//		cout << "Stati eliminati: " << mystates << endl;
// Here we are back to multiphonon
			for (dc=0; dc<N_MOL; dc++)
			{
					mps1.resize(0);
					mps2.resize(0);
					mp=mympp[0];
					mp.set_mu(0);
					for (alpha=0; alpha<N_MOL; alpha++)
					{
						mp.set_mol(alpha);
						nmp=bmap[mp.print_string()];
			//			cout << "Stato BASIS: " << nmp << endl; 
						for (mu=0; mu<=MAX_HI_VIB; mu++)
							mps1.push_back(nmp+mu*shift, udc_k[0][dc]*
														 get_udc(mp.get_k())[alpha][pdc]*
														 get_pol_eigenvec(mp.get_k(),pdc, lev, mu) );
					}
					mp=mympp[1];
					mp.set_mu(0);
					for (alpha=0; alpha<N_MOL; alpha++)
					{
						mp.set_mol(alpha);
						nmp=bmap[mp.print_string()];
		//				cout << "Stato BASIS: " << nmp << endl; 
						for (mu=0; mu<=MAX_HI_VIB; mu++)
							mps2.push_back(nmp+mu*shift, udc_k[1][dc]*
							get_udc(get_sym_k(1,mp.get_k()))[get_sym_mol(1,alpha)][pdc]*
						                  get_pol_eigenvec(mp.get_k(),pdc, lev, mu) );
					}
					mps=mps1+mps2;
					if (mps.normalize())
					{
	//					mps.set_energy(get_pol_eigenval(mp.get_k(), pdc, lev)
	//						           +vib_En_0*mp.get_nph() );
						mps.set_energy(Hint_slow(mps,mps));
						BD_BASIS[dc].push_back(mps);
					}
				}

		}	// end if (flag[i]==0)
		}	// end for (i=0; i< BASIS.size(); i++)

	}	// end if ((*latticePtr).get_N_SYM_OP()>1)

	return;
}*/

#endif