/* Written by Leonardo Silvestri 2007-2013 */

#include <boost/math/special_functions/erf.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <Complex>
#include <math.h>
#include <vector>
#include <OPutils.h>
#include <N_Functions.h>
#include <physics.h>
#include <lattice.h>

namespace ublas = boost::numeric::ublas;
namespace math = boost::math;
using namespace std;

#ifndef Pi
#define Pi 3.141592654
#endif 
#ifndef ZERO_TOL
#define ZERO_TOL 1e-10
#endif 


//////////////////////////////////////////////////
/////////     Lattice class functions	 /////////
//////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////
//	Functions that set internal variables.

template<typename LATFLOAT>
void lattice<LATFLOAT>::initialize(string input_file)
	{ 
		lattice_name=input_file.substr(0,input_file.find('.'));
		cout << "Reading lattice data from file: " << input_file << endl;

		int i,j,k,n;
		string line;
		ifstream infile(input_file.c_str());
		if (!infile.is_open()) nrerror("File not found!");

		Vector				temp(3);
		sym_op<LATFLOAT>	sym_temp;
		Matrix				mat_temp(3,3);	
		LATFLOAT			d_temp;
		int					i_temp;
		pos<LATFLOAT>		pos_temp;

		getline(infile, line, ':');
		CELL.resize(3);
		for (i=0;i<3;i++) 
		{
			CELL[i].resize(3);
			for (j=0;j<3;j++) infile >> CELL[i](j);
		}
		getline(infile, line, ':');
		infile >> N_MOL;
		getline(infile, line, ':');
		infile >> N_SP;
		getline(infile, line, ':');
		infile >> N_SYM_OP;

		if (N_SP*N_SYM_OP!=N_MOL)
			nrerror("lattice<LATFLOAT>::read_crystal_data: Symmetry operations generate more molecules than needed");	
		SYM.resize(0);
		for (k=0;k<N_SYM_OP;k++) 
		{
			getline(infile, line, ':');
			for (i=0;i<3;i++) 
				for (j=0;j<3;j++) 
					infile >> mat_temp(i,j);
			sym_temp.set_rot(mat_temp);
			for (j=0;j<3;j++) infile >> temp[j];
			sym_temp.set_ax(CELL[0]*temp[0]+CELL[1]*temp[1]+CELL[2]*temp[2]);

			for (j=0;j<3;j++) infile >> temp[j];
			sym_temp.set_tr(CELL[0]*temp[0]+CELL[1]*temp[1]+CELL[2]*temp[2]);
			sym_temp.set_map_size(N_MOL);
			for (j=0;j<N_MOL;j++) 
			{
				infile >> i_temp;
				sym_temp.set_map(j,i_temp);
			}
			SYM.push_back(sym_temp);
		}

		OtoF.resize(3,3);
		for (i=0;i<3;i++) 
			for (j=0;j<3;j++) 
				OtoF(i,j)=CELL[j][i];
		InvertMatrix(OtoF,OtoF);

		getline(infile, line, ':');
		RHO.resize(0);
		for (i=0;i<N_SP;i++) 
		{
			for (j=0;j<3;j++) infile >> temp[j];	// Input data are in fractional coordinates
			FracToOrth(temp,temp);					// temp is converted into orthogonal coordinates
			move_into_primitive_cell_orth(temp);
			RHO.push_back(temp);
			for (k=1; k< N_MOL/N_SP; k++)
			{
				pos_temp.alpha=i*N_MOL/N_SP;
				pos_temp.n=temp;
				SYM[k].apply(pos_temp, pos_temp);
				move_into_primitive_cell_orth(pos_temp.n);
				RHO.push_back(pos_temp.n); 
			}
		}	

		getline(infile, line, ':');
		MOL_DIPOLE.resize(0);
		for (i=0;i<N_SP;i++) 
		{
			for (j=0;j<3;j++) infile >> temp[j];
			MOL_DIPOLE.push_back(temp);
			for (k=1; k< N_MOL/N_SP; k++)
				MOL_DIPOLE.push_back(prod(SYM[k].get_rot(),temp));  
		}	

		getline(infile, line, ':');
		for (i=0;i<N_SP;i++) 
		{
			infile >> d_temp;
			for (k=0; k<N_MOL/N_SP; k++) 
				for (j=0;j<3;j++)
					MOL_DIPOLE[k+i*N_MOL/N_SP][j]*=d_temp;
		};

		ATOM.resize(N_MOL);
		CHARGE.resize(N_MOL);
		for (i=0;i<N_SP;i++) 
		{
			getline(infile, line, ':');
			infile >> i_temp;
			for (k=0; k< N_MOL/N_SP; k++)
			{
					ATOM[k+i*N_MOL/N_SP].resize(i_temp);
					CHARGE[k+i*N_MOL/N_SP].resize(i_temp);
			}
			for (n=0; n<i_temp; n++)
			{
				for (j=0;j<3;j++) infile >> temp[j];
				infile >> d_temp;
				CHARGE[i*N_MOL/N_SP][n]=d_temp;
				ATOM[i*N_MOL/N_SP][n]=temp;
				for (k=1; k< N_MOL/N_SP; k++)
				{
					ATOM[k+i*N_MOL/N_SP][n]=prod(SYM[k].get_rot(),temp);
					CHARGE[k+i*N_MOL/N_SP][n]=d_temp;
				}
			}
		}	

		getline(infile, line, ':');
		epsilon_0.resize(3,3);
		for (i=0;i<3;i++) 
			for (j=0;j<3;j++) 
				infile >> epsilon_0(i,j);

		getline(infile, line, ':');	
		infile.close();
		cout << "Finished reading!" << endl;
		
		cell_volume=inner_prod(CELL[0],cross_prod(CELL[1],CELL[2]));
		KCELL.resize(3);
		KCELL[0]=cross_prod(CELL[1],CELL[2])*(2.0*Pi/cell_volume);
		KCELL[1]=cross_prod(CELL[2],CELL[0])*(2.0*Pi/cell_volume);
		KCELL[2]=cross_prod(CELL[0],CELL[1])*(2.0*Pi/cell_volume);

		int flag_2D=0;
		LATFLOAT rmax=max(max(norm_2(CELL[0]),norm_2(CELL[1])),norm_2(CELL[2]));
		compute_NN(rmax,flag_2D);
		string filename = lattice_name + ".log";
		ofstream log_file(filename.c_str());
		print(log_file);
		log_file.close();
		return;
	};
	
// Computes the nearest neighbours for each inequivalent molecule
// Check if it exploits all the symmetry operations 
// If flag_2D==1 it considers only NN on the same xy plane
template<typename LATFLOAT>	
void lattice<LATFLOAT>::compute_NN(LATFLOAT rmax, int flag_2D)
	{
		cout << "Computing N*N with rmax=" << rmax << " and flag_2D=" << flag_2D << endl;
		int alpha;
		int i,j,k,n,m,l;
		int max_a, max_b,max_c;
		int sp, mol;
		LATFLOAT PLANAR_TOL=norm_2(CELL[2])/2.0;
		pos<LATFLOAT> position;
		typename vector<pos<LATFLOAT> >::iterator it;

		if (rmax<0)
			rmax=ZERO_TOL;

		NN.resize(N_MOL);
		for (i=0; i<N_MOL; i++)
			NN[i].resize(0);

		max_a=int(rmax/norm_2(CELL[0]))+1;
		max_b=int(rmax/norm_2(CELL[1]))+1;
		if (flag_2D!=1)
			max_c=int(rmax/norm_2(CELL[2]))+1;
		else
			max_c=1;

		for (sp=0;sp<N_SP;sp++) 
		{
			mol=sp*N_MOL/N_SP;

	// Loop over all the molecules in the crystal
			for (alpha=0; alpha<N_MOL; alpha++) 
				for (n=-max_a; n<=max_a; n+=1) 
					for (m=-max_b; m<=max_b; m+=1)  
						for (l=-max_c; l<=max_c; l+=1) 
							{
								position.alpha=alpha;
								position.n=n*CELL[0]+m*CELL[1]+l*CELL[2]+RHO[alpha]-RHO[mol];
								if (flag_2D!=1 || abs(position.n[2])<PLANAR_TOL)
								if ( norm_2(position.n)<rmax )
								{
									for (it=NN[mol].begin(); it!=NN[mol].end(); it++ ) 
										if ( norm_2((*it).n)>norm_2(position.n) ) break;
									NN[mol].insert(it,position);
								}
							}

// Fixes other molecules basis states so that symmetry operations leave the phonon cloud the same
			for (k=1; k< N_MOL/N_SP; k++)
			{	
				n=SYM[k].get_map()[0];
				NN[mol+n].resize(NN[mol].size());
				for (j=0; j<NN[mol+n].size(); j++)	
				{
					NN[mol+n][j].alpha=SYM[k].get_map(NN[mol][j].alpha);
					NN[mol+n][j].n=prod(SYM[k].get_rot(),NN[mol][j].n);
				}
			}
		} // end loop over species	
		NNrmax=rmax;		// sets internal variable NNrmax
		NNflag_2D=flag_2D;	// sets internal variable NNflag_2D
		return;
	}

template<typename LATFLOAT>	
void lattice<LATFLOAT>::print(ofstream& log_file)
	{
		cout << "Printing lattice data ... " << endl;
//		double energy; 
		unsigned int i,j; //beta, nn, mol;
//		string filename = lattice_name + ".log";
//		ofstream log_file(filename.c_str());
		log_file << "==================================================================================================" << endl;
		log_file << "     Lattice  " << lattice_name << endl;
		log_file << "==================================================================================================" << endl;
		log_file << endl;		
		log_file << "Static epsilon: " << epsilon_0 << endl;
		log_file << "--------------------------- Crystal cell " << endl; 
		for (i=0; i<CELL.size(); i++)
			log_file << CELL[i] << endl;
		log_file << "Cell volume: " << cell_volume << " A^3" << endl;
		log_file << endl << "--------------------------- Crystal reciprocal lattice cell " << endl; 
		for (i=0; i<KCELL.size(); i++)
			log_file << KCELL[i] << endl;
		log_file << endl << "--------------------------- OtoF (transformation matrix from Orthogonal to Fractional coordinates) " << endl; 
		log_file << OtoF << endl ;
		log_file << endl << "--------------------------- Molecule position inside the cell (orthogonal coordinates)" << endl;
		for (i=0; i<RHO.size(); i++)
			log_file << RHO[i] << endl;
		log_file << endl << "--------------------------- Molecular dipole moments (in Debyes)" << endl;
		for (i=0; i<MOL_DIPOLE.size(); i++)
			log_file << MOL_DIPOLE[i] << endl;
		log_file << endl << "------------- Atomic positions with respect to the CoM (orthogonal coordinates in A) and charges (in e)" << endl;
		for (i=0; i<N_MOL; i++)
		{
			log_file << "MOLECULE: " << i << endl; 
			for (j=0; j<ATOM[i].size(); j++)
				log_file << ATOM[i][j] << "\t" << CHARGE[i][j] << endl;
		}
		log_file << endl << "--------------------------- Symmetry operations " << endl;
		for (i=0; i<SYM.size(); i++)
			log_file << SYM[i] << endl;
		log_file << endl << "--------------------------- Nearest neighbours " << endl;
		for (i=0; i<NN.size(); i++)
		{
			log_file << "MOLECULE: " << i << endl; 
			for (j=1; j<NN[i].size(); j++)
			{
				log_file << j << ": " << NN[i][j] 
				<< "\tT: " << distributed_charge_interaction(i, NN[i][j].alpha, -NN[i][j].n) 
				<< "\tD: " << dipole_dipole(MOL_DIPOLE[i], MOL_DIPOLE[NN[i][j].alpha], NN[i][j].n)
				<< "\tST: " << screened_distributed_charge_interaction(i, NN[i][j].alpha, -NN[i][j].n, epsilon_0(0,0), epsilon_0(1,1), epsilon_0(2,2)) 
				<< "\tSD: " << screened_dipole_dipole(MOL_DIPOLE[i], MOL_DIPOLE[NN[i][j].alpha], NN[i][j].n, epsilon_0)
				<< endl;
			}
		}
//		log_file.close();
		return;
	}

////////////////////////////////////////////////////////////////////////////////////
// Functions to set parameters

template<typename LATFLOAT>	
void lattice<LATFLOAT>::set_NN(vector<vector<pos<LATFLOAT> > > new_NN)
	{
		if (NN.size()!=new_NN.size())
			nrerror("lattice::set_NN: wrong NN size!");
		int i,j;
		for (i=0; i<NN.size(); i++)
		{
			NN[i].resize(new_NN[i].size());
			for (j=0; j<new_NN[i].size(); j++)
				NN[i][j]=new_NN[i][j];
		}
		return;
	}

////////////////////////////////////////////////////////////////////////////////////
// Interactions

template<typename LATFLOAT>	
LATFLOAT lattice<LATFLOAT>::interaction( int alpha, int beta, Vector dr, 
					  string mol_mol_int, LATFLOAT NON_SCR_RADIUS, vector<vector<LATFLOAT> >& NN_int)
{
		switch(mol_mol_int[0])
		{
			case 'N':  // Nearest neighbors interaction
				for (int i=0; i<min(NN[alpha].size(),NN_int[alpha].size()); i++)
					if (norm_2(NN[alpha][i].n-dr)<ZERO_TOL)
					{
						if (NN[alpha][i].alpha==beta)
							return(NN_int[alpha][i]);
						else
							return(0.0);
					}
				return(0.0);
			case 'T':  // Transition charge distribution
				if (norm_2(dr)>NON_SCR_RADIUS)
					return(screened_distributed_charge_interaction(alpha, beta, dr, 
					       epsilon_0(0,0), epsilon_0(1,1), epsilon_0(2,2)));
				else
					return(distributed_charge_interaction(alpha, beta, dr));

			default: // Dipole-dipole interaction
				if (norm_2(dr)>NON_SCR_RADIUS)
					return(screened_dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta],dr, epsilon_0));
				else
					return(dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta],dr));
		}
}

template<typename LATFLOAT>	
LATFLOAT lattice<LATFLOAT>::interaction( int alpha, int nn, 
					  string mol_mol_int, LATFLOAT NON_SCR_RADIUS, vector<vector<LATFLOAT> >& NN_int)
{
	if (nn<NN[alpha].size())
	{
		Vector dr=NN[alpha][nn].n;
		int beta=NN[alpha][nn].alpha;
		switch(mol_mol_int[0])
		{
			case 'N':  // Nearest neighbors interaction
				if (nn<NN_int[alpha].size())
						return(NN_int[alpha][nn]);
					else
						return(0.0);
			case 'T':  // Transition charge distribution
				if (norm_2(dr)>NON_SCR_RADIUS)
					return(screened_distributed_charge_interaction(alpha, beta, dr, 
					       epsilon_0(0,0), epsilon_0(1,1), epsilon_0(2,2)));
				else
					return(distributed_charge_interaction(alpha, beta, dr));

			default: // Dipole-dipole interaction
				if (norm_2(dr)>NON_SCR_RADIUS)
					return(screened_dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta],dr, epsilon_0));
				else
					return(dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta],dr));
		}
	}
	else
		return(0.0);
}

// CONTROLLARE SCREENING (anche le nuove cose dello screening parziale) !!!!!!!!!!!!!!!!!!!!
template<typename LATFLOAT>	
typename lattice<LATFLOAT>::Complex	lattice<LATFLOAT>::Ewald_sums(Vector k, int alpha, int beta, int flag_macro)
{
//	cout << "Computing Ewald sums ... " << endl;
	if (k.size()!=3)
		nrerror("lattice::Ewald_sums: wrong input k vector!");
	int i,j,n,m,l;
	LATFLOAT ksq;
	LATFLOAT Gksq,csi,G1,G2;
	int nmax_rec;
	int mmax_rec;
	int lmax_rec;
	int nmax_dir;
	int mmax_dir;
	int lmax_dir;

	LATFLOAT exp_factor=5.0;
	LATFLOAT gamma0=0.2; 
	LATFLOAT cell_volume;
	Complex sum1=0.0, sum2=0.0;
	Complex long_term=0.0;
	Complex I(0.0,1.0);

	vector<Vector>	cell(3);
	vector<Vector>	H(3);
	vector<Vector>	mol_dipole(N_MOL);
	vector<Vector>	rho(N_MOL);
	Vector			R(3); // distance between molecules
	Vector			K(3); // reciprocal lattice vector
	Vector			zero(3);
	zero[0]=0.0;zero[1]=0.0;zero[2]=0.0;

// First we scale cell lengths, dipoles and k vector

	for (i=0; i<3; i++) 
	{
		cell[i].resize(3);
		for (j=0; j<3; j++)
			cell[i][j]=CELL[i][j]/sqrt(epsilon_0(j,j));
	}
	for (i=0;i<N_MOL;i++) 
	{
		rho[i].resize(3);
		for (j=0; j<3; j++)
			rho[i][j]=RHO[i][j]/sqrt(epsilon_0(j,j));
	}

	for (i=0;i<N_MOL;i++) 
	{
		mol_dipole[i].resize(3);
		for (j=0; j<3; j++)
			mol_dipole[i][j]=MOL_DIPOLE[i][j]/(pow(epsilon_0(0,0)*epsilon_0(1,1)*epsilon_0(2,2),LATFLOAT(0.25))
			                              *sqrt(epsilon_0(j,j)));
	}
	for (j=0; j<3; j++)
		k[j]*=sqrt(epsilon_0(j,j));

	cell_volume=inner_prod(cell[0],cross_prod(cell[1],cell[2]));
	H[0]=cross_prod(cell[1],cell[2])*2.0*Pi/cell_volume;
	H[1]=cross_prod(cell[2],cell[0])*2.0*Pi/cell_volume;
	H[2]=cross_prod(cell[0],cell[1])*2.0*Pi/cell_volume;

// Usual Ewald procedure starts here

	nmax_rec=int(2.0*exp_factor*gamma0/norm_2(H[0]))+1;
	mmax_rec=int(2.0*exp_factor*gamma0/norm_2(H[1]))+1;
	lmax_rec=int(2.0*exp_factor*gamma0/norm_2(H[2]))+1;
	nmax_dir=int(exp_factor/gamma0/norm_2(cell[0]))+1;
	mmax_dir=int(exp_factor/gamma0/norm_2(cell[1]))+1;
	lmax_dir=int(exp_factor/gamma0/norm_2(cell[2]))+1;

//	cout << nmax_rec << " " << 	mmax_rec << " " << lmax_rec << " " 
//	<< nmax_dir << " " << mmax_dir << " " <<	lmax_dir << endl; 

	ksq=pow(norm_2(k),2);

// Long-wavelength term (it is included only if macro_flag!=0)

	long_term=0.0;
	if (flag_macro!=0)
	{
		if (norm_2(k)/pow(cell_volume, LATFLOAT(0.3333)) > ZERO_TOL)
			long_term+=4.0*Pi/cell_volume*inner_prod(mol_dipole[alpha],k)*inner_prod(mol_dipole[beta],k)/ksq;
		else
			nrerror("lattice::Ewald_sums: long wave term is non analytical at k=0!");
	}
		
// Sum over the reciprocal lattice

	sum1=0.0;
	for (n=-nmax_rec;n<=nmax_rec;n++)  for (m=-mmax_rec;m<=mmax_rec;m++)  for (l=-lmax_rec;l<=lmax_rec;l++)
	{
		if (n!=0 || m!=0 || l!=0 )
			{
				K=H[0]*n+H[1]*m+H[2]*l;
				Gksq=pow(norm_2(k+K),2);
				sum1+=inner_prod(mol_dipole[alpha],K+k)*inner_prod(mol_dipole[beta],K+k)/Gksq*LATFLOAT(exp(-Gksq/(4.0*gamma0*gamma0)))
					*exp(-I*inner_prod(K,rho[alpha]-rho[beta]));	 // This sign is from my notes
			}
	}
	if (norm_2(k)/pow(cell_volume, LATFLOAT(0.3333)) > ZERO_TOL)
			sum1+=inner_prod(mol_dipole[alpha],k)*inner_prod(mol_dipole[beta],k)/ksq
			      *(exp(-ksq/(4.0*gamma0*gamma0))-1.0);
	sum1*=4.0*Pi/cell_volume;



// Sum over the direct lattice

	sum2=0.0;
	for (n=-nmax_dir;n<=nmax_dir;n++)  for (m=-mmax_dir;m<=mmax_dir;m++)  for (l=-lmax_dir;l<=lmax_dir;l++)
	{
		if (n!=0 || m!=0 || l!=0 || alpha!=beta)
		{
			R=cell[0]*n+cell[1]*m+cell[2]*l+rho[alpha]-rho[beta];  // OK with my notes
			csi=norm_2(R)*gamma0;
			G1=pow(csi,-3)*pow(gamma0,3)*( math::erfc(csi)+2.0/sqrt(Pi)*csi*exp(-csi*csi) );
			G2=4.0/(3.0*sqrt(Pi))*pow(gamma0,3)*exp(-csi*csi);

			sum2+=exp(I*inner_prod(k,R))
				*(
				inner_prod(mol_dipole[alpha],mol_dipole[beta])*G1-LATFLOAT(3.0)*(inner_prod(mol_dipole[alpha],R)*inner_prod(mol_dipole[beta],R)/pow(norm_2(R),2))*(G1+G2)				
				);
		}
	}
	if (alpha==beta) 
		sum2-=4.0/(3.0*sqrt(Pi))*pow(gamma0,3)*inner_prod(mol_dipole[alpha],mol_dipole[beta]);

	if (alpha==0 && beta==3)
	{
		cout << sum1/LATFLOAT(1.602) << endl;
		cout << sum2/LATFLOAT(1.602) << endl;
		cout << long_term/LATFLOAT(1.602) << endl;
	}

	return( (sum1+sum2+long_term) / LATFLOAT(1.602)); // returns energy in eV provided distances are in A and dipoles in Debye
}

template<typename LATFLOAT>	
void lattice<LATFLOAT>::Finite_sums(int flag_2D, string mol_mol_int, LATFLOAT NON_SCR_RADIUS, vector<vector<LATFLOAT> >& NN_int, 
				 LATFLOAT rmax, Vector k, ComplexMatrix& Ltilde)
{
//	cout << "Computing Finite sums ..." << endl;
	unsigned int i,j,alpha, beta;
	int nnmax;
	int non_scr_nn;
	Complex	sum=0.0;
	Complex	I(0.0,1.0);
	Vector V;

	Ltilde.resize(N_MOL,N_MOL);
	for (alpha=0; alpha<N_MOL; alpha++)
		for (beta=0; beta<N_MOL; beta++)
			Ltilde(alpha,beta)=0.0;

	if (rmax>NNrmax || flag_2D!=NNflag_2D)
		compute_NN(rmax, flag_2D);

	for (alpha=0; alpha<N_MOL; alpha++)
	{

// Sets the maximum number of nearest neighbours included in the finite sums
		if (mol_mol_int[0]=='N')
			nnmax=NN_int[alpha].size();
		else
		{
			for (i=1; i<NN[alpha].size(); i++)
				if (norm_2(NN[alpha][i].n)>rmax)
					break;
			nnmax=i;
		}
// Sets the number of non screened nearest neighbours 
		for (i=1; i<NN[alpha].size(); i++)
			if (norm_2(NN[alpha][i].n)>NON_SCR_RADIUS)
					break;
		non_scr_nn=i;

// Non screened interactions
		for (i=1; i<min(non_scr_nn, nnmax); i++)
		{
			beta=NN[alpha][i].alpha;
//			if (flag_2D!=1 || abs(RHO[alpha][2]-RHO[beta][2])<ZERO_TOL)
			if (beta>=alpha)
			{
				V=-NN[alpha][i].n; // OK with Davydov and Philpott
				switch(mol_mol_int[0])
				{
					case 'N':  // Nearest neighbors interaction
						Ltilde(alpha,beta)+=exp(I*inner_prod(k,V))* NN_int[alpha][i];
						break;
					case 'T':  // Non Screened transition charge distribution
						Ltilde(alpha,beta)+=exp(I*inner_prod(k,V))*distributed_charge_interaction(alpha, beta, V) ;
						break;
					default: // Screened dipole interaction
						Ltilde(alpha,beta)+=exp(I*inner_prod(k,V))*dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta], V);
				}
			}
		}

// Screened interactions
		for (i=max(1,non_scr_nn); i<nnmax; i++)
		{
			beta=NN[alpha][i].alpha;
//			if (flag_2D!=1 || abs(RHO[alpha][2]-RHO[beta][2])<ZERO_TOL)
			if (beta>=alpha)
			{
				V=-NN[alpha][i].n; // OK with Davydov and Philpott
				switch(mol_mol_int[0])
				{
					case 'N':  // Nearest neighbors interaction
						Ltilde(alpha,beta)+=exp(I*inner_prod(k,V))* NN_int[alpha][i];
						break;
					case 'T':  // Screened transition charge distribution
						Ltilde(alpha,beta)+=exp(I*inner_prod(k,V))
							              *screened_distributed_charge_interaction(alpha, beta, V, 
										   epsilon_0(0,0), epsilon_0(1,1), epsilon_0(2,2));
						break;
					default: // Screened dipole interaction
						Ltilde(alpha,beta)+=exp(I*inner_prod(k,V))*screened_dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta], V, epsilon_0);
				}
			}
		}
	} // end of for (alpha=0; alpha<N_MOL; alpha++)
	
	for (i=0; i<N_MOL; i++)
		for (j=i+1; j<N_MOL; j++)
			Ltilde(j,i)=conj(Ltilde(i,j));

	return;
}

// This routine DOES NOT exploit inversion symmetry

template<typename LATFLOAT>	
typename lattice<LATFLOAT>::Complex lattice<LATFLOAT>::Correction_to_Ewald_sums(Vector k, int alpha, int beta, LATFLOAT rmax, LATFLOAT NON_SCR_RADIUS)
{
	unsigned int i;
	unsigned int nnmax;
	unsigned int non_scr_nn;
	Complex	I(0.0,1.0);
	Vector V;

	if (rmax>NNrmax || NNflag_2D==1)  // we need 3D sums here
	{
		compute_NN(rmax, 0);
		nnmax=NN[alpha].size();
	}
	else
	{
		for (i=1; i<NN[alpha].size(); i++)
			if (norm_2(NN[alpha][i].n)>rmax)
				break;
		nnmax=i;
	}

// Sets the number of non screened nearest neighbours 
		for (i=1; i<NN[alpha].size(); i++)
			if (norm_2(NN[alpha][i].n)>NON_SCR_RADIUS)
					break;
		non_scr_nn=i;

// Sum over the direct lattice

	Complex	sum=0.0;
	for (i=1; i<nnmax; i++)
	{
		if (NN[alpha][i].alpha==beta)
		{
			V=-NN[alpha][i].n; // OK with Davydov and Philpott
//			V=CELL[0]*n+CELL[1]*m+CELL[2]*l+RHO[alpha]-RHO[beta];
			if (i<non_scr_nn)
				sum+=(exp(I*inner_prod(k,V)))*(distributed_charge_interaction(alpha, beta, V)
				       -dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta], V));
			else
				sum+=(exp(I*inner_prod(k,V)))*(screened_distributed_charge_interaction(alpha, beta, V, 
					epsilon_0(0,0), epsilon_0(1,1), epsilon_0(2,2))-screened_dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta], V, epsilon_0));
		}
	}
	return(sum); // returns energy in eV provided distances are in A and dipoles in Debye
}


// flag_macro is passed  to the routine Ewald_sums
// MODE chooses between Ewald_sums, Finite_sums, or Single_layer_sums
template<typename LATFLOAT>	
void lattice<LATFLOAT>::set_Ltilde_and_u(	Vector kvec, string MODE, string mol_mol_int, LATFLOAT NON_SCR_RADIUS, vector<vector<LATFLOAT> >& NN_int, 
					    LATFLOAT rmax, int flag_macro, ComplexMatrix &Ltilde, ComplexMatrix &udc, vector<LATFLOAT> &Values)
{
	cout << "Computing free exciton states at k=" << kvec << endl;
	int i,j, flag_2D;
	Ltilde.resize(N_MOL,N_MOL);
	udc.resize(N_MOL,N_MOL);
	vector<Complex> a(N_MOL*N_MOL);
	
	switch (MODE[0]) 
	{
				case 'F':
					flag_2D=0;
					Finite_sums(flag_2D, mol_mol_int, NON_SCR_RADIUS, NN_int, rmax, kvec, Ltilde);
					for (i=0; i<N_MOL; i++)	for (j=i; j<N_MOL; j++) a[j*N_MOL+i]=Ltilde(i,j);
					break;
				case 'S':
					flag_2D=1;
					Finite_sums(flag_2D, mol_mol_int, NON_SCR_RADIUS, NN_int, rmax, kvec, Ltilde);
					for (i=0; i<N_MOL; i++)	for (j=i; j<N_MOL; j++)  a[j*N_MOL+i]=Ltilde(i,j);
					break;	
				case 'M':
					for (i=0; i<N_MOL; i++)	for (j=i; j<N_MOL; j++)
					{
						Ltilde(i,j)=Ewald_sums(kvec, i, j, flag_macro) + Correction_to_Ewald_sums(kvec, i, j, rmax, NON_SCR_RADIUS);
						a[j*N_MOL+i]=Ltilde(i,j);
					}
					break;
				default:
					for (i=0; i<N_MOL; i++)	for (j=i; j<N_MOL; j++)
					{
						Ltilde(i,j)=Ewald_sums(kvec, i, j, flag_macro);
						a[j*N_MOL+i]=Ltilde(i,j);
					}
	}

// Completes the Ltilde Hermitian matrix
	for (i=0; i<N_MOL; i++)
		for (j=i+1; j<N_MOL; j++)
		{
			Ltilde(j,i)=conj(Ltilde(i,j));
			a[i*N_MOL+j]=conj(a[j*N_MOL+i]);
		}
			

//	cout << Ltilde << endl;
//	cout << a << endl;
//	getchar();
	heev_cpp('V', 'U', N_MOL, a, Values,'O');

	for (i=0; i<N_MOL; i++)
		for (j=0; j<N_MOL; j++)
		{
			udc(i,j)=a[j*N_MOL+i];
		}

	return;
}


// Computes free exciton bands along given directions

template<typename LATFLOAT>	
void lattice<LATFLOAT>::compute_FE_bands(vector<Vector > &dir, unsigned int npoints, string MODE, 
					  string mol_mol_int, LATFLOAT NON_SCR_RADIUS, vector<vector<LATFLOAT> >& NN_int,
					  LATFLOAT rmax, unsigned int flag_macro, unsigned int n0, unsigned int n1, unsigned int n2, Vector& KLIGHT)
{
	cout << "Computing free exciton bands" << endl;
//	int i,j,m;
	unsigned int i,j,m,n,dc;
	Vector k(3);
	vector<ComplexVector>	band_dipole;
	vector<LATFLOAT>		MyValues(N_MOL);
	ComplexMatrix			udc(N_MOL, N_MOL);
	ComplexMatrix			Ltilde(N_MOL, N_MOL);
	ostringstream ss;
	ss << lattice_name << "_" << MODE << "_" << mol_mol_int << "_" << NON_SCR_RADIUS << "_" 
	   << rmax << "_" << flag_macro << "_FE_bands.txt";
	ofstream band_file(ss.str().c_str());
	
// This part computes bands along specified directions

	if (dir.size()>0)
	{
		for (n=0; n<dir.size(); n++)		
			band_file << dir[n] << "\t" ;
		band_file << endl; 

// This is for k approx 0
		for (n=0; n<dir.size(); n++)		
		{	
			k=dir[n]*0.001;
			set_Ltilde_and_u(k, MODE, mol_mol_int, NON_SCR_RADIUS, NN_int, 
				         rmax, flag_macro, Ltilde, udc, MyValues);
//			band_file << 0.001 << "\t" ;
			band_file << k << "\t" ;
			for (dc=0; dc<MyValues.size(); dc++)
				band_file << MyValues[dc] << "\t";
		}
		band_file << endl;

// This is for k!=0
		for (j=1; j<npoints; j++)
		{
			for (n=0; n<dir.size(); n++)
			{
				k=dir[n]*j/2.0/(npoints-1);
				set_Ltilde_and_u(k, MODE, mol_mol_int, NON_SCR_RADIUS, NN_int, 
				         rmax, flag_macro, Ltilde, udc, MyValues);
						
//				band_file << j/2.0/(npoints-1) << "\t" ;
				band_file << k << "\t" ;
				for (dc=0; dc<MyValues.size(); dc++)
					band_file << MyValues[dc] << "\t";
			}
			band_file << endl;
		}	
	}

// This part computes bands on a user specified 3D grid
	if (n0*n1*n2>0)
	{
		if (n0%2!=0 || n1%2!=0 || n2%2!=0)
			nrerror("compute_FE_bands_3D: grid has an odd number of points!");
		Vector ka, kb, kc;
		Complex I(0.0,1.0);
		for (i=0; i<n0; i++)	
		{
		ka=get_KCELL()[0]*(-n0/2+1+i)/n0;
		for (j=0; j<n1; j++)
		{
		kb=get_KCELL()[1]*(-n1/2+1+j)/n1;
		for (m=0; m<n2; m++)
		{
			kc=get_KCELL()[2]*(-n2/2+1+m)/n2;
			{
				k=ka+kb+kc;
				if (norm_2(k)<ZERO_TOL)	k=KLIGHT;
				band_file << (-n0/2+1+i) << "\t" << (-n1/2+1+j) << "\t" << (-n2/2+1+m) << "\t";
				set_Ltilde_and_u(k, MODE, mol_mol_int, NON_SCR_RADIUS, NN_int, 
			         rmax, flag_macro, Ltilde, udc, MyValues);
				for (dc=0; dc<MyValues.size(); dc++)
					band_file << MyValues[dc] << "\t";
				band_file << endl;
			}
		}}}
	}
	band_file.close();
	
// Here we compute energies, eigenvectors and dipoles at k=KLIGHT
	Complex I(0.0,1.0);
	ss.str("");
	ss << lattice_name << "_" << MODE << "_" << mol_mol_int << "_" << NON_SCR_RADIUS << "_" 
	   << rmax << "_" << flag_macro << "_FE_bands_at_KLIGHT.txt";
	ofstream bands_at_klight(ss.str().c_str());
	k=KLIGHT;
	ComplexVector  dipole(3);
	set_Ltilde_and_u(k, MODE, mol_mol_int, NON_SCR_RADIUS, NN_int, 
			         rmax, flag_macro, Ltilde, udc, MyValues);
	
	bands_at_klight << "---------------------------------------" << endl;
	bands_at_klight << "Bands at k = KLIGHT = " << KLIGHT << endl;
	bands_at_klight << "---------------------------------------" << endl << endl;
	for (dc=0; dc<MyValues.size(); dc++)
	{

		for (j=0; j<3; j++)	
		{
			dipole(j)=0.0;
			for (i=0; i<MyValues.size(); i++)
				dipole(j)+=MOL_DIPOLE[i](j)*udc(i,dc);
		}
		bands_at_klight << 	"Energy: " << MyValues[dc] << "\t Eigenvector: ";
		for (i=0; i<MyValues.size(); i++)
			bands_at_klight << udc(i,dc) << " ";
		bands_at_klight << "\t Dipole at |klight|<<1/a: ";
		for (i=0; i<3; i++)
			bands_at_klight << dipole(i) << " ";
		bands_at_klight << endl;
	}
	bands_at_klight.close();

	return;
}


// Interaction between transition charge distributions 
// r is the distance between the center of mass of the two molecules

template<typename LATFLOAT>
LATFLOAT lattice<LATFLOAT>::distributed_charge_interaction(int alpha, int beta, Vector rba)
{
	if ( norm_2(rba) < (1 - ZERO_TOL)*norm_2(NN[alpha][1].n) )
		return(0.0);
	LATFLOAT sum=0.0;
	unsigned int i,j;
	for (i=0; i<ATOM[alpha].size(); i++)
		for (j=0; j<ATOM[beta].size(); j++)
			sum+=14.397174*CHARGE[alpha][i]*CHARGE[beta][j]/norm_2(ATOM[beta][j]-rba-ATOM[alpha][i]);
	return(sum);
}


template<typename LATFLOAT>
LATFLOAT lattice<LATFLOAT>::screened_distributed_charge_interaction(int alpha, int beta, Vector rba, 
												 LATFLOAT e00, LATFLOAT e11, LATFLOAT e22)
{
	if ( norm_2(rba) < (1 - ZERO_TOL)*norm_2(NN[alpha][1].n) )
		return(0.0);
	LATFLOAT sum=0.0;
	unsigned int i,j;
	LATFLOAT x2,y2,z2;
	for (i=0; i<ATOM[alpha].size(); i++)
		for (j=0; j<ATOM[beta].size(); j++)
		{
			x2=pow(ATOM[beta][j][0]-rba[0]-ATOM[alpha][i][0],2)/e00;
			y2=pow(ATOM[beta][j][1]-rba[1]-ATOM[alpha][i][1],2)/e11;
			z2=pow(ATOM[beta][j][2]-rba[2]-ATOM[alpha][i][2],2)/e22;
			sum+=14.397174*CHARGE[alpha][i]*CHARGE[beta][j]/sqrt(x2+y2+z2);
		}
	return(sum/sqrt(e00*e11*e22));
}

//////////////////////////////////////////////////////////
// This is the explicit instantiation of the template

template class lattice<float>;
template class lattice<double>;
