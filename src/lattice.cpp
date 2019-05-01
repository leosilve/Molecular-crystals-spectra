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

// This function accesses members without get_ functions
template<typename LATFLOAT>
void lattice<LATFLOAT>::initialize(string input_file)
	{ 
		ifstream infile(input_file.c_str());
		if (!infile.is_open()) nrerror("File not found!");
		
		cout << "Initializing lattice from input file: " << input_file << endl;
		
		int i,j,k,n;
		string line;
		Vector				temp(3);
		sym_op<LATFLOAT>	sym_temp;
		Matrix				mat_temp(3,3);	
		LATFLOAT			d_temp;
		int					i_temp;
		pos<LATFLOAT>		pos_temp;

		///////////////////////////////////////////// Initializes lattice_name
		lattice_name=input_file.substr(0,input_file.find('.'));
		
		///////////////////////////////////////////// Initializes CELL and VCELL
		getline(infile, line, ':');
		CELL.resize(3);
		for (i=0;i<3;i++) 
		{
			CELL[i].resize(3);
			for (j=0;j<3;j++) infile >> CELL[i](j);
		}
		VCELL=inner_prod(get_CELL(0),cross_prod(get_CELL(1),get_CELL(2)));
		
		///////////////////////////////////////////// Initializes N_MOL, N_SP and N_SYM_OP
		get_next_from_file(infile, N_MOL);
		get_next_from_file(infile, N_SP);
		get_next_from_file(infile, N_SYM_OP);
		if (get_N_SP()*get_N_SYM_OP()!=get_N_MOL())
			nrerror("lattice<LATFLOAT>::read_crystal_data: Symmetry operations generate more molecules than needed");	

		///////////////////////////////////////////// Initializes SYM
		SYM.resize(0);
		for (k=0;k<get_N_SYM_OP();k++) 
		{
			getline(infile, line, ':');
			for (i=0;i<3;i++) 
				for (j=0;j<3;j++) 
					infile >> mat_temp(i,j);
			sym_temp.set_rot(mat_temp);
			for (j=0;j<3;j++) infile >> temp(j);
			sym_temp.set_ax(get_CELL(0)*temp(0)+get_CELL(1)*temp(1)+get_CELL(2)*temp(2));

			for (j=0;j<3;j++) infile >> temp[j];
			sym_temp.set_tr(get_CELL(0)*temp(0)+get_CELL(1)*temp(1)+get_CELL(2)*temp(2));
			sym_temp.set_map_size(get_N_MOL());
			for (j=0;j<get_N_MOL();j++) 
			{
				infile >> i_temp;
				sym_temp.set_map(j,i_temp);
			}
			SYM.push_back(sym_temp);
		}

		///////////////////////////////////////////// Initializes OtoF
		OtoF.resize(3,3);
		for (i=0;i<3;i++) 
			for (j=0;j<3;j++) 
				OtoF(i,j)=get_CELL(j)(i);
		InvertMatrix(OtoF,OtoF);

		///////////////////////////////////////////// Initializes RHO
		getline(infile, line, ':');
		RHO.resize(0);
		for (i=0;i<get_N_SP();i++) 
		{ 
			for (j=0;j<3;j++) infile >> temp(j);	// Input data are in fractional coordinates
			temp=FracToOrth(temp);					// temp is converted into orthogonal coordinates
			temp=move_into_primitive_cell_orth(temp);
			RHO.push_back(temp);
			for (k=1; k< get_N_MOL()/get_N_SP(); k++)
			{
				pos_temp.alpha=i*get_N_MOL()/get_N_SP();
				pos_temp.n=temp;
				get_SYM(k).apply(pos_temp, pos_temp);
				pos_temp.n=move_into_primitive_cell_orth(pos_temp.n);
				RHO.push_back(pos_temp.n); 
			}
		}	

		///////////////////////////////////////////// Initializes MOL_DIPOLE
		getline(infile, line, ':');
		MOL_DIPOLE.resize(0);
		for (i=0;i<get_N_SP();i++) 
		{
			for (j=0;j<3;j++) infile >> temp(j);
			MOL_DIPOLE.push_back(temp);
			for (k=1; k< get_N_MOL()/get_N_SP(); k++)
				MOL_DIPOLE.push_back(prod(get_SYM(k).get_rot(),temp));  
		}	
		getline(infile, line, ':');
		for (i=0;i<get_N_SP();i++) 
		{
			infile >> d_temp;
			for (k=0; k<get_N_MOL()/get_N_SP(); k++) 
				for (j=0;j<3;j++)
					MOL_DIPOLE[k+i*get_N_MOL()/get_N_SP()](j)*=d_temp;
		};
		
		///////////////////////////////////////////// Initializes ATOM and CHARGE
		ATOM.resize(get_N_MOL());
		CHARGE.resize(get_N_MOL());
		for (i=0;i<get_N_SP();i++) 
		{
			getline(infile, line, ':');
			infile >> i_temp;
			for (k=0; k< get_N_MOL()/get_N_SP(); k++)
			{
					ATOM[k+i*get_N_MOL()/get_N_SP()].resize(i_temp);
					CHARGE[k+i*get_N_MOL()/get_N_SP()].resize(i_temp);
			}
			for (n=0; n<i_temp; n++)
			{
				for (j=0;j<3;j++) infile >> temp(j);
				infile >> d_temp;
				CHARGE[i*get_N_MOL()/get_N_SP()][n]=d_temp;
				ATOM[i*get_N_MOL()/get_N_SP()][n]=temp;
				for (k=1; k< get_N_MOL()/get_N_SP(); k++)
				{
					ATOM[k+i*get_N_MOL()/get_N_SP()][n]=prod(get_SYM(k).get_rot(),temp);
					CHARGE[k+i*get_N_MOL()/get_N_SP()][n]=d_temp;
				}
			}
		}	
 
		///////////////////////////////////////////// Initializes MODE, NON_SCR_RADIUS, RMAX, MOL_MOL_INT
		get_next_string_from_file(infile, MODE);
		get_next_string_from_file(infile, MOL_MOL_INT);
		get_next_from_file(infile, NON_SCR_RADIUS);
		get_next_from_file(infile, RMAX);
		
		
		///////////////////////////////////////////// Initializes KLIGHT			
		getline(infile, line, ':');
		KLIGHT.resize(3);
		for (i=0;i<3;i++) 
			infile >> KLIGHT[i];
		
		///////////////////////////////////////////// Initializes epsilon_0		
		getline(infile, line, ':');
		epsilon_0.resize(3,3);
		for (i=0;i<3;i++) 
			for (j=0;j<3;j++) 
				infile >> epsilon_0(i,j);
		
		///////////////////////////////////////////// Initializes NN_int
		getline(infile, line, ':');
		infile >> i_temp;
		if (i_temp!=get_N_MOL())
		{
			cout << "No valid NN_interactions, setting NN_int size to 0!" << endl;
			NN_int.resize(0);
			getline(infile, line, ':');
			getline(infile, line, ':');
		}
		else
		{
			NN_int.resize(get_N_MOL());
			getline(infile, line, ':');
			infile >> i_temp;
			for (i=0; i<get_N_MOL(); i++)
				 NN_int[i].resize(i_temp);
			getline(infile, line, ':');
			for (i=0; i<get_N_MOL(); i++)
				for (j=0; j<NN_int[i].size(); j++)
					infile >> NN_int[i][j];
		}	
		
		infile.close();
		cout << "Finished reading!" << endl;
		
		///////////////////////////////////////////// Initializes KCELL
		KCELL.resize(3);
		KCELL[0]=cross_prod(get_CELL(1),get_CELL(2))*(2.0*Pi/get_VCELL());
		KCELL[1]=cross_prod(get_CELL(2),get_CELL(0))*(2.0*Pi/get_VCELL());
		KCELL[2]=cross_prod(get_CELL(0),get_CELL(1))*(2.0*Pi/get_VCELL());

		///////////////////////////////////////////// Initializes NN
		int flag_2D;
		if (MODE[0]=='S')
			flag_2D=1;
		else
			flag_2D=0;
		NN=compute_NN(RMAX,flag_2D);
		
		///////////////////////////////////////////// Prints lattice.log
/*		string filename = lattice_name + ".log";
		ofstream log_file(filename.c_str());
		print(log_file);
		log_file.close();*/
		return;
	};
	
// Computes the nearest neighbours for each inequivalent molecule
// Check if it exploits all the symmetry operations 
// If flag_2D==1 it considers only NN on the same xy plane
// This function accesses NN without get_ functions
template<typename LATFLOAT>	
std::vector<std::vector<pos<LATFLOAT> > > lattice<LATFLOAT>::compute_NN(LATFLOAT rmax, int flag_2D)
	{
		cout << "Computing NN with rmax=" << rmax << " and flag_2D=" << flag_2D << endl;
		int alpha;
		int i,j,k,n,m,l;
		int max_a, max_b,max_c;
		int sp, mol;
		LATFLOAT PLANAR_TOL=norm_2(get_CELL(2))/2.0;
		pos<LATFLOAT> position;
		typename vector<pos<LATFLOAT> >::iterator it;
		std::vector<std::vector<pos<LATFLOAT> > > new_NN(get_N_MOL(), vector<pos<LATFLOAT> >(0));

		if (rmax<1e-3)
			rmax=1e-3; 

/*		NN.resize(get_N_MOL());
		for (i=0; i<get_N_MOL(); i++)
			NN[i].resize(0);*/

		max_a=int(rmax/norm_2(get_CELL(0)))+1;
		max_b=int(rmax/norm_2(get_CELL(1)))+1;
		if (flag_2D!=1)
			max_c=int(rmax/norm_2(get_CELL(2)))+1;
		else
			max_c=1;

		for (sp=0;sp<get_N_SP();sp++) 
		{
			mol=sp*get_N_MOL()/get_N_SP();

	// Loop over all the molecules in the crystal
			for (alpha=0; alpha<get_N_MOL(); alpha++) 
				for (n=-max_a; n<=max_a; n+=1) 
					for (m=-max_b; m<=max_b; m+=1)  
						for (l=-max_c; l<=max_c; l+=1) 
							{
								position.alpha=alpha;
								position.n=n*get_CELL(0)+m*get_CELL(1)+l*get_CELL(2)+get_RHO(alpha)-get_RHO(mol);
								if (flag_2D!=1 || abs(position.n[2])<PLANAR_TOL)
								if ( norm_2(position.n)<rmax )
								{
									for (it=new_NN[mol].begin(); it!=new_NN[mol].end(); it++ ) 
										if ( norm_2((*it).n)>norm_2(position.n) ) break;
									new_NN[mol].insert(it,position);
								}
							}

// Fixes other molecules basis states so that symmetry operations leave the phonon cloud the same
			for (k=1; k< get_N_MOL()/get_N_SP(); k++)
			{	
				n=get_SYM(k).get_map()[0];
				new_NN[mol+n].resize(new_NN[mol].size());
				for (j=0; j<new_NN[mol+n].size(); j++)	
				{
					new_NN[mol+n][j].alpha=get_SYM(k).get_map(new_NN[mol][j].alpha);
					new_NN[mol+n][j].n=prod(get_SYM(k).get_rot(),new_NN[mol][j].n);
				}
			}
		} // end loop over species	
		return(new_NN);
	}

template<typename LATFLOAT>	
void lattice<LATFLOAT>::print(ofstream& log_file)
	{
		cout << "Printing lattice data ... " << endl;

		unsigned int i,j;
		log_file << "==================================================================================================" << endl;
		log_file << "     Lattice  " << lattice_name << endl;
		log_file << "==================================================================================================" << endl;
		log_file << endl;		
		log_file << "Static epsilon: " << epsilon_0 << endl;
		log_file << "--------------------------- Crystal cell " << endl; 
		for (i=0; i<3; i++)
			log_file << get_CELL(i) << endl;
		log_file << "Cell volume: " << get_VCELL() << " A^3" << endl;
		log_file << endl << "--------------------------- Crystal reciprocal lattice cell " << endl; 
		for (i=0; i<3; i++)
			log_file << get_KCELL(i) << endl;
		log_file << endl << "--------------------------- OtoF (transformation matrix from Orthogonal to Fractional coordinates) " << endl; 
		log_file << get_OtoF() << endl ;
		log_file << endl << "--------------------------- Molecule position inside the cell (orthogonal coordinates)" << endl;
		for (i=0; i<get_N_MOL(); i++)
			log_file << get_RHO(i) << endl;
		log_file << endl << "--------------------------- Molecular dipole moments (in Debyes)" << endl;
		for (i=0; i<get_N_MOL(); i++)
			log_file << get_MOL_DIPOLE(i) << endl;
		log_file << endl << "------------- Atomic positions with respect to the CoM (orthogonal coordinates in A) and charges (in e)" << endl;
		for (i=0; i<get_N_MOL(); i++)
		{
			log_file << "MOLECULE: " << i << endl; 
			for (j=0; j<get_ATOM_size(i); j++)
				log_file << get_ATOM(i,j) << "\t" <<  get_CHARGE(i,j) << endl;
		}
		log_file << endl << "--------------------------- Symmetry operations " << endl;
		for (i=0; i<get_N_SYM_OP(); i++)
			log_file << get_SYM(i) << endl;
		log_file << endl << "--------------------------- Nearest neighbours " << endl;
		for (i=0; i<get_NN_size(); i++)
		{
			log_file << "MOLECULE: " << i << endl; 
			for (j=1; j<get_NN_size(i); j++)
			{
				log_file << j << ": "  << get_NN(i,j) 
				<< "\tT: " << distributed_charge_interaction(i, NN[i][j].alpha, -NN[i][j].n) 
				<< "\tD: " << dipole_dipole(MOL_DIPOLE[i], MOL_DIPOLE[NN[i][j].alpha], NN[i][j].n)
				<< "\tST: " << screened_distributed_charge_interaction(i, NN[i][j].alpha, -NN[i][j].n, epsilon_0(0,0), epsilon_0(1,1), epsilon_0(2,2)) 
				<< "\tSD: " << screened_dipole_dipole(MOL_DIPOLE[i], MOL_DIPOLE[NN[i][j].alpha], NN[i][j].n, epsilon_0)				<< endl;
			}
		}
		return;
	}

////////////////////////////////////////////////////////////////////////////////////
// Functions to set parameters
/*
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
	}*/

////////////////////////////////////////////////////////////////////////////////////
// Interactions

template<typename LATFLOAT>	
LATFLOAT lattice<LATFLOAT>::interaction( int alpha, int beta, Vector dr)
{
		switch(MOL_MOL_INT[0])
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
LATFLOAT lattice<LATFLOAT>::interaction( int alpha, int nn)
{
	if (nn<NN[alpha].size())
	{
		Vector dr=NN[alpha][nn].n;
		int beta=NN[alpha][nn].alpha;
		switch(MOL_MOL_INT[0])
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
	int i,j,n,m,l;
	if (k.size()!=3)
		nrerror("lattice::Ewald_sums: wrong input k vector!");

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
	LATFLOAT vcell;
	Complex sum1=0.0, sum2=0.0;
	Complex long_term=0.0;
	Complex I(0.0,1.0);

	vector<Vector>	cell(3);
	vector<Vector>	H(3);
	vector<Vector>	mol_dipole(N_MOL);
	vector<Vector>	rho(N_MOL);
	Vector			R(3); // distance between molecules
	Vector			K(3); // reciprocal lattice vector
//	Vector			zero(3);
//	zero(0)=0.0;zero(1)=0.0;zero(2)=0.0;

// First we scale cell lengths, dipoles and k vector
	for (i=0; i<3; i++) 
	{
		cell[i].resize(3);
//		cell[i]=zero;
		for (j=0; j<3; j++)
			cell[i](j)=get_CELL(i)(j)/sqrt(epsilon_0(j,j));
	}
	
	for (i=0;i<N_MOL;i++) 
	{
		rho[i].resize(3);
//		rho[i]=zero;
		for (j=0; j<3; j++)
			rho[i](j)=get_RHO(i)(j)/sqrt(epsilon_0(j,j));
	}

	for (i=0;i<N_MOL;i++) 
	{
		mol_dipole[i].resize(3);
//		mol_dipole[i]=zero;
		for (j=0; j<3; j++)
			mol_dipole[i](j)=get_MOL_DIPOLE(i)(j)/(pow(epsilon_0(0,0)*epsilon_0(1,1)*epsilon_0(2,2),LATFLOAT(0.25))
			                              *sqrt(epsilon_0(j,j)));
	}
	
	for (j=0; j<3; j++)
		k(j)*=sqrt(epsilon_0(j,j));

	vcell=inner_prod(cell[0],cross_prod(cell[1],cell[2]));
	H[0]=cross_prod(cell[1],cell[2])*2.0*Pi/vcell;
	H[1]=cross_prod(cell[2],cell[0])*2.0*Pi/vcell;
	H[2]=cross_prod(cell[0],cell[1])*2.0*Pi/vcell;

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
		if (norm_2(k)/pow(vcell, LATFLOAT(0.3333)) > ZERO_TOL)
			long_term+=4.0*Pi/vcell*inner_prod(mol_dipole[alpha],k)*inner_prod(mol_dipole[beta],k)/ksq;
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
	if (norm_2(k)/pow(vcell, LATFLOAT(0.3333)) > ZERO_TOL)
			sum1+=inner_prod(mol_dipole[alpha],k)*inner_prod(mol_dipole[beta],k)/ksq
			      *(exp(-ksq/(4.0*gamma0*gamma0))-1.0);
	sum1*=4.0*Pi/vcell;

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
				inner_prod(mol_dipole[alpha],mol_dipole[beta])*G1-LATFLOAT(3.0)*(inner_prod(mol_dipole[alpha],R)*inner_prod(mol_dipole[beta],R)/pow(norm_2(R),LATFLOAT(2.0)))*(G1+G2)
				) ;
		}
	}
	if (alpha==beta) 
		sum2-=4.0/(3.0*sqrt(Pi))*pow(gamma0,3)*inner_prod(mol_dipole[alpha],mol_dipole[beta]);

/*	if (alpha==0 && beta==3)
	{
		cout << sum1/LATFLOAT(1.602) << endl;
		cout << sum2/LATFLOAT(1.602) << endl;
		cout << long_term/LATFLOAT(1.602) << endl;
	}*/

//	cout << "done" << endl;
	return( (sum1+sum2+long_term) / LATFLOAT(1.602)); // returns energy in eV provided distances are in A and dipoles in Debye
}

template<typename LATFLOAT>	
typename lattice<LATFLOAT>::Complex lattice<LATFLOAT>::Finite_sums(int flag_2D, Vector k, int alpha, int beta)
{
//	cout << "Computing Finite sums ..." << endl;
//	cout << k << endl;
	unsigned int i,j;
	int nnmax;
	int non_scr_nn;
	Complex	sum=0.0;
	Complex	I(0.0,1.0);
	Vector V(3);

//	for (alpha=0; alpha<N_MOL; alpha++)
//	{

// Sets the maximum number of nearest neighbours included in the finite sums
		if (MOL_MOL_INT[0]=='N')
			nnmax=NN_int[alpha].size();
		else
		{
			for (i=1; i<NN[alpha].size(); i++)
				if (norm_2(NN[alpha][i].n)>RMAX)
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
//			if (flag_2D!=1 || abs(RHO[alpha][2]-RHO[beta][2])<ZERO_TOL)
			if (NN[alpha][i].alpha==beta)
			{
				V=-NN[alpha][i].n; // OK with Davydov and Philpott
				switch(MOL_MOL_INT[0])
				{
					case 'N':  // Nearest neighbors interaction
						sum=exp(I*inner_prod(k,V))* NN_int[alpha][i];
						break;
					case 'T':  // Non Screened transition charge distribution
						sum+=exp(I*inner_prod(k,V))*distributed_charge_interaction(alpha, beta, V) ;
						break;
					default: // Screened dipole interaction
						sum+=exp(I*inner_prod(k,V))*dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta], V);
				}
			}
		}

// Screened interactions
		for (i=max(1,non_scr_nn); i<nnmax; i++)
		{
//			if (flag_2D!=1 || abs(RHO[alpha][2]-RHO[beta][2])<ZERO_TOL)
			if (NN[alpha][i].alpha==beta)
			{
				V=-NN[alpha][i].n; // OK with Davydov and Philpott
				switch(MOL_MOL_INT[0])
				{
					case 'N':  // Nearest neighbors interaction
						sum+=exp(I*inner_prod(k,V))* NN_int[alpha][i];
						break;
					case 'T':  // Screened transition charge distribution
						sum+=exp(I*inner_prod(k,V))
							              *screened_distributed_charge_interaction(alpha, beta, V, 
										   epsilon_0(0,0), epsilon_0(1,1), epsilon_0(2,2));
						break;
					default: // Screened dipole interaction
						sum+=exp(I*inner_prod(k,V))*screened_dipole_dipole(MOL_DIPOLE[alpha], MOL_DIPOLE[beta], V, epsilon_0);
				}
			}
		}
//	} // end of for (alpha=0; alpha<N_MOL; alpha++)
//	cout << sum << endl;
	return(sum);
}

// This routine DOES NOT exploit inversion symmetry

template<typename LATFLOAT>	
typename lattice<LATFLOAT>::Complex lattice<LATFLOAT>::Correction_to_Ewald_sums(Vector k, int alpha, int beta)
{
	unsigned int i;
	unsigned int nnmax;
	unsigned int non_scr_nn;
	Complex	I(0.0,1.0);
	Vector V;

/*	if (RMAX>NN_RMAX || NN_flag2D==1)  // we need 3D sums here
	{
		compute_NN(RMAX, 0);
		nnmax=NN[alpha].size();
	}
	else
	{
		for (i=1; i<NN[alpha].size(); i++)
			if (norm_2(NN[alpha][i].n)>RMAX)
				break;
		nnmax=i;
	}*/
	
	nnmax=get_NN_size(alpha);

// Sets the number of non screened nearest neighbours 
	for (i=1; i<nnmax; i++)
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
//			V=get_CELL(0)*n+get_CELL(1)*m+get_CELL(2)*l+RHO[alpha]-RHO[beta];
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

template<typename LATFLOAT>	
typename lattice<LATFLOAT>::Complex lattice<LATFLOAT>::compute_Ltilde(Vector kvec, int i, int j, int flag_macro)
{
//	cout << "Computing Ltilde at k=" << kvec << endl;
	int flag_2D;
	
	switch (MODE[0]) 
	{
		case 'F':
			flag_2D=0;
			return(Finite_sums(flag_2D, kvec, i,j));
			break;
		case 'S':
			flag_2D=1;
			return(Finite_sums(flag_2D, kvec, i,j));
			break;	
		case 'M':
			if (norm_2(kvec)<ZERO_TOL && flag_macro==1)	
			{
				kvec=get_KLIGHT();
				cout << "lattice<LATFLOAT>::compute_Ltilde: Set kvec=KLIGHT!" << endl;
			}
			return(Ewald_sums(kvec, i, j, flag_macro) + Correction_to_Ewald_sums(kvec, i, j));
			
			break;
		default:
			if (norm_2(kvec)<ZERO_TOL && flag_macro==1)	
			{
				kvec=get_KLIGHT();
				cout << "lattice<LATFLOAT>::compute_Ltilde: Set kvec=KLIGHT!" << endl;
			}
			return(Ewald_sums(kvec, i, j, flag_macro));
	}
	return(0.0);
}

template<typename MYFLOAT>
void lattice<MYFLOAT>::computes_lattice_bands(string& band_input_file)
{
	ifstream infile(band_input_file.c_str());
	if (!infile.is_open()) nrerror("computes_lattice_bands: input file not found!");
	
	istringstream  iss;
	string line;
	vector<MYFLOAT> dir_temp;
	vector<Vector> dir;
	int flag_macro=1;
	int npoints=6;
	int n0,n1,n2;
	int i_temp,i,j,d;
	
	// Reads directions
	getline(infile, line, ':');
	getline(infile, line);
	iss.str( line );
	copy( istream_iterator<MYFLOAT>( iss ), istream_iterator<MYFLOAT>(), back_inserter( dir_temp ) );
	if (dir_temp.size()%3 !=0)
		nrerror("computes_lattice_bands: Wrong input kdirs!");
	dir.resize(dir_temp.size()/3);
	for (d=0; d<dir.size(); d++)
	{
		dir[d].resize(3);
		for (i=0; i<3; i++)
			dir[d](i)=dir_temp[d*3+i];
	}
	iss.clear();

	get_next_from_file(infile, flag_macro);
	get_next_from_file(infile, npoints);
	getline(infile, line, ':');
	infile >> n0;
	infile >> n1;
	infile >> n2;
	
	compute_FE_bands(dir, npoints, n0, n1, n2, flag_macro);
	return;
}
// Computes free exciton bands 
template<typename LATFLOAT>	
void lattice<LATFLOAT>::compute_FE_bands(vector<Vector> dir, int npoints, int n0, int n1, int n2, int flag_macro)
{
	int i,j,m,n,dc,row,col;
	int dim=get_N_MOL();
	LATFLOAT kscale=1.0;
	Vector k(3);
	vector<LATFLOAT>		EigVal(dim);
	vector<Complex>			EigVec(dim*dim,0);
/*	for (i=0; i<dim; i++)
		EigVec[i].resize(dim);
	for (row=0; row<dim; row++)
		for (col=0; col<dim; col++)
			EigVec[col][row]=0.0;*/
	ostringstream ss;
	
// This part computes bands along specified directions

	if (dir.size()>0)
	{
		cout << "Computing free exciton bands along directions: ";
		for (n=0; n<dir.size(); n++)		
			cout << dir[n] << "\t" ;
		
		ss << lattice_name << "_" << MODE << "_" << MOL_MOL_INT << "_" << NON_SCR_RADIUS << "_" 
		<< RMAX << "_" << flag_macro << "_bands_dir.txt";
		ofstream band_file(ss.str().c_str());
		band_file << "# ";
		for (n=0; n<dir.size(); n++)		
			band_file << dir[n] << "\t" ;
		band_file << endl; 

		for (j=1; j<npoints; j++)
		{
			for (n=0; n<dir.size(); n++)
			{
				if (j==0)
					kscale=0.001;
				else
					kscale=LATFLOAT(j)/2.0/LATFLOAT(npoints-1);
				k=dir[n]*kscale;
				for (row=0; row<dim; row++)
					for (col=row; col<dim; col++)
						EigVec[col*dim+row]=compute_Ltilde(k, row, col, flag_macro);
				heev_cpp('V','U',dim,EigVec,EigVal,'O');
				band_file << kscale << "\t" ;
//				band_file << k << "\t" ;
				for (dc=0; dc<EigVal.size(); dc++)
					band_file << EigVal[dc] << "\t";
			}
			band_file << endl;
		}
		band_file.close();
	}
	cout << endl;

// This part computes bands on a user specified 3D grid
	if (n0*n1*n2>0)
	{
		if (n0%2!=0 || n1%2!=0 || n2%2!=0)
			nrerror("compute_FE_bands_3D: grid has an odd number of points!");
		cout << "Computing free exciton bands on a grid " << n0 << "x" << n1 << "x" << n2 << endl;
		ss.str("");
		ss << lattice_name << "_" << MODE << "_" << MOL_MOL_INT << "_" << NON_SCR_RADIUS << "_" 
		<< RMAX << "_" << flag_macro << "_bands_grid.txt";
		ofstream bands_grid(ss.str().c_str());
		Vector ka(3), kb(3), kc(3);
		Complex I(0.0,1.0);
		for (i=0; i<n0; i++)	
		{
			ka=get_KCELL(0)*(-n0/2+1+i)/n0;
			for (j=0; j<n1; j++)
			{
				kb=get_KCELL(1)*(-n1/2+1+j)/n1;
				for (m=0; m<n2; m++)
				{
					kc=get_KCELL(2)*(-n2/2+1+m)/n2;
					k=ka+kb+kc;
					
					if (norm_2(k)<ZERO_TOL)	k=KLIGHT;
					bands_grid << (-n0/2+1+i) << "\t" << (-n1/2+1+j) << "\t" << (-n2/2+1+m) << "\t";
					for (row=0; row<dim; row++)
						for (col=row; col<dim; col++)
							EigVec[col*dim+row]=compute_Ltilde(k, row, col, flag_macro);
					heev_cpp('V','U',dim,EigVec,EigVal,'O');
					for (dc=0; dc<EigVal.size(); dc++)
						bands_grid << EigVal[dc] << "\t";
					bands_grid << endl;
				}
			}
		}
		bands_grid.close();
	}

	
// Here we compute energies, eigenvectors and dipoles at k=KLIGHT
	ss.str("");
	ss << lattice_name << "_" << MODE << "_" << MOL_MOL_INT << "_" << NON_SCR_RADIUS << "_" 
	   << RMAX << "_" << flag_macro << "_bands_KLIGHT.txt";
	ofstream bands_at_klight(ss.str().c_str());
	ComplexVector  dipole(3);
	k=KLIGHT;
	for (row=0; row<dim; row++)
		for (col=row; col<dim; col++)
			EigVec[col*dim+row]=compute_Ltilde(k, row, col, flag_macro);
//	cout << EigVec << endl;
	heev_cpp('V','U',dim,EigVec,EigVal,'O');
	bands_at_klight << "---------------------------------------" << endl;
	bands_at_klight << "Bands at k = KLIGHT = " << k << endl;
	bands_at_klight << "---------------------------------------" << endl << endl;
//	bands_at_klight << 	"Energy (eV)\tdx\tdy\tdz (D) " << endl;
	for (dc=0; dc<dim; dc++)
	{
		bands_at_klight << 	"En: " << EigVal[dc] << "\t EigVec: ";
//		bands_at_klight << 	EigVal[dc] << "\t";
		for (i=0; i<dim; i++)
			bands_at_klight << EigVec[dc*dim+i] << " ";
		for (j=0; j<3; j++)	
		{
			dipole(j)=0.0;
			for (i=0; i<dim; i++)
				dipole(j)+=get_MOL_DIPOLE(i)(j)*EigVec[dc*dim+i];
		}

		bands_at_klight << "\t Dipole: ";
		for (i=0; i<3; i++)
			bands_at_klight << dipole(i) << " ";
/*		{
			if (abs(dipole(i))>0.00001)
				bands_at_klight << dipole(i).real() << "\t";
			else
				bands_at_klight << 0 << "\t";
		}*/
			
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
