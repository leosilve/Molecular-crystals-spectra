/* Written by Leonardo Silvestri 2007-2013 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <complex>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>
#include <OPutils.h>
#include <N_Functions.h> 
#include <ME.h>
#include <physics.h>
#include <lattice.h>
#include <QuantumState.h>
#include <SymmetricDoubleWell.h>
#include <PhononCloud.h>
 
namespace ublas = boost::numeric::ublas;
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

template <typename FLOAT>
void PhononCloud<FLOAT>::initialize(string input_file, string lattice_file)
{
		prob_name=input_file.substr(0,input_file.find('.'));
		
		latticePtr = new lattice<FLOAT>(lattice_file);

		cout << "Reading PC problem data from file: " << input_file << endl;

		int i,j;
		string line;
		istringstream  iss;
		ifstream infile(input_file.c_str());
		if (!infile.is_open()) nrerror("File not found!");
		FLOAT		d_temp;
		FLOAT		rmax_CLOUD;
		int			i_temp;
		string DW_datafile;

		get_next_from_file(infile, elec_En);
		get_next_from_file(infile, NMODES);
		if (NMODES<1 || NMODES>2)
			nrerror("PhononCloud<FLOAT>::initialize: NMODES must be 1 or 2!");

		get_next_from_file(infile, lambda_0);
		if (NMODES==2)
			infile >> lambda_1;

		get_next_from_file(infile, vib_En_0);
		if (NMODES==2)
			infile >> vib_En_1;
	
		MAX_VIB.resize(NMODES);	
		get_next_from_file(infile, MAX_VIB[0]);
		if (NMODES==2)
			infile >> MAX_VIB[1];

		get_next_from_file(infile, FLAG_2D);
		get_next_from_file(infile, CLOUD_RADIUS);
		
// Mettere un controllo anche su questi parametri?????
		get_next_from_file(infile, NP);
		if (NP<=0)
			nrerror("PhononCloud<FLOAT>::initialize: NP must be >0!");
		if (NP==1)
			CLOUD_RADIUS=0.0;
	
		get_next_from_file(infile, EXP_FAC);
		
		get_next_from_file(infile, FLAG_DW);
		if (FLAG_DW==1)
		{
			getline(infile, line, '"');
			getline(infile, DW_datafile, '"');		
		}
 
/////////////////////////////////////////////////////////////////////////////
// Here starts the real initialization
		cout << "Initialization ... " << endl;
	
		KVEC.resize(3);
		for (i=0; i<3; i++)
			KVEC(i)=0.0;
		
		CLOUD_NN=(*latticePtr).compute_NN(CLOUD_RADIUS, FLAG_2D);
		i_temp=CLOUD_NN[0].size();
		for (i=1; i<(*latticePtr).get_N_MOL(); i++)
			if (CLOUD_NN[i].size()!=i_temp)
				nrerror("PhononCloud::initialize: Clouds have different sizes for nonequivalent molecules!");
	
		if (CLOUD_NN.size()==0)
			nrerror("PhononCloud::initialize: CLOUD_NN.size()==0!");
		else
			MAX_CLOUD=CLOUD_NN[0].size();

		MAX_VIBS.resize(get_NMODES());
		for (i=0; i<get_NMODES(); i++)
		{
			MAX_VIBS[i].resize(get_MAX_CLOUD());
			MAX_VIBS[i][0]=get_MAX_VIB(i);
			for (j=1; j<get_MAX_CLOUD(); j+=2)
				MAX_VIBS[i][j]=max(int(get_MAX_VIB(i)*exp(-get_EXP_FAC()*norm_2(CLOUD_NN[0][j].n)/FLOAT(pow((*latticePtr).get_VCELL(),FLOAT(0.333))))) , 1);
			for (j=2; j<get_MAX_CLOUD(); j+=2)
				MAX_VIBS[i][j]=MAX_VIBS[i][j-1];
		}

// Initializes Symmetric Double Well Parameters
		if (get_FLAG_DW()==1)
			initialize_SDW_parameters(DW_datafile, DW_hw0, DW_alpha, DW_A, DW_Vd,
							   DW_hwe, DW_lambdae, DW_basis_dim, DW_nug_max, DW_nue_max, DW_FCF_table, DW_levels);

		compute_BASIS();
		BD_BASIS.resize(0);
		return;
};

template <typename FLOAT>
void PhononCloud<FLOAT>::print(ofstream& log)
{
	// Writes log file
		cout << "Printing PhononCloud Basis data ... " << endl;
		int i,j;
		log << endl;
		log << "==================================================================================================" << endl;
		log << "PhononCloud problem: " << prob_name << endl;
		log << "==================================================================================================" << endl << endl;
		log << "NMODES : " << get_NMODES() << endl;
		log << "lambda : " << get_lambda_0();
		if (NMODES==2)
			log << "\t" << get_lambda_1();
		log << endl;
		log << "vib_En : " << get_vib_En_0();
		if (NMODES==2)
			log << "\t" << get_vib_En_1();
		log << endl;
		log << "CLOUD_RADIUS : " << get_CLOUD_RADIUS() << endl;
		log << "MAX_CLOUD : " << get_MAX_CLOUD() << endl;
		log << "MAX_VIB : " << get_MAX_VIB() << endl;
		log << "MAX_VIBS : " << get_MAX_VIBS() << endl;
		log << "NP : " << get_NP() << endl;
		log << "k : " << get_KVEC() << endl;
		log << "FLAG_DW : " << get_FLAG_DW() << endl;
		log << "Number of basis states: " << get_BASIS_size() << endl;
		for (i=0; i<get_BASIS_size(); i++)
			log << "State n." << i << " : " << get_BASIS(i)	<< endl;
		if (get_BD_BASIS_size()>0)
		{
			log << endl << "Number of symmetric BD_BASIS states: " << get_BD_BASIS_size(0) << endl;
			for (i=0; i<get_BD_BASIS_size(0) ; i++)
				log << "State n." << i << " : " << get_BD_BASIS(0,i) << endl;
			log << endl << "Number of antisymmetric BD_BASIS states: " << get_BD_BASIS_size(1) << endl;
			for (i=0; i<get_BD_BASIS_size(1) ; i++)
				log << "State n." << i << " : " << get_BD_BASIS(1,i)<< endl;
		}
		log << endl << "--------------------------- Cloud Nearest neighbours " << endl;
		for (i=0; i<CLOUD_NN.size(); i++)
		{
			log << "MOLECULE: " << i << endl; 
			for	(j=1; j<CLOUD_NN[i].size(); j++)
				log << j << ": "  << CLOUD_NN[i][j] << endl;
		}
		return;
}

template <typename FLOAT>
bool PhononCloud<FLOAT>::change_KVEC(Vector new_KVEC)
{
	if (new_KVEC.size()==3)
	{
		cout << "New KVEC = " << new_KVEC << endl; 
		KVEC=new_KVEC;
		compute_BASIS();
		BD_BASIS.resize(0);
	}
	else
		return(false);
	return(true);
}

template <typename FLOAT>
bool PhononCloud<FLOAT>::change_KVEC(FLOAT k0, FLOAT k1, FLOAT k2)
{
	KVEC(0)=k0;
	KVEC(1)=k1;
	KVEC(2)=k2;
	cout << "New KVEC = " << KVEC << endl;
	compute_BASIS();
	BD_BASIS.resize(0);
	return(true);
}

template <typename FLOAT>
void PhononCloud<FLOAT>::set_n_vib(int cloud_pos, int v0_max, int v1_max, ME<FLOAT>& exc)
{
	int i, j, np; 
	if ( cloud_pos == get_MAX_CLOUD() - 1 )
	{	
		if (v0_max <= get_MAX_VIBS(0,cloud_pos) && v1_max <= get_MAX_VIBS(1,cloud_pos))
		{
			exc.set_cloud(0, cloud_pos, v0_max);
			exc.set_cloud(1, cloud_pos, v1_max);
//      Checks number of particles
			if (exc.get_np()<=get_NP())
					BASIS.push_back(exc);
		} // end if v_max
		return;
	}  // end if cloud_pos == MAX_CLOUD - 1 
	else
	{
		for (i=0; i<=min(v0_max, get_MAX_VIBS(0,cloud_pos)); i++)
			for (j=0; j<=min(v1_max, get_MAX_VIBS(1,cloud_pos)); j++)
		{
			exc.set_cloud(0, cloud_pos,i);
			exc.set_cloud(1, cloud_pos,j);
			set_n_vib(cloud_pos+1, v0_max-i, v1_max-j, exc);
		}
		return;
	}
} 

template <typename FLOAT>
void PhononCloud<FLOAT>::set_n_vib_one_mode(int cloud_pos, int v0_max, ME<FLOAT>& exc)
{
	int i, np; 
	if ( cloud_pos == get_MAX_CLOUD() - 1 )
	{	
		if (v0_max <= get_MAX_VIBS(0,cloud_pos))
		{
			exc.set_cloud(0, cloud_pos, v0_max);
//      Checks number of particles
			if (exc.get_np()<=get_NP())
					BASIS.push_back(exc);
		} // end if v_max
		return;
	}  // end if cloud_pos == MAX_CLOUD - 1 
	else
	{
		for (i=0; i<=min(v0_max, get_MAX_VIBS(0,cloud_pos)); i++)
		{
			exc.set_cloud(0, cloud_pos,i);
			set_n_vib_one_mode(cloud_pos+1, v0_max-i, exc);
		}
		return;
	}
} 

template <typename FLOAT>
ME<FLOAT> PhononCloud<FLOAT>::compute_symmetric_state(int sym, ME<FLOAT>& exc)
{
	ME<FLOAT> temp;
	temp.set_nmodes(get_NMODES());
	temp.set_max_cloud(get_MAX_CLOUD());
	temp.set_k_3d(get_KVEC());
	int m,j;
	for (m=0; m<get_NMODES(); m++)
			for (j=0; j<get_MAX_CLOUD(); j++)
				temp.set_cloud(m,j,0);
	if (sym < (*latticePtr).get_N_SYM_OP())
	{
		Vector dr(3);
		Matrix invrot(3,3);
		int newmol=(*latticePtr).get_SYM(sym).get_map(exc.get_mol());
		temp.set_mol(newmol);
		InvertMatrix((*latticePtr).get_SYM(sym).get_rot(),invrot);
		for (m=0; m<get_NMODES(); m++)
			for (j=0; j<get_MAX_CLOUD(); j++)
			{
				dr=prod(invrot,CLOUD_NN[newmol][j].n);
				temp.set_cloud(m,j,exc.get_vib(m, latticePtr, dr));
			}
	}
	else
		nrerror("compute_symmetric_state: called with an invalid sym!");
	return(temp);
}

template <typename FLOAT>
void PhononCloud<FLOAT>::compute_BASIS()
{
	int nmol,i,j,n;
	int n_states=0;
	int n_odd=0;
	int n_exc=0;

	cout << "Computing BASIS set ..." << endl;

	BASIS.resize(0);

	ME<FLOAT> exc;
	nmol=0;
	int mode=0;
	int cloud_pos=0;
	exc.set_nmodes(get_NMODES());
	exc.set_max_cloud(get_MAX_CLOUD());
	exc.set_k_3d(get_KVEC());
	exc.set_mol(nmol);

	if (get_NMODES()==2)
	{
		for (i=0; i<=get_MAX_VIB(0)*get_MAX_CLOUD(); i++)		
			for (j=0; j<=get_MAX_VIB(1)*get_MAX_CLOUD(); j++)	
				set_n_vib(cloud_pos, i, j, exc);
	}
	else	
	{	for (i=0; i<=get_MAX_VIB(0)*get_MAX_CLOUD(); i++)		
			set_n_vib_one_mode(cloud_pos, i, exc);
	}

	n_states=BASIS.size();
	BASIS.resize((*latticePtr).get_N_MOL()*n_states);
	for (nmol=1; nmol<(*latticePtr).get_N_MOL(); nmol++)
	{
		for (i=0; i<n_states; i++)
		{
			BASIS[i+nmol*n_states]=BASIS[i];
			BASIS[i+nmol*n_states].set_mol(nmol);
		}
	} 
	cout << "BASIS.size(): " << BASIS.size() << endl;
	return;
}

template <typename FLOAT>
void PhononCloud<FLOAT>::compute_BD_BASIS()
{
	//  Computes sym basis set
	cout << "Computing BD_BASIS ..." << endl;
	if (!(norm_2(get_KVEC())<ZERO_TOL || norm_2(get_KVEC()-get_latticePtr()->get_KLIGHT())<ZERO_TOL))
		nrerror("PhononCloud<FLOAT>::compute_BD_BASIS:: works only for |KVEC|==0!");
	BD_BASIS.resize(2);
	MPSTATE<FLOAT> mps;
	map<string, int> BasisMap;
	map<string, int>::iterator sym_iter, main_iter;
	int i,isym;
	for (i=0; i<2; i++)
		BD_BASIS[i].resize(0);
	for (i=0; i<BASIS.size(); i++)
		BasisMap[BASIS[i].get_string()]=i;
	for ( main_iter=BasisMap.begin(); main_iter!= BasisMap.end(); main_iter++ )
	{
		i=(*main_iter).second;
		sym_iter = BasisMap.find(compute_symmetric_state(1,BASIS[i]).get_string());
		isym=(*sym_iter).second;
		if (isym!=i)
		{
			mps.resize(0);
			mps.push_back(i,1.0/sqrt(2.0));
			mps.push_back(isym,1.0/sqrt(2.0));
			BD_BASIS[0].push_back(mps);
			mps.resize(0);
			mps.push_back(i,1.0/sqrt(2.0));
			mps.push_back(isym,-1.0/sqrt(2.0));
			BD_BASIS[1].push_back(mps);
			BasisMap.erase(sym_iter);
		}
		else
		{
			mps.resize(0);
			mps.push_back(i,1.0);
			BD_BASIS[0].push_back(mps);
		}
	}
	return;
}

template <typename FLOAT>
typename PhononCloud<FLOAT>::ComplexVector PhononCloud<FLOAT>::get_dipole(ME<FLOAT>& exc)
{
	ComplexVector	dipole(3);
	if (exc.get_np()==1 && (norm_2(get_KVEC())<ZERO_TOL || norm_2(get_KVEC()-get_latticePtr()->get_KLIGHT())<ZERO_TOL))
	{
			dipole=(*latticePtr).get_MOL_DIPOLE(exc.get_mol());
			for (int mode=0; mode<get_NMODES(); mode++)
				dipole*=FC_overlap(mode, exc.get_cloud(mode, 0), 0);
	}
	else
	{
		dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
	}
	return(dipole);
}

template <typename FLOAT>
typename PhononCloud<FLOAT>::ComplexVector PhononCloud<FLOAT>::get_dipole(MPSTATE<FLOAT>& mps)
{
	ComplexVector	dipole(3);
	dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
	for (int i=0; i<mps.size(); i++)
		dipole+=mps.get_coeff(i)*get_dipole(BASIS[mps.get_state(i)]);
	return(dipole);
}
/*
template <typename FLOAT>
vector<vector<FLOAT> > PhononCloud<FLOAT>::compute_emission(QuantumState<FLOAT>&	emitting_state, Vector& em_dir)
{
	if (FLAG_DIMER==1)
			compute_DIMER_emission(emitting_state, em_dir);
	else
	return(compute_crystal_emission(emitting_state, em_dir));
}*/

// It works for (NP==1 and NMODES<=2) and for (NP==2 and NMODES==1)
template <typename FLOAT>
vector<vector<FLOAT> > PhononCloud<FLOAT>::compute_emission(QuantumState<FLOAT>&	emitting_state, Vector& em_dir)
{
//	cout << "Computing crystal emission in vibronic approximation ... " << endl;
	if ( get_NP()>2 )
		nrerror("compute_crystal_emission: it works only with one- or two-particle states! Try setting NP=1 or NP=2");
	if (get_NMODES()>2)
		nrerror("compute_crystal_emission: it works only with 1 or 2 vibrational modes! Try setting NMODES<=2");
	
	int i, j, m, mode, state, dwtot, hitot, mol, vib_pos, v1;
	vector<int>	FinVibs(2);
	Complex I(0.0,1.0);
	int dwmax, himax;	

	ComplexVector dipole(3), temp_dipole(3);
	vector<vector<FLOAT> > emission;

// The size of "emission" should be always as below. If we are not computing some of the peaks, we leave it filled with 0
	dwmax=get_MAX_VIB(0);
	if (get_NMODES()>1)
		himax=get_MAX_VIB(1);
	else
		himax=0; 
	
	emission.resize(dwmax+1);
	for (i=0; i<=dwmax; i++)
		emission[i].resize(himax+1);
	
	for (i=0; i<emission.size(); i++) for (j=0; j<emission[i].size(); j++)
		emission[i][j]=0.0;
	
///////////////////////////////////////////////////////////////
// This is the routine for 2 modes and only one particle states	
	if (get_NMODES()==2)
	{
		if ( get_NP()>1 )
			nrerror("compute_crystal_emission: if there are two modes it works only with one-particle states! Try setting NP=1!");
		
		for (dwtot=0; dwtot<=dwmax; dwtot++)	for (hitot=0; hitot<=himax; hitot++)
		{
//			cout << dwtot << hitot << endl;
			FinVibs[0]=dwtot;
			FinVibs[1]=hitot;
			
//---------- 0-0 emission
			if (dwtot+hitot==0)
			{
				dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
				if (norm_2(get_KVEC())<ZERO_TOL || norm_2(get_KVEC()-get_latticePtr()->get_KLIGHT())<ZERO_TOL)
				{
					for (state=0; state<emitting_state.size(); state++)
					{
						temp_dipole[0]=0.0; temp_dipole[1]=0.0;temp_dipole[2]=0.0;
						for (i=0; i<3; i++)
							temp_dipole(i)=(*latticePtr).get_MOL_DIPOLE(BASIS[state].get_mol())(i)*emitting_state.get_coeff(state)
							         *exp(-I*inner_prod(get_KVEC(),(*latticePtr).get_RHO(BASIS[state].get_mol())));
						for (mode=0; mode<get_NMODES(); mode++)
							temp_dipole*=FC_overlap(mode, BASIS[state].get_cloud(mode, 0), FinVibs[mode]);
						dipole+=temp_dipole;
					}
				}
				emission[dwtot][hitot]=pow(abs(inner_prod(dipole,em_dir)),2);
			}
			else
//---------- 0-n emission
			{
				emission[dwtot][hitot]=0.0;
				for (mol=0; mol<(*latticePtr).get_N_MOL(); mol++)
				{
					dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
					for (state=0; state<emitting_state.size(); state++)
					{
						if (BASIS[state].get_mol()==mol)
						{
							temp_dipole[0]=0.0; temp_dipole[1]=0.0;temp_dipole[2]=0.0;
							for (i=0; i<3; i++)
								temp_dipole(i)=(*latticePtr).get_MOL_DIPOLE(mol)(i)*emitting_state.get_coeff(state)
								         *exp(-I*inner_prod(get_KVEC(),(*latticePtr).get_RHO(mol)));
							for (mode=0; mode<get_NMODES(); mode++)
								temp_dipole*=FC_overlap(mode, BASIS[state].get_cloud(mode, 0), FinVibs[mode]);
							dipole+=temp_dipole;
						}
					}
					emission[dwtot][hitot]+=pow(abs(inner_prod(dipole,em_dir)),2);
				}
				
			}
		}
	}
	
///////////////////////////////////////////////////////////////
// This is the routine for 1 mode and up to two-particle states	
	else  
	{	
		hitot=0;
//--------- 0-0 emission it only contains one-particle states
		dwtot=0;			
		emission[dwtot][hitot]=compute_emission_I0(emitting_state, em_dir);
//--------- 0-n emission (n>0)
		for (dwtot=1; dwtot<=dwmax; dwtot++)
		{
			emission[dwtot][hitot]=compute_emission_I1(emitting_state, em_dir, dwtot);
				
			for (v1=1; v1<dwtot; v1++)		
				emission[dwtot][hitot]+=compute_emission_I2(emitting_state, em_dir, v1, dwtot-v1);									   	
		}
	}
	return(emission);
}

template <typename FLOAT>
FLOAT PhononCloud<FLOAT>::compute_emission_I0(QuantumState<FLOAT>&	emitting_state, Vector& em_dir)
{
	int i, state, mode;
	FLOAT I0=0.0;
	Complex I(0.0,1.0);
	ComplexVector dipole(3), temp_dipole(3);
		
	if (norm_2(get_KVEC())<ZERO_TOL || norm_2(get_KVEC()-get_latticePtr()->get_KLIGHT())<ZERO_TOL)
	{
		dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
		for (state=0; state<emitting_state.size(); state++)
		{
			temp_dipole[0]=0.0; temp_dipole[1]=0.0;temp_dipole[2]=0.0;
			for (i=0; i<3; i++)
				temp_dipole(i)=(*latticePtr).get_MOL_DIPOLE(BASIS[state].get_mol())(i)*emitting_state.get_coeff(state)
				*exp(-I*inner_prod(get_KVEC(),(*latticePtr).get_RHO(BASIS[state].get_mol())));
			for (mode=0; mode<get_NMODES(); mode++)
				temp_dipole*=FC_overlap(mode, BASIS[state].get_cloud(mode, 0), 0);
			dipole+=temp_dipole;
		}
		I0=pow(abs(inner_prod(dipole,em_dir)),2);
	}
	return(I0);
}

template <typename FLOAT>
FLOAT PhononCloud<FLOAT>::compute_emission_I1(QuantumState<FLOAT>&	emitting_state, Vector& em_dir, int vtot)
{
	int i, mol, state, mode, m, vib_pos;
	FLOAT I1=0.0;
	Complex I(0.0,1.0);
	ComplexVector dipole(3), temp_dipole(3);
	
	for (mol=0; mol<(*latticePtr).get_N_MOL(); mol++)
	{
		dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
		for (state=0; state<emitting_state.size(); state++)
		{
			// This is the NP==1 term
			if (BASIS[state].get_np()==1 && BASIS[state].get_mol()==mol)
			{
				temp_dipole[0]=0.0; temp_dipole[1]=0.0;temp_dipole[2]=0.0;
				for (i=0; i<3; i++)
					temp_dipole(i)=(*latticePtr).get_MOL_DIPOLE(mol)(i)*emitting_state.get_coeff(state)
					*exp(-I*inner_prod(get_KVEC(),(*latticePtr).get_RHO(mol)));
				for (mode=0; mode<get_NMODES(); mode++)
					temp_dipole*=FC_overlap(mode, BASIS[state].get_cloud(mode, 0), vtot);
				dipole+=temp_dipole;
			}
			
			// This is the NP==2 term
			if ( BASIS[state].get_np()==2 )
			{
				mode=0;
				vib_pos=0;
				temp_dipole[0]=0.0; temp_dipole[1]=0.0;temp_dipole[2]=0.0;
				for (m=1; m<BASIS[state].get_max_cloud(); m++)   // Finds the position of the (non electronically excited) molecular vibration
					if (BASIS[state].get_cloud(mode,m)>0)
						vib_pos=m;	
				if (vib_pos==0) nrerror("compute_crystal_emission_I1: could not find vib_pos in a two-particle state");
				if (BASIS[state].get_cloud(mode,vib_pos)==vtot && get_CLOUD_NN_alpha(BASIS[state].get_mol(),vib_pos)==mol)
				{
					for (i=0; i<3; i++)
						temp_dipole(i)=(*latticePtr).get_MOL_DIPOLE(BASIS[state].get_mol())(i)*emitting_state.get_coeff(state)
						*exp(-I*inner_prod(get_KVEC(),(*latticePtr).get_RHO(BASIS[state].get_mol())))
						*exp(-I*inner_prod(get_KVEC(),get_CLOUD_NN_n(BASIS[state].get_mol(),vib_pos)));
					
					temp_dipole*=FC_overlap(mode, BASIS[state].get_cloud(mode, 0), 0);
					dipole+=temp_dipole;
				}							
			}
			
		}
		I1+=pow(abs(inner_prod(dipole,em_dir)),2);
	}
	return(I1);
}

template <typename FLOAT>
FLOAT PhononCloud<FLOAT>::compute_emission_I2(QuantumState<FLOAT>&	emitting_state, Vector& em_dir, int v1, int v2)
{
	int i, mol1,mol2, nn, state, mode, m, vib_pos, nninv;
	FLOAT I2=0.0;
	Complex I(0.0,1.0);
	ComplexVector dipole(3), temp_dipole(3);
	Vector r;
	
	for (mol1=0; mol1<(*latticePtr).get_N_MOL(); mol1++) for (nn=0; nn<get_MAX_CLOUD(); nn++) 
	{
		r=get_CLOUD_NN_n(mol1,nn);
		mol2=get_CLOUD_NN_alpha(mol1,nn);
		
		vib_pos=0;
		for (vib_pos=0; vib_pos<get_MAX_CLOUD(); vib_pos++)
			if (get_CLOUD_NN_alpha(mol2,vib_pos)==mol1 && get_CLOUD_NN_n(mol2,vib_pos)==-r)
				nninv=vib_pos;
		if (vib_pos==0) nrerror("compute_crystal_emission_I2: could not find the inverse nearest neighbours");
		
		dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
		for (state=0; state<emitting_state.size(); state++)
		{			
			// This is the NP==2 term
			if ( BASIS[state].get_np()==2 )
			{
				mode=0;
				vib_pos=0;
				temp_dipole[0]=0.0; temp_dipole[1]=0.0;temp_dipole[2]=0.0;
				for (m=1; m<BASIS[state].get_max_cloud(); m++)   // Finds the position of the (non electronically excited) molecular vibration
					if (BASIS[state].get_cloud(mode,m)>0)
						vib_pos=m;	
				if (vib_pos==0) nrerror("compute_crystal_emission_I2: could not find vib_pos in a two-particle state");
				
				if (BASIS[state].get_mol()==mol1 && vib_pos==nn && BASIS[state].get_cloud(mode,vib_pos)==v2 )
				{
					for (i=0; i<3; i++)
						temp_dipole(i)=(*latticePtr).get_MOL_DIPOLE(mol1)(i)*emitting_state.get_coeff(state)
						*exp(-I*inner_prod(get_KVEC(),(*latticePtr).get_RHO(mol1)));
					
					temp_dipole*=FC_overlap(mode, BASIS[state].get_cloud(mode, 0), v1);
					dipole+=temp_dipole;
				}	
				else
				{

					if (BASIS[state].get_mol()==mol2 && vib_pos==nninv && BASIS[state].get_cloud(mode,vib_pos)==v1)
					{
						for (i=0; i<3; i++)
							temp_dipole(i)=(*latticePtr).get_MOL_DIPOLE(mol2)(i)*emitting_state.get_coeff(state)
							*exp(-I*inner_prod(get_KVEC(),(*latticePtr).get_RHO(mol2)))
							*exp(-I*inner_prod(get_KVEC(),r));
					
						temp_dipole*=FC_overlap(mode, BASIS[state].get_cloud(mode, 0), v2);
						dipole+=temp_dipole;
					}
				}
			}
			
		}
		I2+=pow(abs(inner_prod(dipole,em_dir)),2);
	}
	return(I2);
}

// nue is the number of quanta on excited molecules
template <typename FLOAT>
FLOAT PhononCloud<FLOAT>::FC_overlap(int mode, int nue, int nug)
{
	FLOAT result=0.0;
	if (mode==0)
	{
		if (FLAG_DW==1)
			result=DW_FCF_table[nug][nue];
		else
			result=S_factor(nue, nug, lambda_0);
	}
	if (mode==1)
	{
		result=S_factor(nue, nug, lambda_1);
	}
	return(result);
}

template <typename FLOAT>
FLOAT PhononCloud<FLOAT>::Vibrational_energies(int mode, int flag_exc, int n)
{
	FLOAT result=0.0;
	if (mode==0)
	{
		if (FLAG_DW==1)
		{
			if (flag_exc==1)	
				result=(n+0.5)*DW_hwe-DW_levels[0];
			else					
				result=DW_levels[n]-DW_levels[0]; 
		}
		else
			result=(n)*vib_En_0;
	}
	if (mode==1)
		result=(n)*vib_En_1;
	return(result);
}



template <typename FLOAT>
typename PhononCloud<FLOAT>::Complex PhononCloud<FLOAT>::Hint(ME<FLOAT> &exc1, ME<FLOAT> &exc2, int flag_macro)
{
	int mode;
	int i,j, n,m,l;
	int nmax, mmax, lmax; 

	Vector dr,r;
	Complex I(0.0,1.0);
	Complex temp_sum;
	Complex sum=0.0;
	
/*	if (FLAG_DIMER==1)
	{
		nmax=0;
		mmax=0;
		lmax=0;
	}
	else
	{*/
		nmax=int(2*get_CLOUD_RADIUS()/norm_2((*latticePtr).get_CELL(0)))+1; 
		mmax=int(2*get_CLOUD_RADIUS()/norm_2((*latticePtr).get_CELL(1)))+1; 
		if (FLAG_2D==1)    // If the mode is "Single layer" then there is no summation in the z direction
			lmax=0;
		else
			lmax=int(2*get_CLOUD_RADIUS()/norm_2((*latticePtr).get_CELL(2)))+1;
//	}

	if (exc1.get_k()==exc2.get_k())
	{
				
// Vibronic states (there are no vibrations on electronically unexcited molecules)
		if (exc1.get_np()==1 && exc2.get_np()==1)
		{
			sum=(*latticePtr).compute_Ltilde(get_KVEC(), exc1.get_mol(), exc2.get_mol(), flag_macro);  // I use LTILDE because it is an infinite sum over dr 
			for (i=0; i<get_NMODES(); i++)
				if (get_MAX_VIB(i)>0)  ////////////////////  Check if this can be removed by putting lambda=0 for purely excitonic states
					sum*=FC_overlap(i, exc1.get_cloud(i,0),0)*FC_overlap(i, exc2.get_cloud(i, 0),0);
			if (exc1==exc2) 
			{
				sum+=elec_En;
				for (i=0; i<get_NMODES(); i++)
					sum+=Vibrational_energies(i, 1, exc1.get_cloud(i,0));
			}
		}
		else
// General case
		{	
			// Diagonal terms
			if (exc1==exc2) 
			{
				sum+=elec_En;
				for (mode=0; mode<get_NMODES(); mode++)
				{
					sum+=Vibrational_energies(mode, 1, exc1.get_cloud(mode,0));
					for (i=1;i<get_MAX_CLOUD();i++)
						sum+=Vibrational_energies(mode, 0, exc1.get_cloud(mode,i));
				}
			}
			else
			{
			//	Off diagonal terms
				for (n=-nmax;n<=nmax;n++)  for (m=-mmax;m<=mmax;m++)  for (l=-lmax;l<=lmax;l++)
				{
						dr=(*latticePtr).get_RHO(exc1.get_mol())-
							(
							(*latticePtr).get_RHO(exc2.get_mol())+(n*(*latticePtr).get_CELL(0)+m*(*latticePtr).get_CELL(1)+l*(*latticePtr).get_CELL(2))
						    );
		
						// I'm excluding the self-interaction and states too far
						if (norm_2(dr)>ZERO_TOL && norm_2(dr)<2.0*get_CLOUD_RADIUS()+ZERO_TOL)   
						{
							temp_sum=2.0;  // 2.0 is just a number used to check if temp_sum has been set to 0
						// This two loops check if clouds are orthogonal and then set temp_sum=0.0
							for (mode=0; mode<get_NMODES(); mode++)
							{
								for (j=1; j<get_MAX_CLOUD(); j++)	
								{		
									if (exc1.get_cloud(mode,j)!=exc2.get_vib(mode, latticePtr, CLOUD_NN[exc1.get_mol()][j].n+dr) && 
										norm_2((*latticePtr).get_NN_n(exc1.get_mol(),j)+dr) > 1e-8)
									{
										temp_sum=0.0;
										break;
									}
									if (exc2.get_cloud(mode,j)!=exc1.get_vib(mode, latticePtr, CLOUD_NN[exc2.get_mol()][j].n-dr) 
										&& norm_2(CLOUD_NN[exc2.get_mol()][j].n-dr) > 1e-8 )
									{
										temp_sum=0.0;
										break;
									}
								}
								if (abs(temp_sum)<ZERO_TOL)
									break;
							}
						// Here we actually compute off-diagonal elements
							if (temp_sum.real() > 1.0)   // if temp_sum has not been set to 0
							{
								temp_sum=exp(I*inner_prod(get_KVEC(),dr))*(*latticePtr).interaction(exc1.get_mol(), exc2.get_mol(), dr);
								for (mode=0; mode<get_NMODES(); mode++)
									if (get_MAX_VIB(mode)>0)
										temp_sum*=FC_overlap(mode, exc1.get_cloud(mode, 0),exc2.get_vib(mode, latticePtr,  dr))
										         *FC_overlap(mode, exc2.get_cloud(mode, 0),exc1.get_vib(mode, latticePtr, -dr));
							}					
							sum+=temp_sum;			
					}
				}				// end of the for loop for the general case
			}					// end else (exc1==exc2)			
		}						// end else (exc1.get_np()==1 && exc2.get_np()==1)			
					
/*		if (abs(imag(sum))/abs(sum) > 1e-8 ) 
		{	
			cout << sum << endl;
			nrerror("Hint: returning an imaginary value!");
		}
		return(sum.real());*/
		return(sum);
	}
	else
		return(0.0);
}

template <typename FLOAT>
typename PhononCloud<FLOAT>::Complex PhononCloud<FLOAT>::Hint(MPSTATE<FLOAT>  &mps1, MPSTATE<FLOAT>  &mps2, int flag_macro)
{
	Complex energy = 0.0;
	for (int i = 0; i< mps1.size(); i++)
		for (int j = 0; j< mps2.size(); j++)
			energy += mps1.get_coeff(i)*mps2.get_coeff(j)*Hint(BASIS[mps1.get_state(i)], BASIS[mps2.get_state(j)],flag_macro);
	return(energy);
}


template <typename FLOAT>
typename PhononCloud<FLOAT>::Complex PhononCloud<FLOAT>::Hint(int i, int j, int DC, int flag_macro)
{
	if (DC<0)
	{
/*		if (FLAG_DIMER==1)
			return(Hint_dimer_DW(BASIS[i], BASIS[j]));
		else*/
			return(Hint(BASIS[i], BASIS[j],flag_macro));
	}
	else
	{
		return(Hint(BD_BASIS[DC][i], BD_BASIS[DC][j],flag_macro));
	}
}


template <typename FLOAT>
bool PhononCloud<FLOAT>::check_BD_decomposition(int dc1, int dc2, int flag_macro, FLOAT tol)
{
	int i,j;
	cout << "Checking BD decomposition: " << endl;
	for (i=0; i< BD_BASIS[dc1].size(); i++)
		for (j=0; j< BD_BASIS[dc2].size(); j++)
		{
	//		cout << i << "\t" << j << endl;
			if (abs(Hint(BD_BASIS[dc1][i], BD_BASIS[dc2][j],flag_macro)) > tol)
			{
				cout << i << "\t" << j << "\t" << Hint(BD_BASIS[dc1][i], BD_BASIS[dc2][j],flag_macro) << endl;
				return(false);
			}
		}
	return(true);
}

/*
 // It works only for 2 vibrational modes
 template <typename FLOAT>
 void PhononCloud<FLOAT>::compute_DIMER_emission(QuantumState<FLOAT>&	emitting_state, Vector& em_dir)
 {
 if (get_NMODES()!=2)
 nrerror("PhononCloud::compute_DIMER_emission: emission can be computed only for NMODES=2!");
 int i, j, mode, state, dwtot, hitot;
 vector<vector<int> >	FinVibs(2);
 for (i=0; i<2; i++)
 FinVibs[i].resize(2);
 ComplexVector dipole(3), temp_dipole(3);
 FLOAT sum;
 vector<vector<FLOAT> > emission;
 emission.resize(get_MAX_VIB(0)+1);
 for (i=0; i<get_MAX_VIB(0)+1; i++)
 emission[i].resize(get_MAX_VIB(1)+1);
 for (i=0; i<emission.size(); i++) for (j=0; j<emission[i].size(); j++)
 emission[i][j]=0.0;
 
 for (dwtot=0; dwtot<=get_MAX_VIB(0); dwtot++)		for (hitot=0; hitot<=get_MAX_VIB(1); hitot++)
 {
 sum=0.0;
 for (FinVibs[0][0]=0; FinVibs[0][0]<=dwtot; FinVibs[0][0]++)	for (FinVibs[1][0]=0; FinVibs[1][0]<=hitot; FinVibs[1][0]++)
 {
 FinVibs[0][1]=dwtot-FinVibs[0][0];
 FinVibs[1][1]=hitot-FinVibs[1][0];
 dipole[0]=0.0; dipole[1]=0.0; dipole[2]=0.0;
 for (state=0; state<emitting_state.size(); state++)
 {
 temp_dipole[0]=0.0; temp_dipole[1]=0.0;temp_dipole[2]=0.0;
 if (BASIS[state].get_cloud(0, 1)==FinVibs[0][(BASIS[state].get_mol()+1)%2]  &&
 BASIS[state].get_cloud(1, 1)==FinVibs[1][(BASIS[state].get_mol()+1)%2] )
 {
 //						cout << temp_dipole << endl;
 for (i=0; i<3; i++)
 temp_dipole(i)=(*latticePtr).get_MOL_DIPOLE(BASIS[state].get_mol())(i)*emitting_state.get_coeff(state);
 for (mode=0; mode<2; mode++)
 temp_dipole*=FC_overlap(mode, BASIS[state].get_cloud(mode, 0), FinVibs[mode][BASIS[state].get_mol()]);
 //						cout << temp_dipole << endl;
 }
 dipole+=temp_dipole;
 }
 sum+=pow(abs(inner_prod(dipole,em_dir)),2);
 }
 emission[dwtot][hitot]=sum;
 }
 //	}
 emitting_state.set_emission(emission);
 return;
 }*/
/*
 template <typename FLOAT>
 FLOAT PhononCloud<FLOAT>::Hint_dimer_DW(ME<FLOAT> &exc1, ME<FLOAT> &exc2)
 {
 int mode;
 int i,j, n,m,l;
 FLOAT	sum=0.0;
 
 // Interaction between two Single Particle states
 if (exc1.get_np()+exc2.get_np()==2)
 {
 if (exc1==exc2)	// Diagonal term for single particle states
 sum=get_LTILDE(exc1.get_mol(),exc2.get_mol()).real() + elec_En+(exc1.get_cloud(0,0)+0.5)*DW_hwe-DW_levels[0]
 +(exc1.get_cloud(1,0))*vib_En_1;
 else // Off-Diagonal term for single particle states
 {
 sum=get_LTILDE(exc1.get_mol(),exc2.get_mol()).real();   
 for (i=0; i<get_NMODES(); i++)
 sum*=FC_overlap(i, exc1.get_cloud(i,0),0)*FC_overlap(i, exc2.get_cloud(i, 0),0);
 }
 }
 else
 // General case
 {	
 // Diagonal term for non single particle states
 if (exc1==exc2) 
 {
 sum=elec_En;
 for (mode=0; mode<get_NMODES(); mode++)
 {
 sum+=Vibrational_energies(mode, 1, exc1.get_cloud(mode,0));
 for (i=1;i<get_MAX_CLOUD();i++)
 sum+=Vibrational_energies(mode, 0, exc1.get_cloud(mode,i));
 }
 }
 else
 {
 //	Off diagonal terms
 if (exc1.get_mol()==exc2.get_mol())  
 return(0.0);
 else
 {			
 sum=get_LTILDE(exc1.get_mol(),exc2.get_mol()).real();
 for (mode=0; mode<get_NMODES(); mode++)
 sum*=FC_overlap(mode, exc1.get_cloud(mode, 0), exc2.get_cloud(mode, 1))
 *FC_overlap(mode, exc2.get_cloud(mode, 0), exc1.get_cloud(mode, 1));
 }				
 }					// end else (exc1==exc2)			
 }						// end else (exc1.get_np()+exc2.get_np()==2)							
 return(sum);
 }*/


// Explicit instantiation

template class PhononCloud<float>;
template class PhononCloud<double>;
