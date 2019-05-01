/* Written by Leonardo Silvestri 2007-2013 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <N_Functions.h>
#include <OPutils.h>
#include <ME.h>
#include <lattice.h>

namespace ublas = boost::numeric::ublas;
using namespace std;

#ifndef ZERO_TOL
#define ZERO_TOL 1e-10
#endif

///////////////////////////////////////////////////////////////////////////
/////////          Molecular Exciton with phonon cloud class
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Long exciton constructor with multiple modes

template<typename FLOAT>	
ME<FLOAT>::ME(FLOAT _kx, FLOAT _ky, FLOAT _kz, int _mol, std::vector<std::vector<int> >& _cloud) 
{
	k.resize(3);
	k[0]=_kx;k[1]=_ky;k[2]=_kz; 
	mol=_mol; 
	int i,j;
	nmodes=_cloud.size();
	max_cloud=0;
	for (i=0; i<nmodes; i++)
		if (cloud[i].size()>max_cloud)
			max_cloud=cloud[i].size();
	cloud.resize(nmodes);
	for (i=0; i<nmodes; i++)
	{
		cloud[i].resize(max_cloud);
		for (j=0; j<max_cloud; j++) 
			if (j<_cloud[i].size() && _cloud[i][j]>=0) 
				cloud[i][j]=_cloud[i][j];
			else
				cloud[i][j]=0;
	}
	set_np();
}; 

///////////////////////////////////////////////////////////////////////////
// Functions to set parameters
	
// Computes the number of particles for the exciton and modifies np accordingly
template<typename FLOAT>	
void ME<FLOAT>::set_np()
	{
		int i,j;
		np=1;
		for (j=1; j<max_cloud; j++) 
			for (i=0; i<nmodes; i++)	
				if (cloud[i][j]>0) 
				{
					np++;
					break;
				}
		return;
	};
	
template<typename FLOAT>	
void ME<FLOAT>::set_k_3d(Vector param) 
{
	if (param.size()!=3) 
		nrerror("ME<FLOAT>::set_k_3d: wrong size, cannot set k!");
	k=param;
	return;
}

template<typename FLOAT>	
void ME<FLOAT>::set_k_3d(FLOAT _kx, FLOAT _ky, FLOAT _kz) 
{
	k.resize(3);
	k[0]=_kx;k[1]=_ky;k[2]=_kz;
	return;
};

template<typename FLOAT>	
void ME<FLOAT>::set_cloud(int _mode, int _npos, int _nvib) 
{
	if (_mode>=0 && _mode<nmodes)
		if (_npos<max_cloud && _npos>=0 && _nvib>=0) 
		{
			cloud[_mode][_npos]=_nvib;
			set_np();
		}
	return;
};

template<typename FLOAT>	
void ME<FLOAT>::set_nmodes(int _nmodes) 
{
	if (_nmodes>=0)
	{
		nmodes=_nmodes;
		cloud.resize(nmodes);
	}
	return;
};

template<typename FLOAT>	
void ME<FLOAT>::set_max_cloud(int _max_cloud) 
{
	if (_max_cloud>=0)
	{
		max_cloud=_max_cloud;
		for (int i=0; i<nmodes; i++)
			cloud[i].resize(max_cloud);
	}
	return;
};

///////////////////////////////////////////////////////////////////////////
// Functions to get parameters
	
template<typename FLOAT>	
int ME<FLOAT>::get_cloud(int mode, int i)  const
{
		if (mode>=0 && mode<get_nmodes())
			if (i>=0 && i<get_max_cloud())
				return(cloud[mode][i]);
			else
				return(0);
		nrerror("ME<FLOAT>::get_cloud: 'mode' out of bounds!");
		return(0);
};

	
// Returns the number of vibrational quanta at relative position r
template<typename FLOAT>	
int ME<FLOAT>::get_vib(int mode, lattice<FLOAT>*	latticePtr, Vector r)
	{	
		for (int i=0; i<get_max_cloud(); i++)
			if (norm_2((*latticePtr).get_NN_n(get_mol(),i)-r)<ZERO_TOL)
				return(get_cloud(mode,i));
		return(0);
	}

template<typename FLOAT>	
string ME<FLOAT>::get_string() const
	{
		ostringstream ss;
		ss << get_mol();
		for (int i=0; i<get_cloud_size(); i++)
			for (int j=0; j<get_cloud_size(i); j++)
				ss << get_cloud(i,j);
		return(ss.str());
	}

template<typename FLOAT>	
bool ME<FLOAT>::operator==(const ME<FLOAT>& exc) const
{
	if (get_nmodes()!= exc.get_nmodes() || get_max_cloud()!= exc.get_max_cloud())
		return(false);
	for (int i=0; i<get_nmodes(); i++)		for (int j=0; j<get_max_cloud(); j++)
		if (get_cloud(i,j)!=exc.get_cloud(i,j)) return(false);
	if (get_k()==exc.get_k() && get_mol()==exc.get_mol())
		return(true);
	return(false);
} 

// Screen output of a molecular exciton 
template <typename FLOAT>
ostream &operator<<(ostream& stream, const ME<FLOAT>& exc)   
{
	int i,j;
	stream << "|" << exc.get_k() << "," << exc.get_mol() ;
	for (i=0; i<exc.get_nmodes(); i++) 
	{
		stream <<  ";" ;
		for (j=0; j<exc.get_max_cloud(); j++)
			stream << exc.get_cloud(i,j) ;	
	}
	stream << "> (np=" << exc.get_np() << ")" ;
	return stream;
};

// Explicit instantiation
template class ME<float>;
template class ME<double>;
template ostream &operator<<(ostream& stream, const ME<float>& exc);
template ostream &operator<<(ostream& stream, const ME<double>& exc); 
