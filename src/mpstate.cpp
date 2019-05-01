/* Written by Leonardo Silvestri 2007-2013 */

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <mpstate.h>

#ifndef ZERO_TOL
#define ZERO_TOL 1e-10
#endif

///////////////////////////////////////////////////////////////////////////
/////////          Block diagonal Basis set

using namespace std;

template <class FLOAT>
const MPSTATE<FLOAT> MPSTATE<FLOAT>::operator+ (const MPSTATE& param) const {
	MPSTATE temp(size());
	int i, j;
	for (i=0; i<size(); i++)
	{
		temp.set_state(i,get_state(i));
		temp.set_coeff(i,get_coeff(i));
	}
	for (i=0; i<param.size(); i++)
	{
		for (j=0; j<size(); j++)
		{
			if (temp.get_state(j)==param.get_state(i))
			{
				temp.set_coeff(j,get_coeff(j)+param.get_coeff(i));
				break;
			}
			
		}
		if (j==size())
			temp.push_back(param.get_state(i),param.get_coeff(i));
	}
	return (temp);
};

template <class FLOAT>
const MPSTATE<FLOAT> MPSTATE<FLOAT>::operator* (const FLOAT param) const {
	
	MPSTATE temp(size());
	int i;
	for (i=0; i<size(); i++)
	{
		temp.set_state(i,get_state(i));
		temp.set_coeff(i,get_coeff(i)*param);
	}
	return (temp);
};

template <class FLOAT>
const MPSTATE<FLOAT> MPSTATE<FLOAT>::operator* (const complex<FLOAT> param) const {
	MPSTATE temp(size());
	int i;
	for (i=0; i<size(); i++)
	{
		temp.set_state(i,get_state(i));
		temp.set_coeff(i,get_coeff(i)*param);
	}
	return (temp);
};

template <class FLOAT>
bool MPSTATE<FLOAT>::normalize() 
{
	int i;
	FLOAT sum=0.0;
	for (i=0;i<size();i++)
		sum+=pow(abs(get_coeff(i)),2);
	if (sum< ZERO_TOL)
		return(false);
	for (i=0;i<size();i++)
		if (abs(get_coeff(i))/sum < ZERO_TOL)
		{
			sum-=pow(abs(get_coeff(i)),2);
			erase(i);
			i--;
		}
	for (i=0;i<size();i++)
		set_coeff(i, get_coeff(i)/sqrt(sum));
	return(true);
}
	
template <typename FLOAT>
ostream &operator<<(ostream& stream, const MPSTATE<FLOAT>& mps)   // Screen output of a MPSTATE 
{
	if (mps.size()==0)
		stream << "| >" ;
	else
	{
		stream << "| " ;
		for (int i=0; i<mps.size()-1; i++)
			stream  << std::fixed << std::setprecision(2) << mps.get_coeff(i) << "(" << mps.get_state(i) << ")" << "+" ;
		stream  << std::fixed << std::setprecision(2) << mps.get_coeff(mps.size()-1)  << "(" << mps.get_state(mps.size()-1) << ") >" ;
		stream  << " en: " << std::fixed << std::setprecision(4) << mps.get_energy(); 
	}
	return stream;
};

// Explicit instantiation
template class MPSTATE<float>;
template class MPSTATE<double>; 
template ostream &operator<<(ostream& stream, const MPSTATE<float>& mps);
template ostream &operator<<(ostream& stream, const MPSTATE<double>& mps);
