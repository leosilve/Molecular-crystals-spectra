/* Written by Leonardo Silvestri 2007-2013 */

#include <boost/numeric/ublas/vector.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <QuantumState.h>
#include <OPutils.h>

using namespace std;


#ifndef Pi
#define Pi 3.141592654
#endif 
#ifndef ZERO_TOL
#define ZERO_TOL 1e-10
#endif

	
template<class FLOAT>
void QuantumState<FLOAT>::set_coeff(int _pos, complex<FLOAT> _coeff) 
{ 
	if (_pos>=0 && _pos < coeff.size() ) 
		coeff[_pos]=_coeff;
	return;
};

template<class FLOAT>
bool QuantumState<FLOAT>::normalize() 
{
		int i;
		FLOAT sum=0.0;
		for (i=0;i<coeff.size();i++)
			sum+=pow(abs(coeff[i]),2);
		if (sum<ZERO_TOL)
			return(false);
		else
			for (i=0;i<size();i++)
			{
				if (abs(coeff[i])>ZERO_TOL)
					coeff[i]/=sqrt(sum);
				else
					coeff[i]=0.0;
			}
		return(true);
}

template<class FLOAT>
const QuantumState<FLOAT> QuantumState<FLOAT>::operator+ (const QuantumState<FLOAT>& param)  const
{
		if (size()!=param.size())
			nrerror("QuantumState<FLOAT> operator+ : QuantumStates have different sizes!");
		QuantumState<FLOAT> temp(size());
		for (int i=0; i<size(); i++)
			temp.set_coeff(i,get_coeff(i)+param.get_coeff(i));
		return (temp);
};


template<class FLOAT>
const QuantumState<FLOAT> QuantumState<FLOAT>::operator* (const complex<FLOAT> param) const 
{
	QuantumState<FLOAT> temp(size());
	for (int i=0; i<size(); i++)
		temp.set_coeff(i,get_coeff(i)*param);
	return (temp);
};

template <typename FLOAT>
ostream &operator<<(ostream& stream, const QuantumState<FLOAT>& qs)   // Screen output of a QuantumState
{
	if (qs.size()==0)
		stream << "| >" ;
	else
	{
		stream << "| " ;
		for (int i=0; i<qs.size()-1; i++)
			stream  << std::fixed << std::setprecision(2) << qs.get_coeff(i) << "(" << i << ")" << "+" ;
		stream  << std::fixed << std::setprecision(2) << qs.get_coeff(qs.size()-1)  << "(" << qs.size()-1 << ") >" ;
		stream  << " en: " << std::fixed << std::setprecision(4) << qs.get_energy(); 
	}
	return stream;
};

// Explicit instantiation

template class QuantumState<float>;
template class QuantumState<double>;

template ostream &operator<<(ostream& stream, const QuantumState<float>& qs);   // Screen output of a QuantumState
template ostream &operator<<(ostream& stream, const QuantumState<double>& qs);