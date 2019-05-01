/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _QUANTUMSTATE_H
#define _QUANTUMSTATE_H

#include <boost/numeric/ublas/vector.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <complex>
#include <vector>


///////////////////////////////////////////////////////////////////////////
/////////          Quantum state expressed in BASIS set

template <class FLOAT>
class QuantumState {
private:
	std::vector<std::complex<FLOAT> > coeff;
	FLOAT energy;
public:	
	QuantumState() {
		coeff.resize(0);
		energy=0.0;
	};

	QuantumState(int _size) {
		coeff.resize(_size);
		for (int i=0; i<_size; i++)
			coeff[i]=0.0;
		energy=0.0;
	}
	
	int size() const { return(coeff.size()); };
	std::complex<FLOAT> get_coeff(int _i) const { return(coeff[_i]); };
	FLOAT get_energy() const { return(energy); };
		
	void resize(int _size) { coeff.resize(_size); };
	void set_coeff(int _pos, std::complex<FLOAT> _coeff);
	void set_energy(FLOAT _energy) { energy=_energy; };

	bool normalize();

	const QuantumState<FLOAT> operator+ (const QuantumState<FLOAT>& param) const;
	const QuantumState<FLOAT> operator* (const std::complex<FLOAT> param) const;

}; // end of class QuantumState


inline bool operator< ( const QuantumState<float>& qs1, const QuantumState<float>& qs2 )
{
	return(qs1.get_energy()<qs2.get_energy());
};
inline bool operator< ( const QuantumState<double>& qs1, const QuantumState<double>& qs2 )
{
	return(qs1.get_energy()<qs2.get_energy());
};

template <typename FLOAT>
std::ostream &operator<<(std::ostream& stream, const QuantumState<FLOAT>& qs);   // Screen output of a QuantumState


#endif