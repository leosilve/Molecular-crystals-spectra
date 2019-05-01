/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _MPSTATE_H
#define _MPSTATE_H

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <complex>
#include <vector>

///////////////////////////////////////////////////////////////////////////
/////////          Block diagonal Basis set

template <class FLOAT>
class MPSTATE {
private:
	std::vector<int> state;
	std::vector<std::complex<FLOAT> > coeff;
	FLOAT energy;
public:	
	MPSTATE() {
		state.resize(0);
		coeff.resize(0);
		energy=0.0;
	};

	MPSTATE(int _size) {
		state.resize(_size);
		coeff.resize(_size);
		energy=0.0;
	}

	void push_back(int _i, std::complex<FLOAT> _c) {state.push_back(_i); coeff.push_back(_c); };
	void resize(int _size) { state.resize(_size); coeff.resize(_size); };
	void set_energy(FLOAT _energy) { energy=_energy; };
	void set_state(int _pos, int _state) { 
		if (_pos>=0 && _pos < state.size() ) 
			state[_pos]=_state;
	};
	void set_coeff(int _pos, std::complex<FLOAT> _coeff) { 
		if (_pos>=0 && _pos < coeff.size() ) 
			coeff[_pos]=_coeff;
	};
	void erase(int _pos) {
		if (_pos>=0 && _pos < coeff.size() )
		{
			state.erase(state.begin()+_pos);
			coeff.erase(coeff.begin()+_pos);
		}
	};
	// Normalizes MPSTATE and discards negligible components
	bool normalize();
	
	int					get_state(int _i)	const { return(state[_i]); };
	std::complex<FLOAT> get_coeff(int _i)	const { return(coeff[_i]); };
	FLOAT				get_energy()		const { return(energy); };
	int size() const { return(state.size()); };

	const MPSTATE operator + (const MPSTATE&) const;
	const MPSTATE operator * (const FLOAT) const;
	const MPSTATE operator * (const std::complex<FLOAT>) const;
};

template <typename FLOAT>
std::ostream &operator<<(std::ostream& stream, const MPSTATE<FLOAT>& mps);

#endif
