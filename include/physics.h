/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace ublas = boost::numeric::ublas;

template <typename FLOAT>
FLOAT Boltzmann_factor(FLOAT dE, FLOAT T);

// If lengths are in Angstroems and dipoles in Debye the result is in eV
	template <typename FLOAT>
FLOAT dipole_dipole(ublas::vector<FLOAT> d1, ublas::vector<FLOAT> d2, ublas::vector<FLOAT> r);

template <typename FLOAT>
FLOAT screened_dipole_dipole(ublas::vector<FLOAT> d1, ublas::vector<FLOAT> d2, ublas::vector<FLOAT> r, 
		                           ublas::matrix<FLOAT, ublas::column_major> &epsilon_0);

// FC overlap between shifted harmonic oscillator wavefunctions with the same energy, nu refers to the shifted potential
template <typename FLOAT>
FLOAT S_factor(int nu, int mu, FLOAT lambda);

// Computes the overlap between the nug(th) state of a harmonic oscillator of energy hwg and 
// the nue(th) state of a harmonic oscillator of energy hwe displaced by lambdag*sqrt(2*hbar/(m*wg)), 
// i.e. the displacement lambdag is expressed in the natural units of the excited oscillator
template <typename FLOAT>
double FCF(int nug, FLOAT hwg, int nue, FLOAT hwe, FLOAT lambdae0);

#endif /* _PHYSICS_H_ */