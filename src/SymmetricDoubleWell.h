/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _SYMMETRICDOUBLEWELL_H_
#define _SYMMETRICDOUBLEWELL_H_

#include <DEIntegrator.h>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/factorials.hpp>  

namespace math = boost::math;

#ifndef Pi
#define Pi 3.141592654
#endif 

class FunctionObject
{
public:
	FunctionObject(unsigned _n, unsigned _m, double _beta)
	{
		n = _n;
		m = _m;
		beta = _beta;
	}
	double operator()(double x) const
	{
		return(sqrt(sqrt(beta)/(Pi*pow(2.0,int(n+m-1))*math::factorial<double>(n)*math::factorial<double>(m)))
			   *math::hermite(n,sqrt(2.0*beta)*x)*exp(-beta*x*x)
			   *math::hermite(m,sqrt(2.0)*x)*exp(-x*x));
	}
private:
	unsigned n;
	unsigned m;
	double beta;
};

double find_genFCF(unsigned nug, unsigned nue, double beta, double precision);

// Computes the symmetric double well hamiltonian matrix element between two harmonic oscillator eigenfunctions 
// with quantum numbers m,n. The perturbation to the harmonic oscillator of energy hw0 is 
// a Gaussian with peak height A (in units of energy) and width inversely proportional to alpha (in units of energy)
template <typename FLOAT>
double Hmn(FLOAT hw0, FLOAT A, FLOAT alpha, int m, int n);

template <typename FLOAT>
void find_SDW_eigenvectors(FLOAT hw0, FLOAT A, FLOAT alpha, int n, 
						   std::vector<FLOAT> &eigenvalues, std::vector<std::vector<FLOAT> >  &eigenvectors);

void initialize_SDW_parameters(std::string DW_datafile, double& DW_hw0, double& DW_alpha, double& DW_A, double& DW_Vd,
							   double& DW_hwe, double& DW_lambdae, unsigned int& DW_basis_dim, unsigned int& DW_nug_max,
							   unsigned int& DW_nue_max, std::vector<std::vector<double> >& DW_FCF_table,
							   std::vector<double>&	DW_levels);

#endif /* _SYMMETRICDOUBLEWELL_H_ */