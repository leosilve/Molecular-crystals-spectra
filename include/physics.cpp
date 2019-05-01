/* Written by Leonardo Silvestri 2007-2013 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/special_functions/factorials.hpp>  
#include <boost/math/special_functions/hermite.hpp>  
#include <boost/math/special_functions/binomial.hpp>  
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include <N_Functions.h>   

namespace ublas = boost::numeric::ublas;
namespace math = boost::math;

#ifndef ZERO_TOL
#define ZERO_TOL 1e-10
#endif

template <typename FLOAT>
FLOAT Boltzmann_factor(FLOAT dE, FLOAT T)
{
	if (T<1e-5)
	{
		if (std::abs(dE)<1e-5)
			return(1.0);
		else
			return(0.0);
	}
//	FLOAT kB=8.617343*1e-5;  // in ev K^-1
	return(exp(-dE/(8.617343*1e-5*T)));
}

// If lengths are in Angstroems and dipoles in Debye the result is in eV
template <typename FLOAT>
FLOAT dipole_dipole(ublas::vector<FLOAT> d1, ublas::vector<FLOAT> d2, ublas::vector<FLOAT> r)
{
		return(
				pow(norm_2(r),-5)*(inner_prod(d1,d2)*pow(norm_2(r),2)-3.0*inner_prod(d1,r)*inner_prod(d2,r))/1.602
			);
}

template <typename FLOAT>
FLOAT screened_dipole_dipole(ublas::vector<FLOAT> d1, ublas::vector<FLOAT> d2, ublas::vector<FLOAT> r, 
                             ublas::matrix<FLOAT, ublas::column_major> &epsilon_0)
{
		for (int i=0; i<3; i++)
		{
			r[i]/=sqrt(epsilon_0(i,i));
			d1[i]/=(pow(epsilon_0(0,0)*epsilon_0(1,1)*epsilon_0(2,2),FLOAT(0.25))
				                             *sqrt(epsilon_0(i,i)));
			d2[i]/=(pow(epsilon_0(0,0)*epsilon_0(1,1)*epsilon_0(2,2),FLOAT(0.25))
				                             *sqrt(epsilon_0(i,i)));
		}
		return(
					pow(norm_2(r),-5)*(inner_prod(d1,d2)*pow(norm_2(r),2)-3.0*inner_prod(d1,r)*inner_prod(d2,r))/1.602
		);
}

// FC overlap between shifted harmonic oscillator wavefunctions with the same energy, nu refers to the shifted potential
template <typename FLOAT>
FLOAT S_factor(int nu, int mu, FLOAT lambda)
{
		FLOAT sum=0.0;
	for (int i=0; i<=std::min(mu,nu); i++)
			sum+=pow(-1.0,nu-i)*pow(lambda,mu+nu-2*i)*math::factorial<FLOAT>(mu)*math::factorial<FLOAT>(nu)/(math::factorial<FLOAT>(i)*math::factorial<FLOAT>(mu-i)*math::factorial<FLOAT>(nu-i));
		return(
					exp(-lambda*lambda/2.0)/sqrt(math::factorial<FLOAT>(mu)*math::factorial<FLOAT>(nu))*sum
			);
}


// Computes the overlap between the nug(th) state of a harmonic oscillator of energy hwg and 
// the nue(th) state of a harmonic oscillator of energy hwe displaced by lambdag*sqrt(2*hbar/(m*wg)), 
// i.e. the displacement lambdag is expressed in the natural units of the excited oscillator
template <typename FLOAT>
double FCF(int nug, FLOAT hwg, int nue, FLOAT hwe, FLOAT lambdae)
{
	FLOAT sum=0.0;
	int kg, ke;

		for (kg=0; kg<=nug; kg++)	for (ke=0; ke<=nue; ke++)
		{
			if ((kg+ke)%2==0)
				sum+=math::binomial_coefficient<FLOAT>(nug,kg)/sqrt(math::factorial<FLOAT>(nug))*
				         math::binomial_coefficient<FLOAT>(nue,ke)/sqrt(math::factorial<FLOAT>(nue))*pow(2.0,kg+ke-(nug+nue)/2.0)/pow(hwg+hwe,(kg+ke)/2)*
				         pow(sqrt(hwg),kg)*math::hermite<FLOAT>(nug-kg, -hwe*lambdae*sqrt(2.0*hwg/hwe)/(hwg+hwe))*
				         pow(sqrt(hwe),ke)*math::hermite<FLOAT>(nue-ke, hwg*lambdae*sqrt(2.0)/(hwg+hwe))*
						 math::double_factorial<FLOAT>(kg+ke-1);
		}
		sum*=sqrt(2.0*sqrt(hwg*hwe)/(hwg+hwe)*exp(-hwg*2.0*lambdae*lambdae/(hwg+hwe)));
	return(sum);
}

///////////////////////////////
// Explicit instantiation

template double Boltzmann_factor(double dE, double T);
template double dipole_dipole(ublas::vector<double> d1, ublas::vector<double> d2, ublas::vector<double> r);
template double screened_dipole_dipole(ublas::vector<double> d1, ublas::vector<double> d2, ublas::vector<double> r, 
                             ublas::matrix<double, ublas::column_major> &epsilon_0);
template double S_factor(int nu, int mu, double lambda);
template double FCF(int nug, double hwg, int nue, double hwe, double lambdae);

template float Boltzmann_factor(float dE, float T);
template float dipole_dipole(ublas::vector<float> d1, ublas::vector<float> d2, ublas::vector<float> r);
template float screened_dipole_dipole(ublas::vector<float> d1, ublas::vector<float> d2, ublas::vector<float> r, 
									  ublas::matrix<float, ublas::column_major> &epsilon_0);
template float S_factor(int nu, int mu, float lambda);
template double FCF(int nug, float hwg, int nue, float hwe, float lambdae);


/* routines used by the Double Well problem
 double product(int start, int end, int step)
 {
 double result=1.0;
 for (int i=start; i<=end; i+=step)
 result*=i;
 return(result);
 }
 
 double product_sqrts(int start, int end, int step)
 {
 double result=1.0;
 for (int i=start; i<=end; i+=step)
 result*=sqrt(i);
 return(result);
 }
 void store_prod(int start, int end, int step, vector<int> &vec)
 {
 for (int i=start; i<=end; i+=step)
 vec.push_back(i);
 return;
 }
 
 double simplify_prods(vector<int> &num, vector<int> &den)
 {
 vector<int>::iterator it_num, it_den;
 sort(num.rbegin(), num.rend());
 sort(den.rbegin(), den.rend());
 //	cout << "num " << num << endl;
 //	cout << "den " << den << endl;
 int flag_begin=0;
 int flag_end=0;
 for (it_num=num.begin(); it_num!=num.end(); it_num++)
 {
 if (flag_begin==1)
 {
 it_num--;
 flag_begin=0;
 }
 for (it_den=den.begin(); it_den!=den.end(); it_den++)
 if (*it_num==*it_den)
 {
 num.erase(it_num);
 if (it_num!=num.begin())	
 it_num--;
 else
 flag_begin=1;
 if (it_num==num.end())	flag_end=1;
 den.erase(it_den);
 //				cout << "num " << num << endl;
 //				cout << "den " << den << endl;
 break;
 }
 if (flag_end==1) break;
 }
 
 double prod=1.0;
 int i;
 for (i=0; i<min(num.size(),den.size()); i++)
 prod*=double(num[i])/double(den[i]);
 if (num.size()>den.size())
 {
 for (i=den.size(); i<num.size(); i++)
 prod*=double(num[i]);
 }
 else
 {
 for (i=num.size(); i<den.size(); i++)
 prod/=double(den[i]);
 }
 return(prod);
 }*/


