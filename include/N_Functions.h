/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _N_FUNCTIONS_H
#define _N_FUNCTIONS_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>
#include <complex>

namespace ublas = boost::numeric::ublas;

// Numero di combinazioni di n elementi a k a k con ripetizioni
double CRnk_leo(int n, int k);

// Numero di disposizioni di n elementi a k a k senza ripetizioni
double Dnk_leo(int n, int k);

// Linear interpolation 
template<typename MYFLOAT, typename MYDATA>
void mylinint(std::vector<MYFLOAT>& xa, std::vector<MYDATA>& ya, int n, MYFLOAT x, MYDATA *y);


// On exit:
// singular values are stored in s;
// range contains U, whose columns are left eigenvectors and (for positive singular values) form a basis for the range;
// nullspace contains V,  whose columns are right eigenvectors and (for zero singular values) form a basis for the null space;

int CSVD(ublas::matrix<std::complex<double>,ublas::column_major> &a, std::vector<double> &s, 
		 ublas::matrix<std::complex<double>,ublas::column_major> &range,
		 ublas::matrix<std::complex<double>,ublas::column_major> &nullspace);
/*
template<typename T>
void ComplexMatrixEigendecomposition(ublas::matrix<complex<T>,ublas::column_major> &a, vector<complex<T> > &w, 
									 ublas::matrix<complex<T>,ublas::column_major> *vr);

// w contains the diagonal of D, above. w must always be complex. (output)
// vr is an N x N matrix containing the right ("usual") eigenvectors of a in its columns.
template<typename T>
void RealMatrixEigendecomposition(ublas::matrix<T,ublas::column_major> &a, vector<complex<T> > &w, 
								  ublas::matrix<complex<T>,ublas::column_major> *vr);
*/
void ComplexMatrixEigendecomposition(ublas::matrix<std::complex<double>,ublas::column_major> &a, std::vector<std::complex<double> > &w, 
									 ublas::matrix<std::complex<double>,ublas::column_major> *vr);
void geev_cpp(int dim, std::vector<std::complex<double> > &a, std::vector<std::complex<double> > &w, std::vector<std::complex<double> > &vr);
void gesvd_cpp(int dim1, int dim2, std::vector<std::complex<double> > &a, std::vector<double> &s, std::vector<std::complex<double> > &u, 
			    std::vector<std::complex<double> > &vt);
void gesv_cpp(int dim, int dimrhs, std::vector<std::complex<double> > &a, std::vector<std::complex<double> > &b);

// Eigendecomposition of a complex Hermitian matrix A = Q * D * Q'
void heev_cpp(char jobz, char uplo, int dim, std::vector<std::complex<float> >&a,  std::vector<float> &w, char workspace);
void heev_cpp(char jobz, char uplo, int dim, std::vector<std::complex<double> >&a,  std::vector<double> &w, char workspace);

//Eigendecomposition of a real symmetric matrix A = Q * D * Q'
void syev_cpp(char jobz, char uplo, int dim, std::vector<float>&a, std::vector<float> &w, char workspace);
void syev_cpp(char jobz, char uplo, int dim, std::vector<double>&a, std::vector<double> &w, char workspace);

template<typename T>
T  Gaussian_dist(T dx, T width);

 /* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const ublas::matrix<T, ublas::column_major>& input, ublas::matrix<T, ublas::column_major>& inverse);

bool InvertComplexMatrix(const ublas::matrix<std::complex<double>, ublas::column_major>& input, ublas::matrix<std::complex<double>, ublas::column_major>& inverse); 

template<class T>
bool Polynomial_regression(std::vector<T>& x, std::vector<T>& y, std::vector<T>& par);

template<class T>
const ublas::vector<T>	cross_prod(const ublas::vector<T>& v1, const ublas::vector<T>& v2);

inline bool operator<(const ublas::vector<float>& v1, const ublas::vector<float>& v2)
{
	if (v1.size() != v2.size())
		return(false);
	else
		if (norm_2(v1)<norm_2(v2))
			return(true);
	return(false);
}
inline bool operator<(const ublas::vector<double>& v1, const ublas::vector<double>& v2)
{
	if (v1.size() != v2.size())
		return(false);
	else
		if (norm_2(v1)<norm_2(v2))
			return(true);
	return(false);
}

inline bool operator==(const ublas::vector<float>& v1, const ublas::vector<float>& v2)
{
	if (v1.size() != v2.size())
		return(false);
	else
		for (unsigned int i=0; i<v1.size(); i++)
			if (v1[i]!=v2[i])
				return(false);
	return(true);
}

inline bool operator==(const ublas::vector<double>& v1, const ublas::vector<double>& v2)
{
	if (v1.size() != v2.size())
		return(false);
	else
		for (unsigned int i=0; i<v1.size(); i++)
			if (v1[i]!=v2[i])
				return(false);
	return(true);
}


#endif
