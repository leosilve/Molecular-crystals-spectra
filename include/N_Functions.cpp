/* Written by Leonardo Silvestri 2007-2013 */

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

#include <blaswrap.h>
#include <f2c.h>
#include <clapack.h>

#include <boost/math/special_functions/factorials.hpp> 
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>

#include <OPutils.h> 
#include <N_Functions.h>

#ifndef Pi
#define Pi 3.141592654
#endif
#ifndef ZERO_TOL
#define ZERO_TOL 1e-10  // Numbers smaller than this are considered equal to 0
#endif

typedef std::complex<double>	Complex_d;
typedef std::complex<float>		Complex_f;

namespace ublas = boost::numeric::ublas; 
namespace math = boost::math;

// Numero di combinazioni di n elementi a k a k con ripetizioni
double CRnk_leo(int n, int k)
{
	return( math::factorial<double>(n+k-1)/math::factorial<double>(k)/math::factorial<double>(n-1) );
}

// Numero di disposizioni di n elementi a k a k senza ripetizioni
double Dnk_leo(int n, int k)
{
	if (n<k)	return(0.0);
	return( math::factorial<double>(n)/math::factorial<double>(n-k) );
}

template<typename MYFLOAT, typename MYDATA>
void mylinint(std::vector<MYFLOAT>& xa, std::vector<MYDATA>& ya, int n, MYFLOAT x, MYDATA *y)
{
	if (x<xa[0] || x>xa[n-1]) nrerror("mylinint: x out of bounds!");
	if (x==xa[0])
	{
		*y=ya[0];
		return;
	}
	if (x==xa[n-1])
	{
		*y=ya[n-1];
		return;
	}
	int i,ns=0;
	MYFLOAT dif,dift;
	dif=fabs(x-xa[0]);
	for (i=1;i<n;i++) 
	{
		if ( (dift=fabs(x-xa[i])) < dif) 
		{
			ns=i;
			dif=dift;
		}
	}
	if (x>xa[ns]) 	
		*y=ya[ns]+(ya[ns+1]-ya[ns])/(xa[ns+1]-xa[ns])*(x-xa[ns]);
	else
		*y=ya[ns-1]+(ya[ns]-ya[ns-1])/(xa[ns]-xa[ns-1])*(x-xa[ns-1]);
	return;
}

/*
// On exit:
// singular values are stored in s;
// range contains U, whose columns are left eigenvectors and (for positive singular values) form a basis for the range;
// nullspace contains V,  whose columns are right eigenvectors and (for zero singular values) form a basis for the null space;
template<typename T>
int CSVD(ublas::matrix<complex<T>,ublas::column_major> &a, vector<T> &s, 
		  ublas::matrix<complex<T>,ublas::column_major> &range,
		  ublas::matrix<complex<T>,ublas::column_major> &nullspace)
{
	char opt='O'; // O: optimal? M=minimal?
	char jobu='S';
	char jobvt='S';
	const int ierr = lapack::gesvd(opt,jobu, jobvt, a, s, range, nullspace) ;
	nullspace=herm(nullspace);
	if (ierr != 0)
	{
		std::cout << "CSVD: routine lapack::gesvd gave error code " << ierr << std::endl;
		nrerror("lapack::gesvd gave error code!");
	}
	int rank=0;
	for (unsigned int i=0; i<s.size(); i++)
		if (s[i]>ZERO_TOL)
			rank++;
	return(rank);
}


// w contains the diagonal of D, above. w must always be complex. (output)
// vr is an N x N matrix containing the right ("usual") eigenvectors of a in its columns.
template<typename T>
void RealMatrixEigendecomposition(ublas::matrix<T,ublas::column_major> &a, vector<complex<T> > &w, 
								  ublas::matrix<complex<T>,ublas::column_major> *vr)
{
	ublas::matrix<complex<T>,ublas::column_major>*	vl=NULL;  // This is because we don't want left eigenvectors
	const int ierr = lapack::geev(a, w, vl, vr, lapack::optimal_workspace());
	if (ierr != 0)
	{
		std::cout << "geev_cpp: routine gave error code " << ierr << std::endl;
		nrerror("geev_cpp gave error code!");
	}
	return;
}*/

#define min std::min

// On exit:
// singular values are stored in s;
// range contains U, whose columns are left eigenvectors and (for positive singular values) form a basis for the range;
// nullspace contains V,  whose columns are right eigenvectors and (for zero singular values) form a basis for the null space;
int CSVD(ublas::matrix<std::complex<double>,ublas::column_major> &a, std::vector<double> &s, 
		 ublas::matrix<std::complex<double>,ublas::column_major> &range,
		 ublas::matrix<std::complex<double>,ublas::column_major> &nullspace)
{

	int i, j;
	int m = a.size1();
	int n = a.size2();
	int minsize = min(n,m);
	
	std::vector<Complex_d> avec(n*m);
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			avec[j*m+i]=a(i,j);
	
	int ldu = a.size1();
	
	std::vector<Complex_d> u(m*minsize);
	std::vector<Complex_d> vt(n*minsize);
	
	gesvd_cpp(m, n, avec, s, u, vt);

	range.resize(m,minsize);
	nullspace.resize(minsize,n);
	
				 
	for (i=0; i<m; i++)
		for (j=0; j<minsize; j++)
			range(i,j)=u[j*m+i];
	for (i=0; i<minsize; i++)
		for (j=0; j<n; j++)
			nullspace(i,j)=vt[j*minsize+i];

	nullspace=herm(nullspace);
	int rank=0;
	for (unsigned int i=0; i<s.size(); i++)
		if (s[i]>ZERO_TOL)
			rank++;
	return(rank);
}

	/*  ZGESVD computes the singular value decomposition (SVD) of a complex */
	/*  M-by-N matrix A, optionally computing the left and/or right singular */
	/*  vectors. The SVD is written */
		
	/*       A = U * SIGMA * conjugate-transpose(V) */
		
	/*  where SIGMA is an M-by-N matrix which is zero except for its */
	/*  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and */
	/*  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA */
	/*  are the singular values of A; they are real and non-negative, and */
	/*  are returned in descending order.  The first min(m,n) columns of */
	/*  U and V are the left and right singular vectors of A. */
		
	/*  Note that the routine returns V**H, not V. */
		/*  Arguments */
	/*  ========= */
		
	/*  JOBU    (input) CHARACTER*1 */
	/*          Specifies options for computing all or part of the matrix U: */
	/*          = 'A':  all M columns of U are returned in array U: */
	/*          = 'S':  the first min(m,n) columns of U (the left singular */
	/*                  vectors) are returned in the array U; */
	/*          = 'O':  the first min(m,n) columns of U (the left singular */
	/*                  vectors) are overwritten on the array A; */
	/*          = 'N':  no columns of U (no left singular vectors) are */
	/*                  computed. */
		
	/*  JOBVT   (input) CHARACTER*1 */
	/*          Specifies options for computing all or part of the matrix */
	/*          V**H: */
	/*          = 'A':  all N rows of V**H are returned in the array VT; */
	/*          = 'S':  the first min(m,n) rows of V**H (the right singular */
	/*                  vectors) are returned in the array VT; */
	/*          = 'O':  the first min(m,n) rows of V**H (the right singular */
	/*                  vectors) are overwritten on the array A; */
	/*          = 'N':  no rows of V**H (no right singular vectors) are */
	/*                  computed. */
		
	/*          JOBVT and JOBU cannot both be 'O'. */
		
	/*  M       (input) INTEGER */
	/*          The number of rows of the input matrix A.  M >= 0. */
		
	/*  N       (input) INTEGER */
	/*          The number of columns of the input matrix A.  N >= 0. */
		
	/*  A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
	/*          On entry, the M-by-N matrix A. */
	/*          On exit, */
	/*          if JOBU = 'O',  A is overwritten with the first min(m,n) */
	/*                          columns of U (the left singular vectors, */
	/*                          stored columnwise); */
	/*          if JOBVT = 'O', A is overwritten with the first min(m,n) */
	/*                          rows of V**H (the right singular vectors, */
	/*                          stored rowwise); */
	/*          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A */
	/*                          are destroyed. */
		
	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,M). */
		
	/*  S       (output) DOUBLE PRECISION array, dimension (min(M,N)) */
	/*          The singular values of A, sorted so that S(i) >= S(i+1). */
		
	/*  U       (output) COMPLEX*16 array, dimension (LDU,UCOL) */
	/*          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'. */
	/*          If JOBU = 'A', U contains the M-by-M unitary matrix U; */
	/*          if JOBU = 'S', U contains the first min(m,n) columns of U */
	/*          (the left singular vectors, stored columnwise); */
	/*          if JOBU = 'N' or 'O', U is not referenced. */
		
	/*  LDU     (input) INTEGER */
	/*          The leading dimension of the array U.  LDU >= 1; if */
	/*          JOBU = 'S' or 'A', LDU >= M. */
		
	/*  VT      (output) COMPLEX*16 array, dimension (LDVT,N) */
	/*          If JOBVT = 'A', VT contains the N-by-N unitary matrix */
	/*          V**H; */
	/*          if JOBVT = 'S', VT contains the first min(m,n) rows of */
	/*          V**H (the right singular vectors, stored rowwise); */
	/*          if JOBVT = 'N' or 'O', VT is not referenced. */
		
	/*  LDVT    (input) INTEGER */
	/*          The leading dimension of the array VT.  LDVT >= 1; if */
	/*          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N). */
		
	/*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
	/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
		
	/*  LWORK   (input) INTEGER */
	/*          The dimension of the array WORK. */
	/*          LWORK >=  MAX(1,2*MIN(M,N)+MAX(M,N)). */
	/*          For good performance, LWORK should generally be larger. */
		
	/*          If LWORK = -1, then a workspace query is assumed; the routine */
	/*          only calculates the optimal size of the WORK array, returns */
	/*          this value as the first entry of the WORK array, and no error */
	/*          message related to LWORK is issued by XERBLA. */
		
	/*  RWORK   (workspace) DOUBLE PRECISION array, dimension (5*min(M,N)) */
	/*          On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains the */
	/*          unconverged superdiagonal elements of an upper bidiagonal */
	/*          matrix B whose diagonal is in S (not necessarily sorted). */
	/*          B satisfies A = U * B * VT, so it has the same singular */
	/*          values as A, and singular vectors related by U and VT. */
		
	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit. */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
	/*          > 0:  if ZBDSQR did not converge, INFO specifies how many */
	/*                superdiagonals of an intermediate bidiagonal form B */
	/*                did not converge to zero. See the description of RWORK */
/*                above for details. */


void gesvd_cpp(int dim1, int dim2, std::vector<Complex_d> &a, std::vector<double> &s, std::vector<Complex_d> &u, 
			   std::vector<Complex_d> &vt)
{
	char jobu = 'S';
	char jobvt = 'S';
	long int m = dim1;
	long int n = dim2;
	long int lda = m;
	
	long int ldvt = min(n,m);	//	The leading dimension of the array VT.  LDVT >= 1; if */
									//  JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
	long int ldu = m;				//  The leading dimension of the array U.  LDU >= 1; if */
									//          JOBU = 'S' or 'A', LDU >= M. */
	long int info;
	long int lwork = 4*n;							// The dimension of the array WORK.  LWORK >=  MAX(1,2*MIN(M,N)+MAX(M,N)). 
	std::vector<Complex_d> work(lwork);				// (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
	std::vector<double> rwork(5*min(n,m));     //	DOUBLE PRECISION array, dimension (5*min(M,N))
	
	s.resize(min(n,m));
	u.resize(ldu*min(n,m));
	vt.resize(ldvt*n);
	
	const int ierr = zgesvd_(&jobu, &jobvt, &m, &n, 
							 reinterpret_cast<doublecomplex*>(&a[0]), &lda, &s[0], reinterpret_cast<doublecomplex*>(&u[0]), 
							 &ldu, reinterpret_cast<doublecomplex*>(&vt[0]), &ldvt, reinterpret_cast<doublecomplex*>(&work[0]), 
							 &lwork, &rwork[0], &info);
	if (ierr != 0)
	{
		std::cout << "ComplexMatrixEigendecomposition: routine zgesvd_ gave error code " << ierr << std::endl;
		nrerror("gesvd_cpp gave error code!");
	}
	return;
}

#undef min

/*  ZGESV computes the solution to a complex system of linear equations */
/*     A * X = B, */
/*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices. */

/*  The LU decomposition with partial pivoting and row interchanges is */
/*  used to factor A as */
/*     A = P * L * U, */
/*  where P is a permutation matrix, L is unit lower triangular, and U is */
/*  upper triangular.  The factored form of A is then used to solve the */
/*  system of equations A * X = B. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The number of linear equations, i.e., the order of the */
/*          matrix A.  N >= 0. */

/*  NRHS    (input) INTEGER */
/*          The number of right hand sides, i.e., the number of columns */
/*          of the matrix B.  NRHS >= 0. */

/*  A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*          On entry, the N-by-N coefficient matrix A. */
/*          On exit, the factors L and U from the factorization */
/*          A = P*L*U; the unit diagonal elements of L are not stored. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  IPIV    (output) INTEGER array, dimension (N) */
/*          The pivot indices that define the permutation matrix P; */
/*          row i of the matrix was interchanged with row IPIV(i). */

/*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS) */
/*          On entry, the N-by-NRHS matrix of right hand side matrix B. */
/*          On exit, if INFO = 0, the N-by-NRHS solution matrix X. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B.  LDB >= max(1,N). */
void gesv_cpp(int dim, int dimrhs, std::vector<Complex_d> &a, std::vector<Complex_d> &b)
{
	long int n = dim;
	long int nrhs=dimrhs;
	long int lda = n;
	long int ldb = n;
	std::vector<long int> ipiv(n); 

	long int info;

	if (a.size()!=n*n || b.size()!=n*nrhs)
		nrerror("gesv_cpp: matrices have wrong dimensions!");
	
	const int ierr = zgesv_(&n, &nrhs, reinterpret_cast<doublecomplex*>(&a[0]), 
							&lda, &ipiv[0], reinterpret_cast<doublecomplex*>(&b[0]), &ldb, &info);
	if (ierr != 0)
	{
		std::cout << "ComplexMatrixEigendecomposition: routine zgesv_ gave error code " << ierr << std::endl;
		nrerror("gesv_cpp gave error code!");
	}
	return;
}

bool InvertComplexMatrix(const ublas::matrix<std::complex<double>, ublas::column_major>& input, ublas::matrix<std::complex<double>, ublas::column_major>& inverse) 
{
	if (input.size1()!=input.size2())
		nrerror("InvertComplexMatrix: cannot invert rectangular matrices!");
	int n =	input.size1();
	std::vector<Complex_d> a(n*n);
	std::vector<Complex_d> b(n*n, 0.0);
	int i,j;
	for (i=0; i<n; i++)
		b[i*n+i]=1.0;

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			a[j*n+i]=input(i,j);
	gesv_cpp(n, n, a, b);

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			inverse(i,j)=b[j*n+i];
	return true;
}

void ComplexMatrixEigendecomposition(ublas::matrix<Complex_d , ublas::column_major> &a, std::vector<Complex_d> &w, 
									 ublas::matrix<Complex_d , ublas::column_major> *vr)
{
	int n = a.size1();
	std::vector<Complex_d> avec(n*n);
	std::vector<Complex_d> vrvec(n*n);

	int i, j;
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			avec[j*n+i]=a(i,j);
	
	geev_cpp(n, avec, w, vrvec);
	
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			(*vr)(i,j)=vrvec[j*n+i];
	return;
}

void geev_cpp(int dim, std::vector<Complex_d> &a, std::vector<Complex_d> &w, std::vector<Complex_d> &vr)
{
	char jobvl = 'N';
	char jobvr = 'V';
	long int n = dim;
	long int lda = n;
	long int ldvl = 1;
	long int ldvr = n;
	long int info;
	w.resize(n); 
	doublecomplex*	vl=NULL;						// This is because we don't want left eigenvectors
	long int lwork = 3*n;							// The dimension of the array WORK.  LWORK >= max(1,2*N). 
	std::vector<Complex_d> work(lwork);					// (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
	std::vector<double> rwork(max(1,3*dim-2));     //	DOUBLE PRECISION array, dimension (max(1, 3*N-2))
	
	const int ierr = zgeev_(&jobvl, &jobvr, &n, reinterpret_cast<doublecomplex*>(&a[0]), &lda, 
							reinterpret_cast<doublecomplex*>(&w[0]), vl, &ldvl, 
							reinterpret_cast<doublecomplex*>(&vr[0]), &ldvr, 
							reinterpret_cast<doublecomplex*>(&work[0]), 
							&lwork, &rwork[0], &info);
	if (ierr != 0)
	{
		std::cout << "ComplexMatrixEigendecomposition: routine zgeev_ gave error code " << ierr << std::endl;
		nrerror("geev gave error code!");
	}
	return;
}
	
// Eigendecomposition of a complex Hermitian matrix A = Q * D * Q'
// If jobz='V' on exit matrix a contains the eigenvectors (one in each column)
void heev_cpp(char jobz, char uplo, int dim, std::vector<Complex_d> &a, 
			  std::vector<double> &w, char workspace)
{
	long int n = dim;
	long int lda = n;
	long int info;
	w.resize(n);                            //	DOUBLE PRECISION array, dimension (N)
											//	If INFO = 0, the eigenvalues in ascending order.
	std::vector<Complex_d> work(1);			//	COMPLEX*16 array, dimension (MAX(1,LWORK))
	std::vector<double> rwork(max(1,3*dim-2));     //	DOUBLE PRECISION array, dimension (max(1, 3*N-2))

	long int lwork;							//	The length of the array WORK.  LWORK >= max(1,2*N-1).
											//	If LWORK = -1, then a workspace query is assumed; the routine
											//	only calculates the optimal size of the WORK array, returns
											//	this value as the first entry of the WORK array, and no error
											//	message related to LWORK is issued by XERBLA.
	
	//	This is to determine the optimal size	
	if (workspace=='O' || workspace=='o')
	{
		lwork = -1;
		zheev_(&jobz, &uplo, &n, reinterpret_cast<doublecomplex*>(&a[0]), &lda, &w[0], 
			   reinterpret_cast<doublecomplex*>(&work[0]), &lwork, &rwork[0], &info);     
		lwork = (int)work[0].real();
	}
	else
		lwork = max(1,2*dim-1);

	
	//	This is to compute eigenvalues and eigenvectors
	work.resize(lwork);
	const int ierr = zheev_(&jobz, &uplo, &n, reinterpret_cast<doublecomplex*>(&a[0]), &lda, &w[0], 
							reinterpret_cast<doublecomplex*>(&work[0]), &lwork, &rwork[0], &info);     
	if (ierr != 0)
	{
		std::cout << "heev_cpp: routine gave error code " << ierr << std::endl;
		nrerror("heev_cpp gave error code!");
	}
	return;
}

void heev_cpp(char jobz, char uplo, int dim, std::vector<Complex_f>&a, 
			  std::vector<float> &w, char workspace)
{
	long int n = dim;
	long int lda = n;
	long int info;
	w.resize(n);                            //	float PRECISION array, dimension (N)
											//	If INFO = 0, the eigenvalues in ascending order.
	std::vector<Complex_f> work(1);				//	COMPLEX*16 array, dimension (MAX(1,LWORK))
	std::vector<float> rwork(max(1,3*n-2));     //	float PRECISION array, dimension (max(1, 3*N-2))

	long int lwork;							//	The length of the array WORK.  LWORK >= max(1,2*N-1).
											//	If LWORK = -1, then a workspace query is assumed; the routine
											//	only calculates the optimal size of the WORK array, returns
											//	this value as the first entry of the WORK array, and no error
											//	message related to LWORK is issued by XERBLA.
	
	//	This is to determine the optimal size	
	if (workspace=='O' || workspace=='o')
	{
		lwork = -1;
		cheev_(&jobz, &uplo, &n, reinterpret_cast<complex*>(&a[0]), &lda, &w[0], 
			   reinterpret_cast<complex*>(&work[0]), &lwork, &rwork[0], &info);     
		lwork = (int)work[0].real();
	}
	else
		lwork = max(1,2*n-1);
	
	
	//	This is to compute eigenvalues and eigenvectors
	work.resize(lwork);
	const int ierr = cheev_(&jobz, &uplo, &n, reinterpret_cast<complex*>(&a[0]), &lda, &w[0], 
							reinterpret_cast<complex*>(&work[0]), &lwork, &rwork[0], &info);     
	if (ierr != 0)
	{
		std::cout << "heev_cpp: routine gave error code " << ierr << std::endl;
		nrerror("heev_cpp gave error code!");
	}
	return;
}


void syev_cpp(char jobz, char uplo, int dim, std::vector<float>&a, 
			  std::vector<float> &w, char workspace)
{
	long int n = dim;
	long int lda = n;
	long int info;
	w.resize(n);                            //	float PRECISION array, dimension (N)
	//	If INFO = 0, the eigenvalues in ascending order.
	std::vector<float> work(1);  //	COMPLEX*16 array, dimension (MAX(1,LWORK))

	long int lwork;							//	The length of the array WORK.  LWORK >= max(1,2*N-1).
	//	If LWORK = -1, then a workspace query is assumed; the routine
	//	only calculates the optimal size of the WORK array, returns
	//	this value as the first entry of the WORK array, and no error
	//	message related to LWORK is issued by XERBLA.
	
	//	This is to determine the optimal size	
	if (workspace=='O' || workspace=='o')
	{
		lwork = -1;
		ssyev_(&jobz, &uplo, &n, &a[0], &lda, &w[0], &work[0], &lwork, &info);     
		lwork = (int)work[0];
	}
	else
		lwork = max(1,2*n-1);
	
	
	//	This is to compute eigenvalues and eigenvectors
	work.resize(lwork);
	const int ierr = ssyev_(&jobz, &uplo, &n, &a[0], &lda, &w[0], &work[0], &lwork, &info); 
	
	if (ierr != 0)
	{
		std::cout << "syev_cpp: routine gave error code " << ierr << std::endl;
		nrerror("syev_cpp gave error code!");
	}
	return;
}


void syev_cpp(char jobz, char uplo, int dim, std::vector<double>&a, 
			  std::vector<double> &w, char workspace)
{
	long int n = dim;
	long int lda = n;
	long int info;
	w.resize(n);                            //	double PRECISION array, dimension (N)
	//	If INFO = 0, the eigenvalues in ascending order.
	std::vector<double> work(1);  //	COMPLEX*16 array, dimension (MAX(1,LWORK))
//	std::vector<double> rwork(max(1,3*n-2));     //	double PRECISION array, dimension (max(1, 3*N-2))
	
	long int lwork;							//	The length of the array WORK.  LWORK >= max(1,2*N-1).
	//	If LWORK = -1, then a workspace query is assumed; the routine
	//	only calculates the optimal size of the WORK array, returns
	//	this value as the first entry of the WORK array, and no error
	//	message related to LWORK is issued by XERBLA.
	
	//	This is to determine the optimal size	
	if (workspace=='O' || workspace=='o')
	{
		lwork = -1;
		dsyev_(&jobz, &uplo, &n, &a[0], &lda, &w[0], &work[0], &lwork,  &info);     
		lwork = (int)work[0];
	}
	else
		lwork = max(1,2*n-1);
	
	
	//	This is to compute eigenvalues and eigenvectors
	work.resize(lwork);
	const int ierr = dsyev_(&jobz, &uplo, &n, &a[0], &lda, &w[0], &work[0], &lwork,  &info); 
	
	if (ierr != 0)
	{
		std::cout << "syev_cpp: routine gave error code " << ierr << std::endl;
		nrerror("syev_cpp gave error code!");
	}
	return;
}

template<typename T>
T Gaussian_dist(T dx,T width)
{
	return(
		   1.0/(width * sqrt(2.0*Pi) ) * exp( - pow ( (dx)/width , 2) /2.0 )
		   );
}

/* Matrix inversion routine.
   Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const ublas::matrix<T, ublas::column_major>& input, ublas::matrix<T, ublas::column_major>& inverse) 
{
 	using namespace boost::numeric::ublas;
 	typedef permutation_matrix<std::size_t> pmatrix;
 	// create a working copy of the input
 	matrix<T> A(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());

 	// perform LU-factorization
 	int res = lu_factorize(A,pm);
        if( res != 0 ) return false;

 	// create identity matrix of "inverse"
 	inverse.assign(ublas::identity_matrix<T>(A.size1()));

 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);
 	return true;
 }

template<class T>
const ublas::vector<T>	cross_prod(const ublas::vector<T>& v1, const ublas::vector<T>& v2) 
{
	if (v1.size() != 3 || v2.size() != 3)
		nrerror("cross_prod: wrong vector size!");
	ublas::vector<T> res(3);
	res(0)=v1(1)*v2(2)-v2(1)*v1(2);
	res(1)=-v1(0)*v2(2)+v2(0)*v1(2);
	res(2)=v1(0)*v2(1)-v2(0)*v1(1);
	return(res);
}

/*
template<class T>
inline bool operator<(const ublas::vector<T>& v1, const ublas::vector<T>& v2) const
{
	if (v1.size() != v2.size())
		nrerror("operator <: wrong comparison between ublas vectors!");
	else
		if (norm_2(v1)<norm_2(v2))
			return(true);
	return(false);
}

template<class T>
inline bool operator==(const ublas::vector<T>& v1, const ublas::vector<T>& v2) const
{
	if (v1.size() != v2.size())
		nrerror("operator ==: wrong comparison between ublas vectors!");
	else
		for (unsigned int i=0; i<v1.size(); i++)
			if (v1[i]!=v2[i])
				return(false);
	return(true);
}*/

template<class T>
bool Polynomial_regression(std::vector<T>& x, std::vector<T>& y, std::vector<T>& par)
{
	int i,j;
	long int dim=par.size();
	long int ndata=x.size();
	if (y.size()!=ndata)
		nrerror("Polynomial_regression: x and y have different sizes!");
//	std::cout << ndata << " " << dim << std::endl;
	ublas::matrix<T, ublas::column_major> xmat(ndata,dim);
	ublas::matrix<T, ublas::column_major> invxmat(dim,dim);
	ublas::matrix<T, ublas::column_major> tempmat(dim,dim);
	ublas::matrix<T, ublas::column_major> transxmat(dim,ndata);
	for (i=0; i<ndata; i++) for (j=0; j<dim; j++)
		xmat(i,j)=pow(x[i],j);
//	std::cout << xmat << std::endl;
	transxmat=ublas::trans(xmat);
	tempmat=ublas::prod(transxmat,xmat);
	InvertMatrix(tempmat,invxmat);
	transxmat=ublas::prod(invxmat,transxmat);
	for (j=0; j<dim; j++) 
	{
		par[j]=0;
		for (i=0; i<ndata; i++)	par[j]+=transxmat(j,i)*y[i];
	}
	return true;
}
///////////////////////////////////////////////
// Explicit instantiation

//// double

template bool InvertMatrix(const ublas::matrix<double, ublas::column_major>& input, ublas::matrix<double, ublas::column_major>& inverse);
template const ublas::vector<double> cross_prod(const ublas::vector<double>& v1, const ublas::vector<double>& v2);
template double Gaussian_dist(double dx,double width);
template void mylinint(std::vector<double>& xa, std::vector<double>& ya, int n, double x, double *y);
template void mylinint(std::vector<double>& xa, std::vector<std::complex<double> >& ya, int n, double x, std::complex<double> *y);
template bool Polynomial_regression(std::vector<double>& x, std::vector<double>& y, std::vector<double>& par);
/*
template
int CSVD(ublas::matrix<complex<double>,ublas::column_major> &a, vector<double> &s, 
		 ublas::matrix<complex<double>,ublas::column_major> &range,
		 ublas::matrix<complex<double>,ublas::column_major> &nullspace);

template
void ComplexMatrixEigendecomposition(ublas::matrix<complex<double>,ublas::column_major> &a, vector<complex<double> > &w, 
									 ublas::matrix<complex<double>,ublas::column_major> *vr);

// w contains the diagonal of D, above. w must always be complex. (output)
// vr is an N x N matrix containing the right ("usual") eigenvectors of a in its columns.
template
void RealMatrixEigendecomposition(ublas::matrix<double,ublas::column_major> &a, vector<complex<double> > &w, 
								  ublas::matrix<complex<double>,ublas::column_major> *vr);*/

/*template
bool operator<(const ublas::vector<double>& v1, const ublas::vector<double>& v2) const;

template
bool operator==(const ublas::vector<double>& v1, const ublas::vector<double>& v2) const;
*/

//// float

template bool InvertMatrix(const ublas::matrix<float, ublas::column_major>& input, ublas::matrix<float, ublas::column_major>& inverse);
template const ublas::vector<float>	cross_prod(const ublas::vector<float>& v1, const ublas::vector<float>& v2);
template float Gaussian_dist(float dx,float width);
template void mylinint(std::vector<float>& xa, std::vector<float>& ya, int n, float x, float *y);
template void mylinint(std::vector<float>& xa, std::vector<std::complex<float> >& ya, int n, float x, std::complex<float> *y);
template bool Polynomial_regression(std::vector<float>& x, std::vector<float>& y, std::vector<float>& par);

/*
template
int CSVD(ublas::matrix<complex<float>,ublas::column_major> &a, vector<float> &s, 
		 ublas::matrix<complex<float>,ublas::column_major> &range,
		 ublas::matrix<complex<float>,ublas::column_major> &nullspace);

template
void ComplexMatrixEigendecomposition(ublas::matrix<complex<float>,ublas::column_major> &a, vector<complex<float> > &w, 
									 ublas::matrix<complex<float>,ublas::column_major> *vr);

template
void RealMatrixEigendecomposition(ublas::matrix<float,ublas::column_major> &a, vector<complex<float> > &w, 
								  ublas::matrix<complex<float>,ublas::column_major> *vr);*/

/*
template
bool operator<(const ublas::vector<float>& v1, const ublas::vector<float>& v2) const;

template
bool operator==(const ublas::vector<float>& v1, const ublas::vector<float>& v2) const;

*/



