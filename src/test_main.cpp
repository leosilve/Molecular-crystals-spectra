#include <lattice.h>
//#include <N_Functions.h>
//#include <MultiPhonon.h>
#include <PhononCloud.h>
#include <OP_model.h>
//#include <OP_epsilon.h>
//#include <OP_absorption.h>
#include <OPutils.h>
#include <iostream>
#include <vector>
#include <complex>

#ifndef Pi
#define Pi 3.141592654
#endif 
#define MYFLOAT float

namespace ublas = boost::numeric::ublas;
using namespace std;

typedef ublas::matrix<MYFLOAT, ublas::column_major>	Matrix;
typedef ublas::vector<MYFLOAT>	Vector;

typedef std::complex<double>							dComplex;
typedef ublas::matrix<dComplex,	ublas::column_major>	dComplexMatrix;

int main (int argc, char* argv[]) {

	int i,j,m,rank;
	cout << "Testing ComplexMatrixEigendecomposition" << endl;
	
	dComplexMatrix CMtemp(4,4);  // My ComplexMatrix to run all the tests
	dComplexMatrix EigVecs(4,4);
	dComplexMatrix* EigVecsPtr=&EigVecs;
	vector<dComplex> kz(4);
	dComplex I(0.0,1.0);
	
	CMtemp(0,0)=1.0;
						CMtemp(1,1)=0.2;	CMtemp(1,2)=2.0*I;
						CMtemp(2,1)=1.0*I;	CMtemp(2,2)=-0.3;	CMtemp(2,3)=-0.1*I;
	CMtemp(3,0)=0.4;											CMtemp(3,3)=2.5;
	
	cout << "CMtemp: " << endl << CMtemp << endl;
	
	ComplexMatrixEigendecomposition(CMtemp, kz, EigVecsPtr);
	cout << "Eigenvalues : " << kz << endl;
	cout << "Eigenvectors (each row is an eigvec): " << endl;
	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
			cout << (*EigVecsPtr)(j,i) << " ";
		cout << endl;
	}
 
	cout << "Testing CSVD with a rectangular input matrix" << endl;
	
	vector<double>	SingularValues(4);
	dComplexMatrix	Eigensystem(4,3);
	dComplexMatrix	range(4,3);
	dComplexMatrix	null_space(3,3);


	for (m=0; m<4; m++) for (j=0; j<3; j++) 
		Eigensystem(m,j)=CMtemp(m,j); 

	cout << "Input Matrix: " << endl;
	for (i=0; i<4; i++)
	{
		for (j=0; j<3; j++)
			cout << Eigensystem(i,j) << " ";
		cout << endl;
	}
	rank=CSVD(Eigensystem, SingularValues, range, null_space);
	cout << "Rank : " << rank << endl;
	cout << "SingularValues : " << SingularValues << endl;
	cout << "Range (each row is a base vector if the corresponding singular value is >0): "  << endl;
	for (i=0; i<3; i++)
	{
		for (j=0; j<4; j++)
			cout << range(j,i) << " ";
		cout << endl;
	}
	cout << "Null space (each row is a base vector if the corresponding singular value is 0): "  << endl;
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
			cout << null_space(j,i) << " ";
		cout << endl;
	}
	
	
	cout << "Testing gesv_cpp to solve CMtemp.X==b" << endl;
	
	vector<dComplex> a(4*4);
	vector<dComplex> b(4*2);
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			a[j*4+i]=CMtemp(i,j);
	for (i=0; i<4; i++)
		for (j=0; j<2; j++)
			b[j*4+i]=Eigensystem(i,j)*(1.0+i+j);
	cout << "b: "  << endl;
	for (i=0; i<4; i++)
	{
		for (j=0; j<2; j++)
			cout << b[j*4+i] << " ";
		cout << endl;
	}
	
	gesv_cpp(4,2,a,b);
	
	cout << "sol : " << endl;
	for (j=0; j<2; j++)
	{
		for (i=0; i<4; i++)
			cout << b[j*4+i] << " ";
		cout << endl;
	}
	
	cout << "Testing InvertComplexMatrix to invert CMtemp" << endl;
	dComplexMatrix	inverse(4,4);
	InvertComplexMatrix(CMtemp, inverse);
	cout << "inverse: " << endl;
	cout << inverse << endl;
	
	
	return(0); 
}

