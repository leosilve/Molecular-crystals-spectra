/* Written by Leonardo Silvestri 2007-2013 */

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

//////////////////////////////////////////////////////////
// Reads input file supplied by user

	if (argc<2)
		nrerror("Please supply input file name!");
	
	string model_input_file=argv[1];
	
	OPmodel<MYFLOAT, PhononCloud<MYFLOAT> >*	MyModel = new OPmodel<MYFLOAT, PhononCloud<MYFLOAT> >(model_input_file);
	
	MyModel->solve();
	
	delete MyModel;
	
	return(0);
}

