/* Written by Leonardo Silvestri 2007-2013 */

#include <lattice.h>
//#include <N_Functions.h>
//#include <MultiPhonon.h>
//#include <PhononCloud.h>
//#include <OP_model.h>
//#include <OP_epsilon.h>
//#include <OP_absorption.h>
#include <OPutils.h>
#include <iostream>
#include <vector>
#include <complex>

#ifndef Pi
#define Pi 3.141592654
#endif 
#define MYFLOAT double

namespace ublas = boost::numeric::ublas;
using namespace std;

typedef ublas::matrix<MYFLOAT, ublas::column_major>	Matrix;
typedef ublas::vector<MYFLOAT>	Vector;

int main (int argc, char* argv[]) {

	
//////////////////////////////////////////////////////////
// Constructs lattice and computes excitonic bands

	if (argc<2)
		nrerror("Please supply a lattice input file name!");

	string lattice_file=argv[1];

	lattice<MYFLOAT>	MyLattice(lattice_file);
	
	ostringstream ss;
	ss.str("");
	ss << MyLattice.get_lattice_name() << ".log";
	ofstream log_file(ss.str().c_str());
	MyLattice.print(log_file);
	log_file.close();
	

	if (argc>2)
	{
		string bands_input_file=argv[2];
		MyLattice.computes_lattice_bands(bands_input_file);
	}
		
	return 0;

}

