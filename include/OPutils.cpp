/* Written by Leonardo Silvestri 2007-2013 */

#include <fstream>
#include <iostream>
#include <istream>
#include <vector>
#include <complex>
#include <string> 
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <ME.h>
//#include <mpstate.h>

/* Numerical Recipes standard error handler */
//extern  void nrerror(char error_text[]);

using namespace std;

/* Numerical Recipes standard error handler */
void nrerror(string error_text)
{
	cout << "Numerical Recipes run-time error..." << endl;
	cout << error_text << endl;
	cout << "...now exiting to system..." << endl;
	getchar();
	exit(1);
}


// Screen output of a std vector
template<class T>
ostream& operator<<(ostream& stream, const vector<T>& v)   
{
	stream << "(";
	for (unsigned int i=0; i< v.size()-1; i++)
		stream << v[i] << "," ;
	stream << v[v.size()-1] << ")"; 
	return stream;
};

/*
// Screen output of a ublas::vector 
template<class T>
ostream &operator<<(ostream& stream, boost::numeric::ublas::vector<T>& v)   
{
	stream << "(";
	for (unsigned int i=0; i< v.size()-1; i++)
		stream << v(i) << "," ;
	stream << v(v.size()-1) << ")"; 
	return stream;
}*/

/*
// Screen output of a ublas::matrix 
template<class T>
ostream &operator<<(ostream& stream, boost::numeric::ublas::matrix<T>& m)   
{
	stream << "[" << m.size1() << "," << m.size2() << "]"; 
	stream << "(";
	for (unsigned int i=0; i< m.size1()-1; i++)
	{
		stream << "(";
		for (unsigned int j=0; j< m.size2()-1; j++)
			stream << m(i,j) << "," ;
		stream << m(i,m.size2()-1) << ")";
		stream << ",";
	}
	stream << "(";
	for (unsigned int j=0; j< m.size2()-1; j++)
		stream << m(m.size1()-1,j) << "," ;
	stream << m(m.size1()-1,m.size2()-1) << "))";
	return stream;
}*/


template<typename T>
void get_next_from_file(std::ifstream& infile, T& data)
{
	string line;
	getline(infile, line, ':');
	infile >> data;
	return;
}

void get_next_string_from_file(ifstream& infile, string& mystring)
{
	getline(infile, mystring, ':');
	getline(infile, mystring, '"');
	getline(infile, mystring, '"');	
	return;
}

// Explicit instantiation of templates

template ostream& operator<<(ostream&, const vector<complex<double> >&);
template ostream& operator<<(ostream&, const vector<complex<float> >&);
template ostream& operator<<(ostream&, const vector<double>&);
template ostream& operator<<(ostream&, const vector<float>&);
template ostream& operator<<(ostream&, const vector<int>&);
template ostream& operator<<(ostream&, const vector<string>&);
template ostream& operator<<(ostream&, const vector<vector<int> >&);
template ostream& operator<<(ostream&, const vector<vector<float> >&);
template ostream& operator<<(ostream&, const vector<vector<double> >&);
template ostream& operator<<(ostream&, const vector<vector<complex<double> > >&);
template ostream& operator<<(ostream&, const vector<vector<complex<float> > >&);
//template ostream& operator<<(ostream&, const vector<ME<float> >&);
//template ostream& operator<<(ostream&, const vector<MPSTATE<float> >&);
//template ostream& operator<<(ostream&, const vector<ME<double> >&);
//template ostream& operator<<(ostream&, const vector<MPSTATE<double> >&);

/*
template ostream& operator<<(ostream&, boost::numeric::ublas::vector<double>&);
template ostream& operator<<(ostream&, boost::numeric::ublas::vector<float>&);
template ostream& operator<<(ostream&, boost::numeric::ublas::vector<int>&);
template ostream& operator<<(ostream&, boost::numeric::ublas::vector<complex<double> >&);
template ostream& operator<<(ostream&, boost::numeric::ublas::vector<complex<float> >&);
*/
/*
template ostream& operator<<(ostream&, boost::numeric::ublas::matrix<double>&);
template ostream& operator<<(ostream&, boost::numeric::ublas::matrix<float>&);
template ostream& operator<<(ostream&, boost::numeric::ublas::matrix<int>&);
template ostream& operator<<(ostream&, boost::numeric::ublas::matrix<complex<double> >&);
template ostream& operator<<(ostream&, boost::numeric::ublas::matrix<complex<float> >&);
*/

template void get_next_from_file(std::ifstream& infile, unsigned int& data);
template void get_next_from_file(std::ifstream& infile, int& data);
template void get_next_from_file(std::ifstream& infile, float& data);
template void get_next_from_file(std::ifstream& infile, double& data);
