/* Written by Leonardo Silvestri 2007-2013 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <pos.h>

namespace ublas = boost::numeric::ublas;
using namespace std;

///////////////////////////////////////////////////////////////////////////
/////////	struct pos   (position of a molecule inside a crystal) ////////
///////////////////////////////////////////////////////////////////////////

template<class FLOAT>
inline bool operator<(const pos<FLOAT>& p1, const pos<FLOAT>& p2)
{
	if (norm_2(p1.n)<norm_2(p2.n))
		return(true);
	else
		return(false);
}

template <typename FLOAT>
ostream &operator<<(ostream& stream, const pos<FLOAT> pos)   // Screen output of a molecular position 
{
	stream << "(" << pos.n << "," << pos.alpha << ")" ;
	return stream;
}

template inline bool operator<(const pos<double>& p1, const pos<double>& p2);
template inline bool operator<(const pos<float>& p1, const pos<float>& p2);
template ostream &operator<<(ostream& stream, const pos<double> pos);
template ostream &operator<<(ostream& stream, const pos<float> pos);
