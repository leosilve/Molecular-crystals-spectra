/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _POS_H
#define _POS_H

#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
#include <iostream>
#include <vector>

namespace ublas = boost::numeric::ublas;

///////////////////////////////////////////////////////////////////////////
/////////	struct pos   (position of a molecule inside a crystal) ////////
///////////////////////////////////////////////////////////////////////////

template <typename FLOAT>
struct pos {
	ublas::vector<FLOAT> n;
	int alpha;
};

template<class FLOAT>
inline bool operator<(const pos<FLOAT>& p1, const pos<FLOAT>& p2);

// Screen output of a molecular position 
template <typename FLOAT>
std::ostream &operator<<(std::ostream& stream, const pos<FLOAT> pos);   

#endif