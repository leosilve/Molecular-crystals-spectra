/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _POS_H
#define _POS_H

#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
#include <iostream>
#include <vector>

namespace ublas = boost::numeric::ublas;


///	Position of a molecule inside a crystal.
///
///	A molecule inside a crystal is identified by 
/// a 3D vector indicating the position of its unit cell
/// and by an integer specifying the molecular position 
/// occupied inside the cell. In symbols |n,alpha>

template <typename FLOAT>
struct pos {
	ublas::vector<FLOAT> n; ///< position of the unit cell 
	int alpha;				///< position inside the unit cell
};

template<class FLOAT>
inline bool operator<(const pos<FLOAT>& p1, const pos<FLOAT>& p2);

/// Screen output of a molecular position 
template <typename FLOAT>
std::ostream &operator<<(std::ostream& stream, const pos<FLOAT> pos);   

#endif