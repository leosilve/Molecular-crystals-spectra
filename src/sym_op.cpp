/* Written by Leonardo Silvestri 2007-2013 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <fstream>
#include <iostream>
//#include <string>
//#include <math.h>
#include <vector>
#include <sym_op.h>
#include <OPutils.h>

namespace ublas = boost::numeric::ublas;
using namespace std;


///////////////////////////////////////////////////////////////////////////
/////////     class sym_op    (symmetry operation)				  /////////
///////////////////////////////////////////////////////////////////////////

template <typename FLOAT>
void sym_op<FLOAT>::apply(pos<FLOAT> &ini, pos<FLOAT> &fin)
	{
		fin.n=ini.n-ax;
		fin.n=prod(rot,fin.n);
		fin.n+=ax+tr;
		fin.alpha=map[ini.alpha];
		return;
	};

template <typename FLOAT>
void sym_op<FLOAT>::apply(pos<FLOAT> &p)
	{
		p.n=p.n-ax;
		p.n=prod(rot,p.n);
		p.n+=ax+tr;
		p.alpha=map[p.alpha];
		return;
	};

// Long constructor
template <typename FLOAT>
sym_op<FLOAT>::sym_op(typename sym_op<FLOAT>::Matrix &_rot, typename sym_op<FLOAT>::Vector &_ax, 
	   typename sym_op<FLOAT>::Vector &_tr, vector<int> &_map)
	{
		if ( _rot.size1()!=3 || _rot.size2()!=3 )
			nrerror("sym_op: cannot set rot!");
		if ( _ax.size()!=3 || _tr.size()!=3 )
			nrerror("sym_op: cannot set ax or tr!");
		set_rot(_rot);
		set_ax(_ax);
		set_tr(_tr);
		set_map(_map);
	};	

template <typename FLOAT>
ostream &operator<<(ostream& stream, const sym_op<FLOAT>& op)   // Screen output of a symmetry operation
{
	stream << "Rotational part: "		<< endl << op.get_rot() << endl; 
	stream << "Axis position: "			<< endl << op.get_ax() << endl;
	stream << "Traslational part : "	<< endl << op.get_tr() << endl;
	stream << "Molecular mapping : "	<< endl; // << op.get_map() << endl;
	for (int i=0; i<op.get_map_size(); i++)
		stream << op.get_map(i) << " ";
	stream << endl;
	return stream;
};

////////////////////////////////////////
// Explicit instantiation

template class sym_op<double>;
template class sym_op<float>;
template ostream &operator<<(ostream&, const sym_op<double>& op);
template ostream &operator<<(ostream&, const sym_op<float>& op);


