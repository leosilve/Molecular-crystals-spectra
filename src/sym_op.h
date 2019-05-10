/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _SYM_OP_H
#define _SYM_OP_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <pos.h>

namespace ublas = boost::numeric::ublas;

///////////////////////////////////////////////////////////////////////////
///     Symmetry operations.
///
///     Symmetry operation of the molecular crystal defined by
///     - Rotation matrix
///     - Rotation axis position (in fractional coordinates)
///     - Translation (in fractional coordinates)
///     - Molecule mapping
///     This class is used by the "lattice" class and in various calculations
///     to take advantage of symmetry operations and reduce the size of the matrices to be diagonalised
///////////////////////////////////////////////////////////////////////////
template <class FLOAT>
class sym_op {

	typedef ublas::matrix<FLOAT, ublas::column_major>	Matrix;
	typedef ublas::vector<FLOAT>						Vector;

private:
	Matrix		rot; ///< Rotation matrix 3x3
	Vector		ax;  ///< Rotation axis position (in fractional coordinates)
	Vector		tr;  ///< Translation (in fractional coordinates)
	std::vector<int>	map; ///< Molecule mapping. This is a vector of integers
                             ///< where the value of the i-th component indicates the non equivalent
                             ///< molecular position inside the unit cell into which molecule i-th is transformed by the symmetry operation

public:

	void set_rot(Matrix _rot)	{rot=_rot; };
	void set_ax(Vector _ax)		{ax=_ax; };
	void set_tr(Vector _tr)		{tr=_tr; };
	void set_map(int i, int n)	{map[i]=n; };
	void set_map(std::vector<int> _map) {map=_map; };
	void set_map_size(int n) {map.resize(n); };
	
	std::vector<int>		get_map()		const { return(map); };
	int				get_map(int i)	const { return(map[i]); };
	int				get_map_size()	const { return(map.size()); };
	Vector			get_ax()		const { return(ax); };
	Vector			get_tr()		const { return(tr); };
	Matrix			get_rot()		const { return(rot); };
	
	void apply(pos<FLOAT> &ini, pos<FLOAT> &fin);
	void apply(pos<FLOAT> &p);

	sym_op()
	{
		Matrix		rot(3,3);
		Vector		ax(3);
		Vector		tr(3);
		std::vector<int>	map(0);
		for (unsigned i=0; i < rot.size1(); ++i)
		{
			ax(i)=0.0;
			tr(i)=0.0;
			for (unsigned j=0; j<rot.size2(); ++j)
				rot(i,j)=0.0;
		}
	};	// Short constructor
	
	sym_op(Matrix &_rot, Vector &_ax, Vector &_tr, std::vector<int> &_map); // Long constructor
};


// Screen output of a symmetry operation
template <typename MYFLOAT>
std::ostream &operator<<(std::ostream& stream, const sym_op<MYFLOAT>& op);   

#endif
