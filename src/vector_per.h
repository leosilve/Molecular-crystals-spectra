/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _VECTOR_PER_H
#define _VECTOR_PER_H

#include <iostream>
#include <vector>

namespace vector_per_mod
{
	int posmod(int a, int b);
	int lmod(int a, int b);

} // end of namespace vector_per_mod


//using namespace std;
//using namespace vector_per_mod;

///////////////////////////////////////////////////////////////////////////
///    Periodic vector 3D class (inheritance from vector class).
///
///    It defines a 3D vector in terms of its integer coefficients
///    n1, n2 and n3 over an orthogonal basis. The period is N1, N2 and N3, 
///    respectively in the three directions. It is used for the 
///	   MultiPhonon basis set.
///////////////////////////////////////////////////////////////////////////

class vector_per : public std::vector<int>
{
	
private:
	static int N1;
	static int N2;
	static int N3;

public:
	template<class> friend class MultiPhonon;
	friend void set_vector_per_N(int _N1, int _N2, int _N3);
	friend int N1() ;
	friend int N2() ;
	friend int N3() ;

	vector_per()	: std::vector<int>(3,0) {}                               // default constructor inherited from class
	vector_per(int _n1, int _n2, int _n3) : std::vector<int>(3) 
	{(*this)[0] =vector_per_mod::lmod(_n1,N1); (*this)[1]=vector_per_mod::lmod(_n2,N2); (*this)[2]=vector_per_mod::lmod(_n3,N3);};		 // alternative constructor
	void set(int _n1, int _n2, int _n3) 
	{(*this)[0] =vector_per_mod::lmod(_n1,N1); (*this)[1]=vector_per_mod::lmod(_n2,N2); (*this)[2]=vector_per_mod::lmod(_n3,N3);};
	void set_3d(vector_per  _v) 
		{(*this)[0] =_v.n1(); (*this)[1]=_v.n2(); (*this)[2]=_v.n3();};
	void set_n1(int _val) 	{(*this)[0] =vector_per_mod::lmod(_val,N1); };
	void set_n2(int _val) 	{(*this)[1] =vector_per_mod::lmod(_val,N2); };
	void set_n3(int _val) 	{(*this)[2] =vector_per_mod::lmod(_val,N3); };
	int n1( ) {	return((*this)[0]); };
	int n2( ) {	return((*this)[1]); };
	int n3( ) {	return((*this)[2]); };

	vector_per operator+ (vector_per param);	
	vector_per operator- (vector_per param);
	vector_per operator* (int param);
	void operator+= (vector_per param);
	void operator-= (vector_per param);
	void operator*= (int param);
	int operator== (vector_per param);
	int operator!= (vector_per param);
	int operator> (vector_per param);
	int operator>= (vector_per param);
	int operator< (vector_per param);
	int operator<= (vector_per param);

	void resize() {std::cout << "Vector_per cannot be resized!"<< std::endl; return;}

};  // end of vector_per class

std::ostream &operator<<(std::ostream& stream, vector_per v);   // Screen output of a vector_per

#endif