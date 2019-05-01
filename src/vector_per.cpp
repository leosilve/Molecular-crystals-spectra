/* Written by Leonardo Silvestri 2007-2013 */

#include <iostream>
#include <vector>
#include <vector_per.h>
#include <OPutils.h>

namespace vector_per_mod
{

int posmod(int a, int b)
{
	if (a%b==0)
		return(0);
	if (a>0)
		return (a-(a/b)*b);
	else
		return (a-(a/b-1)*b);
}

int lmod(int a, int b)
{
	if (b<=0) 
		nrerror("lmod: second argument must be positive !");
	if (a==0)
		return(0);
	if (b%2==1) 
	{
// Second argument is odd
	if (abs(a)<b/2.0)
		return(a);
	if (a>0)
		if (a-(a/b)*b>b/2.0)
			return(a-(a/b)*b-b);
		else
			return(a-(a/b)*b);
	else
		if (a-(a/b-1)*b>b/2.0)
			return(a-(a/b-1)*b-b);
		else
			return(a-(a/b-1)*b);
	}
	else // if (b%2==1)
	{
// Second argument is even
		if (abs(a)<b/2 || a==b/2)
			return(a);
		if (a>0)
			if (a-(a/b)*b>b/2)
				return(a-(a/b)*b-b);
			else
				return(a-(a/b)*b);
		else
			if (a-(a/b-1)*b>b/2)
				return(a-(a/b-1)*b-b);
			else
				return(a-(a/b-1)*b);
	}
}

} // end of namespace vector_per_mod

using namespace std;
using namespace vector_per_mod;

	vector_per vector_per::operator+ (vector_per param) {
		vector_per temp;
		temp[0] = lmod( (n1() + param[0]), N1);
		temp[1] = lmod( (n2() + param[1]), N2);
		temp[2] = lmod( (n3() + param[2]), N3);
	return (temp);
	};
	
	vector_per vector_per::operator- (vector_per param) {
		vector_per temp;
		temp[0] = lmod((n1() - param[0]), N1);
		temp[1] = lmod((n2() - param[1]), N2);
		temp[2] = lmod((n3() - param[2]), N3);
		return (temp);
	};
	
	vector_per vector_per::operator* (int param) {
		vector_per temp;
		temp[0] = lmod((n1()*param), N1);
	temp[1] = lmod((n2()*param), N2);
	temp[2] = lmod((n3()*param), N3);
	return (temp);
	};
	
	void vector_per::operator+= (vector_per param) {
		set_n1( n1() + param[0] );
		set_n2( n2() + param[1] );
	set_n3( n3() + param[2] );
		return;
	};
	
	void vector_per::operator-= (vector_per param) {
		set_n1( n1() - param[0] );
		set_n2( n2() - param[1] );
		set_n3( n3() - param[2] );
		return;
	};
	
	void vector_per::operator*= (int param) {
		set_n1( n1() * param );
		set_n2( n2() * param );
	set_n3( n3() * param );
	return;
	};
	
	int vector_per::operator== (vector_per param) {
		if (n1()==param.n1() && n2()==param.n2() && n3()==param.n3())  
			return(true);
	else
			return(false);
	};
	
	int vector_per::operator!= (vector_per param) {
		if (n1()==param.n1() && n2()==param.n2() && n3()==param.n3())  
			return(false);
		else
			return(true);
	};
	
	int vector_per::operator> (vector_per param) {
		if ( n1()*n1()+n2()*n2()+n3()*n3() > param.n1()*param.n1()+param.n2()*param.n2()+param.n3()*param.n3() ) 
			return(true);
		else
		{
			if ( n1()*n1()+n2()*n2()+n3()*n3() == param.n1()*param.n1()+param.n2()*param.n2()+param.n3()*param.n3() ) 
			{
			if (n1() > param.n1())
					return(true);
				else
					if (n1() == param.n1())
					{
						if (n2() > param.n2())
							return(true);
						else
							if (n2() == param.n2())
							{
								if (n3() > param.n3())
									return(true);
							}
					}
			}
		}			
		return(false);
	};
	
	int vector_per::operator>= (vector_per param) {		
		if ( *this > param || *this == param) 
			return(true);
		else
			return(false);
	};
		
	int vector_per::operator< (vector_per param) {
		if ( *this >= param )  
			return(false);
		else
			return(true);
	};
	
	int vector_per::operator<= (vector_per param) {
		if ( *this > param )  
			return(false);
		else
			return(true);
	};

// Declaration of static members
int vector_per::N1=1;
int vector_per::N2=1;
int vector_per::N3=1;

// Friend functions
int N1() { return(vector_per::N1); };
int N2() { return(vector_per::N2); };
int N3() { return(vector_per::N3); };
void set_vector_per_N(int _N1, int _N2, int _N3)
{
	vector_per::N1=_N1;
	vector_per::N2=_N2;
	vector_per::N3=_N3;
};


ostream &operator<<(ostream& stream, vector_per v)   // Screen output of a vector_per
{
	stream << "(" << v[0] << "," << v[1] << "," << v[2] << ")"; 
	return stream;
};
