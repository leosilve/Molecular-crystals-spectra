/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _OPUTILS_H_
#define _OPUTILS_H_

#include <iostream>
#include <vector>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/matrix.hpp>

/* Numerical Recipes standard error handler */
void nrerror(std::string error_text);

// Screen output of a std vector
template<class T>
std::ostream &operator<<(std::ostream& stream, const std::vector<T>& v);   

/*
// Screen output of a ublas::vector 
template<class T>
std::ostream &operator<<(std::ostream& stream, boost::numeric::ublas::vector<T>& v);   
*/
/*
// Screen output of a ublas::matrix 
template<class T>
std::ostream &operator<<(std::ostream& stream, boost::numeric::ublas::matrix<T>& m);   
*/

template<typename T>
void get_next_from_file(std::ifstream& infile, T& data);

void get_next_string_from_file(std::ifstream& infile, std::string& mystring);

#endif /* _OPUTILS_H_ */
