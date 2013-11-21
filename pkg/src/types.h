/*
 *
 * Definition of types used in the code.
 *
 */

#ifndef TYPES_H_
#define TYPES_H_

#include <set>
#include <vector>
#include <map>
#include <string>

#include <RcppArmadillo.h>


// Integer types ***************************************

// general signed integer
typedef int Int;
// general unsigned integer
typedef unsigned int PosInt;
// general large unsigned integer
typedef unsigned long long int PosLargeInt;



// STL Set types ***************************************

// a set of ints
typedef std::set<Int> IntSet;
// a set of unsigned ints
typedef std::set<PosInt> PosIntSet;


// STL Vector types ************************************

// (1) Double vectors:

// here is a double vector
typedef std::vector<double> MyDoubleVector;
// and a long double vector
typedef std::vector<long double> LongDoubleVector;

// here is a double deque
typedef std::deque<double> DoubleDeque;


// (2) Integer vectors:

// here is an integer vector
typedef std::vector<Int> IntVector;
// and an unsigned int vector
typedef std::vector<PosInt> PosIntVector;

// (3) Logical vector
typedef std::vector<bool> BoolVector;


// Special types ***************************************

// one power set
typedef std::multiset<Int> Powers;
// vector of powers
typedef std::vector<Powers> PowersVector;


// the "index" type
typedef LongDoubleVector::size_type Ind;
// one index set
typedef std::set<Ind> IndSet;

// Armadillo types ************************************

// a vector
typedef arma::colvec AVector;
// a matrix
typedef arma::mat AMatrix;


// Constants ************************************

// the machine precision
static const double EPS = sqrt(DOUBLE_EPS);



#endif /* TYPES_H_ */
