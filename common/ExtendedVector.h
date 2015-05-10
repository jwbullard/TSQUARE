/**
@file
@brief Specifies the ExtendedVector class.
*/
#ifndef EXTENDEDVECTORH
#define EXTENDEDVECTORH

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "global.h"

namespace tsquare
{

/**
@class ExtendedVector
@brief Extension of the STL vector container.

The ExtendedVector class generalizes the C++ STL vector container for
holding complex numbers.  It has all the properties of the STL vector,
but can have any lower and upper indices.  This is specifically
designed for holding spherical harmonic coefficients.

@todo Consider making this a template class so that the extended vector
can hold any kind of object, not just complex numbers.
*/
class ExtendedVector {

private:

/**
@brief The vector container on which ExtendedVector is based.
*/
std::vector<std::complex<Real> > vec_;

int lo_;    /**< The first element index. */
int hi_;    /**< The last element index. */

public:

/**
@brief Default constructor.

Initializes the low and high index to zero and clears the vector.
*/
ExtendedVector();

/**
@brief Overloaded constructor.

Sets the first and last element index, and then resizes the
vector as needed.

@param lo is the value of the first index
@param hi is the value of the last index
*/
ExtendedVector(const int lo, const int hi);

/**
@brief Destructor.
*/
~ExtendedVector() { vec_.clear(); lo_ = 0; hi_ = 0; }

/**
@brief Gets the size of the vector.

@return the number of elements in the vector
*/
int size() const { return (int)(vec_.size()); }

/**
@brief Gets the first element's index.

@return the first element's index
*/
int getLo() const { return lo_; }

/**
@brief Gets the last element's index.

@return the last element's index
*/
int getHi() const { return hi_; }

/**
@brief Resizes the vector with new first and last indices.

@param lo is the new first index
@param hi is the new last index
*/
void resize(const int lo, const int hi);

/**
@brief Clears the vector to have zero elements.
*/
void clear();

/**
@brief Sets an element of the vector.

@param i is the index of the <i>i</i>-th element from the first index
@param val is the value to place in the element
*/
void setElement(const int i, const Real val) { vec_[i - lo_] = val; }

/**
@brief Gets an element of the vector.

@param i is the index of the <i>i</i>-th element from the first index
@return the value stored in the <i>i</i>-th element from the first index
*/
std::complex<Real> getElement(const int i) const { return vec_[i - lo_]; }

/**
@brief Overloaded assignment operator.

@param &rhs is a reference to the right hand side of the assignment
@return the assigned ExtendedVector
*/
ExtendedVector &operator=(const ExtendedVector &rhs);

/**
@brief Overloaded bracket operator.

@param i is the number of elements to go past the first element
@return the value stored in the <i>i</i>-th element past the first one
*/
std::complex<Real> &operator[](int i);

};

} // End of tsquare namespace
#endif
