/**
@file
@brief Method definitions for ExtendedVector class.
*/
#include "ExtendedVector.h"

tsquare::ExtendedVector::ExtendedVector()
{
    lo_ = 0;
    hi_ = 0;
    vec_.clear();
}

tsquare::ExtendedVector::ExtendedVector(const int lo, const int hi)
{
    lo_ = lo;
    hi_ = hi;
    vec_.clear();
    vec_.resize(hi_ - lo_ + 1);
}

void tsquare::ExtendedVector::resize(const int lo, const int hi)
{
    lo_ = lo;
    hi_ = hi;
    vec_.clear();
    vec_.resize(hi_ - lo_ + 1);
}

void tsquare::ExtendedVector::clear()
{
    lo_ = 0;
    hi_ = 0;
    vec_.clear();
}

tsquare::ExtendedVector &tsquare::ExtendedVector::operator=(const tsquare::ExtendedVector &rhs)
{
    if (this == &rhs) {
        return *this;
    }
    lo_ = rhs.lo_;
    hi_ = rhs.hi_;
    vec_.clear();
    for (int i = 0; i < rhs.vec_.size(); i++) {
        vec_.push_back(rhs.vec_[i]);
    }
    return *this;
}

std::complex<Real> &tsquare::ExtendedVector::operator[](int i)
{
    if (i < lo_ || i > hi_) {
        std::cout << "Matrix EOB Error: i = " << i << ", lo_ = "
             << lo_ << ", hi_ = " << hi_ << std::endl;
        return vec_[0];
    }
    return vec_[i - lo_];
}
