/**
@file
@brief Method definitions for the Sphere class.
*/
#include "Sphere.h"

tsquare::Sphere& tsquare::Sphere::operator=(const tsquare::Sphere &rhs)
{
    // Only do the assignment if RHS is a different object

    if (this != &rhs) {
        radius_ = rhs.getRadius();
        xr_ = rhs.getXr();
        yr_ = rhs.getYr();
        zr_ = rhs.getZr();
        volume_ = rhs.getVolume();
        area_ = rhs.getArea();
    }
    return *this;
}

bool tsquare::Sphere::operator==(const tsquare::Sphere &other)
{
    bool answer = false;
    if ((radius_ == other.getRadius())
        && (xr_ == other.getXr())
        && (yr_ == other.getYr())
        && (zr_ == other.getZr()) ) answer = true;
    return answer;
}

bool tsquare::Sphere::operator!=(const tsquare::Sphere &other)
{
    return !(*this == other);
}

bool tsquare::Sphere::operator< (const tsquare::Sphere &other)
{
    return (radius_ < other.getRadius());
}

bool tsquare::Sphere::operator> (const tsquare::Sphere &other)
{
    return (radius_ > other.getRadius());
}

bool tsquare::Sphere::operator<= (const tsquare::Sphere &other)
{
    return (radius_ <= other.getRadius());
}

bool tsquare::Sphere::operator>= (const tsquare::Sphere &other)
{
    return (radius_ >= other.getRadius());
}
