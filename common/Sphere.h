/**
@file
@brief Specifies the Sphere class.

As defined here, a Sphere object is a geometric object with four
degrees of freedom:  radius and three Cartesian coordinate of the center.
*/
#ifndef SPHEREH
#define SPHEREH

#include <vector>
#include <cmath>
#include "global.h"

namespace tsquare
{

/**
@class
@brief The Sphere class.

As defined here, a Sphere object is a geometric object with four
degrees of freedom:  radius and three Cartesian coordinate of the center.
*/
class Sphere {

protected:
    
Real radius_;         /**< sphere radius */
Real xr_;             /**< <i>x</i> coordinate of the sphere center */
Real yr_;             /**< <i>y</i> coordinate of the sphere center */
Real zr_;             /**< <i>z</i> coordinate of the sphere center */
Real volume_;         /**< sphere volume */
Real area_;           /**< sphere surface area */

public:
    
/**
@brief Default constructor.

Sets the sphere radius to 1 and the center to the origin.
*/
Sphere() {
    radius_ = 1.0;
    xr_ = yr_ = zr_ = 0.0;
    setVolume();
    setArea();
}

/**
@brief Overloaded constructor.

Sets the sphere radius to value passed as a parameter, with center at the origin.

@param radius is the prescribed radius of the sphere
*/
Sphere(Real radius) {
    radius_ = radius;
    xr_ = yr_ = zr_ = 0.0;
    setVolume();
    setArea();
}
    
/**
@brief Overloaded constructor.

Sets the sphere radius and center coordinates.

@param radius is the prescribed radius of the sphere
@param xr is the <i>x</i> coordinate of the center
@param yr is the <i>y</i> coordinate of the center
@param zr is the <i>z</i> coordinate of the center
*/
Sphere(Real radius, Real xr, Real yr, Real zr) {
    radius_ = radius;
    xr_ = xr;
    yr_ = yr;
    zr_ = zr;
    setVolume();
    setArea();
}

/**
@brief Sets the sphere radius.

Function ensures that the radius is a positive number, otherwise
forces the radius to zero.

@param radius is the prescribed radius of the sphere
*/
void setRadius(const Real radius) {
    radius_ = (radius >= 0.0) ? radius : 0.0;
    setVolume();
    setArea();
}

/**
@brief Sets the <i>x</i> coordinate of the sphere center.

@param xr is the <i>x</i> coordinate of the center
*/
void setXr(const Real xr) { xr_ = xr; }

/**
@brief Sets the <i>y</i> coordinate of the sphere center.

@param yr is the <i>y</i> coordinate of the center
*/
void setYr(const Real yr) { yr_ = yr; }

/**
@brief Sets the <i>z</i> coordinate of the sphere center.

@param zr is the <i>z</i> coordinate of the center
*/
void setZr(const Real zr) { zr_ = zr; }

/**
@brief Sets the coordinates of the sphere center.

@param xr is the <i>x</i> coordinate of the center
@param yr is the <i>y</i> coordinate of the center
@param zr is the <i>z</i> coordinate of the center
*/
void setPosition(const Real xr, const Real yr, const Real zr) {
    xr_ = xr; yr_ = yr; zr_ = zr;
}

/**
@brief Sets the coordinates of the sphere center.

@param pos is a vector holding the Cartesian coordinates of the center
*/
void setPosition(std::vector<Real> pos) {
    xr_ = pos[0]; yr_ = pos[1]; zr_ = pos[2];
}

/**
@brief Computes and sets the volume of the sphere.

\f[
V = \frac{4}{3} \pi r^3
\f]
*/
void setVolume() { volume_ = (4.0 * Pi / 3.0) * std::pow(radius_,3.0); }

/**
@brief Computes and sets the surface area of the sphere.

\f[
V = 4 \pi r^2
\f]
*/
void setArea() { area_ = (4.0 * Pi) * std::pow(radius_,2.0); }
    
/**
@brief Gets the sphere radius.

@return the sphere radius
*/
Real getRadius(void) const { return radius_; }

/**
@brief Gets the <i>x</i> coordinate of the sphere center.

@return the <i>x</i> coordinate of the sphere center
*/
Real getXr(void) const { return xr_; }

/**
@brief Gets the <i>y</i> coordinate of the sphere center.

@return the <i>y</i> coordinate of the sphere center
*/
Real getYr(void) const { return yr_; }

/**
@brief Gets the <i>z</i> coordinate of the sphere center.

@return the <i>z</i> coordinate of the sphere center
*/
Real getZr(void) const { return zr_; }

/**
@brief Gets the Cartesian coordinates of the sphere center

@return a vector storing the three Cartesian coordinates of the center
    - 0 = <i>x</i> coordinate
    - 1 = <i>y</i> coordinate
    - 2 = <i>z</i> coordinate
*/
std::vector<Real> getPosition() {
    std::vector<Real> pos;
    pos.clear();
    pos.push_back(xr_); pos.push_back(yr_); pos.push_back(zr_);
    return pos;
}

/**
@brief Gets (does not compute) the sphere volume.

@return the sphere volume
*/
Real getVolume() const { return volume_; }

/**
@brief Gets (does not compute) the sphere surface area.

@return the sphere surface area
*/
Real getArea() const { return area_; }

/**
@brief Compute and return the distance of the sphere center form the origin.

@return the distance of the sphere center from the origin
*/
Real getDistance() const { return std::sqrt(xr_*xr_ + yr_*yr_ + zr_*zr_); }

/**
@brief Compute and return the polar angle of the sphere center relative to the origin
       of a spherical polar coordinate system.

\f[
\theta = \text{Arccos}(z_r / r)
\f]
where <i>r</i> is the distance from the origin to the sphere center.

@return the polar angle \f$\theta\f$ (radians)
*/
Real getTheta() const {
    Real r = getDistance();
    if (r >= 0.0) return (std::acos(zr_/r));
    return (0.0);
}

/**
@brief Compute and return the azimuthal angle of the sphere center relative to the origin
       of a spherical polar coordinate system.

\f[
\phi = \pi + \text{Arctan}(y_r / x_r)
\f]

@return the azimuthal angle \f$\phi\f$ (radians)
*/
Real getPhi() const {
    Real r = getDistance();
    if (r >= 0.0) return (Pi + std::atan2(yr_,xr_));
    return (0.0);
}
    
/**
@brief Translate the sphere to a new position.

@param dx is the <i>x</i> displacement of the sphere center
@param dy is the <i>y</i> displacement of the sphere center
@param dz is the <i>z</i> displacement of the sphere center
*/
void doTranslate(const Real dx, const Real dy, const Real dz) {
    xr_ += dx; yr_ += dy; zr_ += dz;
    return;
}

/**
@brief Dilate or contract the sphere

@param dr is the amount byt which to increment or decrement the sphere radius
*/
void doDilate(const Real dr) {
    radius_ += dr;
    return;
}
    
/**
@brief Overloaded assignment operator.

@param &rhs is the sphere object on the right side of the assignment.
@return the assigned sphere object
*/
Sphere & operator=(const Sphere &rhs);

/**
@brief Overloaded boolean equality operator.

@param &other is the sphere object on the right side of the operator.
@return true iff the sphere object on the left side is equal to that on the right.
*/
bool operator==(const Sphere &other);

/**
@brief Overloaded boolean inequality operator.

@param &other is the sphere object on the right side of the operator.
@return true iff the sphere object on the left side is _not_ equal to that on the right.
*/
bool operator!=(const Sphere &other);

/**
@brief Overloaded boolean less-than operator.

@param &other is the sphere object on the right side of the operator.
@return true iff the sphere object on the left side has a smaller radius than that on the right.
*/
bool operator< (const Sphere &other);

/**
@brief Overloaded boolean greater-than operator.

@param &other is the sphere object on the right side of the operator.
@return true iff the sphere object on the left side has a greater radius than that on the right.
*/
bool operator> (const Sphere &other);

/**
@brief Overloaded boolean less-than-or-equal-to operator.

@param &other is the sphere object on the right side of the operator.
@return true iff the sphere object on the left side has a radius less than or equal to
        the radius of the sphere object on the right.
*/
bool operator<= (const Sphere &other);

/**
@brief Overloaded boolean greater-than-or-equal-to operator.

@param &other is the sphere object on the right side of the operator.
@return true iff the sphere object on the left side has a radius greater than or equal to
        the radius of the sphere object on the right.
*/
bool operator>= (const Sphere &other);
};

} // End of tsquare namespace

#endif
