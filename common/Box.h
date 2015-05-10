/**
@file
@brief Specifies class for box-shaped particles.
*/

#ifndef BOXH
#define BOXH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <list>
#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include "Star.h"
#include "ExtendedVector.h"
#include "Sphere.h"
#include "global.h"

namespace tsquare
{

/**
@class
@brief Class of box-shaped particles.

Box shapes are star shapes, so the Box class derives from the Star class.
*/
class Box : public tsquare::Star {

protected:
    
std::vector<Real> halfdim_;    /**< Vector of box's half-dimensions */

public:
    
/**
@brief Default constructor initializes all members to default values.
*/
Box();

/**
@brief Overloaded constructor sets up the box dimensions.

@param thickness is the box thickness
@param width is the box width
@param length is the box length
@param verbose is true iff extra output is desired
*/
Box(Real thickness, Real width, Real length, bool verbose=false);
    
/**
@name Dimensions
@brief Functions for setting and getting the box dimensions <i>L</i>, <i>W</i>,
and <i>T</i>.
*/

/**@{*/
/**
@brief Sets the box dimensions.

After setting the dimensions, the box dimensions are reordered if necessary so that
\f$ \mbox{L} \ge \mbox{W} \ge \mbox{T}\f$.

@param thickness is the new box thickness
@param width is the box width
@param length is the box length
*/
void setDimensions(Real thickness, Real width, Real length);
    
/**
@brief Sets the box dimensions.

After setting the dimensions, the box dimensions are reordered if necessary so that
\f$ \mbox{L} \ge \mbox{W} \ge \mbox{T}\f$.

@param dims is vector of dimensions
*/
void setDimensions(std::vector<Real> dims);
/**@}*/

/**
@name Radius
@brief \f$R(\theta,\phi)\f$ specifies the distance from the object's
center to a point on the surface with angular coordinates \f$(\theta,\phi)\f$.

This group provides functions for finding <i>R</i> and its derivates with
respect to the polar and azimuthal angles.
*/

/**@{*/
/**
@brief Get the position at a given angular coordinate pair \f$(\theta,\phi)\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the length of the position vector \f$R\f$
*/
Real getR(const Real theta, const Real phi);

/**
@brief Gets the distances from the center to each surface point.
@return the vector storing the distances from the center to each point
         in the surface mesh, relative to the center
*/
virtual std::vector<Real> getR(void) const { return r_; }

/**
@brief Get the derivative of position with respect to azimuthal angle, \f$R_{\phi}\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the derivative \f$R_{\phi}\f$
*/
Real getRp(Real theta, const Real phi);

/**
@brief Get the second derivative of position with respect to azimuthal angle, \f$R_{\phi \phi}\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the second derivative \f$R_{\phi \phi}\f$
*/
Real getRpp(Real theta, const Real phi);

/**
@brief Get the derivative of position with respect to polar angle, \f$R_{\theta}\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the derivative \f$R_{\theta}\f$
*/
Real getRt(Real theta, const Real phi);

/**
@brief Get the second derivative of position with respect to polar angle, \f$R_{\theta \theta}\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the derivative \f$R_{\theta \theta}\f$
*/
Real getRtt(Real theta, const Real phi);

/**
@brief Get the mixed second derivative of position, \f$R_{\theta \phi} = R_{\phi \theta}\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the length of the position vector \f$R_{\theta \phi}\f$
*/
Real getRtp(Real theta, const Real phi);
/**@}*/

/**
@brief Dilate the box by a given linear scale factor

@note The scale factor must be positive definite.

@param dr is the linear scale factor (dimensionless)
*/
void doDilate(const Real dr);
    
/**
@brief Compute (but do not return) the box volume.
*/
void computeVolume(void);

/**
@brief Compute (but do not return) the box surface area.
*/
void computeArea(void);

/**
@brief Compute the box's bounding spheres.

The "bounding spheres" are two spheres:
    - the minimum enclosing sphere, and
    - the maximum inscribed sphere
*/
void computeBoundingspheres(void);

/**
@brief Compute the length of the particle.

The length, <i>L</i>, is the length of the longest line segement that
can be drawn entirely within the object. 

@param sc is the vector of \f$(\theta,\phi)\f$ pairs to consider
@param xc is the <i>x</i> coordinate of each point in sc
@param yc is the <i>y</i> coordinate of each point in sc
@param zc is the <i>z</i> coordinate of each point in sc
@param &lengthvec holds the vector of distance from the center to each point
@param &t1 is the estimated lower bound on the polar angle, \f$\theta_1\f$ (radians)
@param &t2 is the estimated upper bound on the polar angle, \f$\theta_2\f$ (radians)
@param &p1 is the estimated lower bound on the azimuthal angle, \f$\phi_1\f$ (radians)
@param &p2 is the estimated upper bound on the azimuthal angle, \f$\phi_2\f$ (radians)
@return the length, <i>L</i>
*/
Real computeLength(std::vector<std::vector<Real> > sc, std::vector<Real> xc, std::vector<Real> yc,
                   std::vector<Real> zc, std::vector<Real> &lengthvec, Real &t1, Real &t2,
                   Real &p1, Real &p2)
{
    return (dim_[0]);
}

/**
@brief Compute the width of the particle.

The width, <i>W</i>, is the length of the longest line segement that
can be drawn entirely within the object _and_ perpendicular to the
length <i>L</i>.

@param sc is the vector of \f$(\theta,\phi)\f$ pairs to consider in the coarse set
@param xc is the <i>x</i> coordinate of each point in sc (the coarse set)
@param yc is the <i>y</i> coordinate of each point in sc (the coarse set)
@param zc is the <i>z</i> coordinate of each point in sc (the coarse set)
@param lvec is the length vector
@param &wvec will hold the width vector
@param &t1 is the estimated lower bound on the polar angle, \f$\theta_1\f$ (radians)
@param &t2 is the estimated upper bound on the polar angle, \f$\theta_2\f$ (radians)
@param &pp1 is the estimated lower bound on the azimuthal angle, \f$\phi_1\f$ (radians)
@param &pp2 is the estimated upper bound on the azimuthal angle, \f$\phi_2\f$ (radians)
@return the width, <i>W</i>
*/
Real computeWidth(std::vector<std::vector<Real> > sc, std::vector<Real> xc,
                          std::vector<Real> yc, std::vector<Real> zc, std::vector<Real> lvec,
                          std::vector<Real> &wvec, Real &t1, Real &t2, Real &pp1, Real &pp2)
{
    return (dim_[1]);
}

/**
@brief Compute the thickness of the particle.

The thickness, <i>T</i>, is the length of the longest line segement that
can be drawn entirely within the object _and_ perpendicular to both the
length <i>L</i> and width <i>W</i>.

@param sc is the vector of \f$(\theta,\phi)\f$ pairs to consider in the coarse set
@param xc is the <i>x</i> coordinate of each point in sc (the coarse set)
@param yc is the <i>y</i> coordinate of each point in sc (the coarse set)
@param zc is the <i>z</i> coordinate of each point in sc (the coarse set)
@param rc is the set of distances to each point in sc (the coarse set)
@param lvec is the length vector already found
@param wvec is the width vector already found
@param &tvec will hold the thickness vector
@param &t1 is the estimated lower bound on the polar angle, \f$\theta_1\f$ (radians)
@param &t2 is the estimated upper bound on the polar angle, \f$\theta_2\f$ (radians)
@param &pp1 is the estimated lower bound on the azimuthal angle, \f$\phi_1\f$ (radians)
@param &pp2 is the estimated upper bound on the azimuthal angle, \f$\phi_2\f$ (radians)
@return the thickness, <i>T</i>
*/
Real computeThickness(std::vector<std::vector<Real> > sc, std::vector<Real> xc,
                          std::vector<Real> yc, std::vector<Real> zc, std::vector<Real> rc,
                          std::vector<Real> lvec, std::vector<Real> wvec, std::vector<Real> &tvec,
                          Real &t1, Real &t2, Real &pp1, Real &pp2)
{
    return (dim_[2]);
}

/**
@brief Compute the triaxial width of the particle.

The triaxial width, <i>D<sub>2</sub></i>, is the length of the longest
line segment that can be drawn in the plane of the object's maximum cross-sectional
area _and_ perpendicular to the length <i>L</i>.

@param numpoints is the number of quadrature points for projected area 
@warning numpoints must be even
@param rstep is the rotation step size (radians)
@param &rotate_angle will hold the rotation angle (radians)
@param &mpa will hold the maximum projected area
@param &ampa will hold the angle rotated from the beginning of rotations at
    the maximum projected area (radians)
@param &amaxwidth will hold the absolute maximum width found, corresponding
    to <i>W</i>
@param &aamw will hold the angle rotated from the beginning of rotations where
    <i>W</i> was found
@return the triaxial width, <i>D<sub>2</sub></i>
*/
Real computeTriaxialwidth(const int numpoints, const Real rstep, Real &rotate_angle,
                          Real &mpa, Real &ampa, Real &amaxwidth, Real &aamw) 
{
    return (dim_[1]);
}

/**
@brief Compute the triaxial thickness of the particle.

The triaxial thickness, <i>D<sub>3</sub></i>, is the length of the
longest line segment that can be drawn inside the object while simultaneously
being perpendicular to both <i>L</i> _and_ <i>D<sub>2</sub></i>.

@param numpoints is the number of quadrature points for projected area 
@warning numpoints must be even
@param init_angle is the starting angle (radians)
@param d2_angle is the angle at which the triaxial width was found
@param w_angle is the angle at which the width was found
@param &absmax_thickness will hold the absolute maximum thickness
@return the triaxial thickness, <i>D<sub>3</sub></i>
*/
Real computeTriaxialthickness(const int numpoints, const Real init_angle,
                              const Real d2_angle, const Real w_angle,
                              Real &absmax_thickness)
{
    return (dim_[2]);
}

/**
@brief Compute (but do not return) the box dimensions.

@param triaxialcalc is true iff we want triaxial dimensions instead of (<i>L</i>,<i>W</i>,<i>T</i>)
*/
void computeDimensions(bool triaxialcalc=false);

/**
@brief Get the moment of inertia tensor

@return the moment of inertia tensor, <i>I</i>
*/
std::vector<Real> getI(void);
    
/**
@brief Writes formatted box properties to standard out.
*/
void printProperties(void);

/**
@brief Writes formatted box properties to standard out.

@return reference to string describing the class type
*/
std::string &getType(void) const { return (std::string &)("ORTHORHOMBIC BOX"); }

/**
@brief Box class assignment operator.

@param &that is the box object on the left hand side
@return a box object equal to that
*/
Box& operator=(const Box& that);
};

} // End of tsquare namespace
#endif
