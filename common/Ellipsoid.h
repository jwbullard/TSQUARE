/**
@file
@brief Specifies the class for axially symmetric ellipsoids.
*/
#ifndef ELLIPSOIDH
#define ELLIPSOIDH

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
@brief The ellipsoid class.

Ellipsoids are star-shaped, so the ellipsoid class derives from the Star class.
*/
class Ellipsoid : public tsquare::Star {

protected:
    
Real adim_;            /**< half the a-axis of the ellipsoid */
Real cdim_;            /**< half the c-axis of the ellipsoid */
Real aspectratio_;     /**< aspect ratio, c/a */

public:
    
/**
@brief Default constructor initializes data.
*/
Ellipsoid();

/**
@brief Constructor sets the ellipsoid dimenions to specified values.

@param a is the a-axis
@param c is the c-axis
@param verbose is true iff extra output is desired
*/
Ellipsoid(Real a, Real c, bool verbose);
    
/**
@brief Sets the semi-axis in the <i>xy</i>-plane.

@param a is the semi-axis to set.
*/
void setA(const Real a);

/**
@brief Gets the semi-axis in the <i>xy</i>-plane.

@return the semi-axis in the <i>xy</i>-plane.
*/
Real getA() const { return adim_; }

/**
@brief Sets the semi-axis along the <i>z</i> coordinate. 

@param c is the semi-axis to set.
*/
void setC(const Real c);

/**
@brief Gets the semi-axis along the <i>z</i> coordinate.

@return the semi-axis along the <i>z</i> coordinate.
*/
Real getC() const { return cdim_; }

/**
@brief Gets the aspect ratio, <i>c</i>/<i>a</i>.

@return the aspect ratio (dimensionless)
*/
Real getAspectratio() const { return aspectratio_; }

/**
@brief Sets the position at a given angular coordinate pair \f$(\theta,\phi)\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
*/
void setR(const Real theta, const Real phi);

/**
@brief Sets the position at a given surface point index.

@param surfaceindex is the index of the required surface point
*/
void setR(const unsigned int surfaceindex);
    
/**
@brief Gets the position at a given angular coordinate pair \f$(\theta,\phi)\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the length of the position vector in the \f$(\theta,\phi)\f$ direction
*/
Real getR(const Real theta, const Real phi);

/**
@brief Gets the distances from the center to each surface point.
@return the vector storing the distances from the center to each point
         in the surface mesh, relative to the center
*/
virtual std::vector<Real> getR(void) const { return r_; }

/**
@brief Dilate the ellipsoid by a given linear scale factor

@note The scale factor must be positive definite.

@param dr is the linear scale factor (dimensionless)
*/
void doDilate(const Real dr);
    
/**
@brief Get the derivative of position with respect to azimuthal angle, \f$R_{\phi}\f$.

This derivative is identically zero for axisymmetric ellipsoids.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the derivative \f$R_{\phi} \equiv 0\f$
*/
Real getRp(Real theta, const Real phi) const { return 0.0; }

/**
@brief Get the second derivative of position with respect to azimuthal angle, \f$R_{\phi \phi}\f$.

This derivative is identically zero for axisymmetric ellipsoids.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the second derivative \f$R_{\phi \phi} \equiv 0\f$
*/
Real getRpp(Real theta, const Real phi) const { return 0.0; }

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

This derivative is identically zero for axisymmetric ellipsoids.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the length of the position vector \f$R_{\theta \phi} \equiv 0\f$
*/
Real getRtp(Real theta, const Real phi) const { return 0.0; }
    
/**
@brief Compute (but do not return) the ellipsoid volume.
*/
void computeVolume(void);

/**
@brief Compute (but do not return) the ellipsoid surface area.
*/
void computeArea(void);

/**
@brief Compute the ellipsoid's bounding spheres.

The "bounding spheres" are two spheres:
    - the minimum enclosing sphere, and
    - the maximum inscribed sphere
*/
void computeBoundingspheres(void);

/**
@brief Compute (but do not return) the ellipsoid dimensions.

@param triaxialcalc is true iff we want triaxial dimensions instead of (<i>L</i>,<i>W</i>,<i>T</i>)
*/
void computeDimensions(bool triaxialcalc);

/**
@brief Get the moment of inertia tensor

@return the moment of inertia tensor, <i>I</i>
*/
std::vector<Real> getI(void);

/**
@brief Writes formatted ellipsoid properties to standard out.
*/
void printProperties();

/**
@brief Returns a string description of the Ellipsoid class.

@return reference to string describing the class type
*/
std::string &getType(void) const { return (std::string &)("AXISYMMETRIC ELLIPSOID"); }
/**@}*/
    
/**
@brief Ellipsoid class assignment operator.

@param &that is the ellipsoid object on the left hand side
@return an ellipsoid object equal to that
*/
Ellipsoid& operator=(const Ellipsoid& that);
};

} // End of namespace tsquare
#endif
