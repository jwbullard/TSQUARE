/**
@file
@brief Header for main differential geometry toolbox program.
*/
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <cstdlib>
#include <locale>
#include <Particle.h>
#include <Star.h>
#include <Ellipsoid.h>
#include <global.h>
#ifdef PARALLEL
#include "mpi.h"
#endif

/**
@brief Master function for detailed analysis of individual particle
*/
int diffgeom(void);

/**
@brief Master function for acquiring shape statistics on a collection of particles
*/
int shapestats(int argc, char *argv[]);

/**
@brief Gets the excluded volume of the object.

@warning Currently does not work
@todo Complete or destroy
*/
Real getExcludedvolume(tsquare::Particle *p);

/**
@brief Gets the sum of two vectors.

@param v1 is the first vector (STL container)
@param v2 is the second vector (STL container)
@param @result will hold the sum of the two vectors
*/
void getVectorsum(std::vector<Real> v1, std::vector<Real> v2, std::vector<Real> &result);

/**
@brief Prints curvature(s) at a particular point on the surface.

@param pointid is the index of the point in the surface mesh
@param npoints is the number of points in the surface mesh
@param *p is a pointer to the particle
@param all is true iff curvatures will be printed for every point in the mesh
*/
void printCurvature(int pointid, int npoints, tsquare::Particle *p, bool all);

/**
@brief Prints curvature(s) for points along a surface path.

@param *p is a pointer to the particle
*/
void printCurvatureOnAPath(tsquare::Particle *p);

/**
@brief Prints curvature on a portion of a "latitude".

For star-shaped particles, the term "latitude" is correctly used because the polar
angle will be held constant and the azimuthal angle varied to create the path.
For non-star-shaped particles, the term "latitude" is misleading; the path will be
one for which the first surface parameter is held constant and the other one varied.

@todo Check if this function will work for non-star-shaped particles.

@param *p is a pointer to the particle
@param ctheta is the constant value of the first surface parameter
@param astart is the beginning value of the second surface parameter
@param aend is the ending value of the second surface parameter
@param ainc is the amount to increment the second surface parameter along the path
@param &frname is the file name to writing the results
*/
void printCurvatureOnALatitude(tsquare::Particle *p, Real ctheta, Real astart, Real aend, Real ainc, std::string &frname);

/**
@brief Prints curvature on a portion of a "longitude".

For star-shaped particles, the term "longitude" is correctly used because the azimuthal
angle will be held constant and the polar angle varied to create the path.
For non-star-shaped particles, the term "longitude" is misleading; the path will be
one for which the first surface parameter varies while the second is held constant.

@todo Check if this function will work for non-star-shaped particles.

@param *p is a pointer to the particle
@param cphi is the constant value of the second surface parameter
@param astart is the beginning value of the first surface parameter
@param aend is the ending value of the first surface parameter
@param ainc is the amount to increment the first surface parameter along the path
@param &frname is the file name to writing the results
*/
void printCurvatureOnALongitude(tsquare::Particle *p, Real cphi, Real astart, Real aend, Real ainc, std::string &frname);

/**
@brief A structure for containing all the "geom" input data for a particle.

The "geom" input data are the items in a single row of the "geom" file
corresponding to the particle being analyzed, including its name,
some macroscopic dimensions, and integrated quantities like volume and
surface area.
*/
typedef struct lineitem {

    std::string name;     /**< The name of the file holding the SH coefficients for the particle. */
    Real xlow;       /**< The minimum <i>x</i> coordinate of the particle */
    Real xhi;        /**< The maximum <i>x</i> coordinate of the particle */
    Real ylow;       /**< The minimum <i>y</i> coordinate of the particle */
    Real yhi;        /**< The maximum <i>y</i> coordinate of the particle */
    Real zlow;       /**< The minimum <i>z</i> coordinate of the particle */
    Real zhi;        /**< The maximum <i>z</i> coordinate of the particle */
    Real volume;     /**< Particle volume */
    Real surfarea;   /**< Particle surface area */

/**
@brief Normalized surface area.

In this case, the normalized quantity is the surface area of the object divided
by the surface area of a sphere having the same volume as the particle.  Therefore,
the normalized surface area is always greater than or equal to unity.
*/
    Real nsurfarea;

/**
@brief Volume equivalent spherical diameter (VESD).

The VESD is the diameter of a sphere having the same volume as the particle.
*/
    Real diam;

    Real Itrace;      /**< Trace of the moment-of-inertia tensor */
    int Nnn;	      /**< SH expansion degree needed to get within 5% of Gaussian curvature */

/**
@brief Normalized Gaussian curvature.

The total curvature of any closed surface (Euler characteristic = 2) in 3D must equal
4\f$\pi\f$ according to the Gauss-Bonnet theorem.  The normalized Gaussian curvature in
this case is the computed total curvature divided by 4\f$\pi\f$, so it should always be
near unity.
*/
    Real NGC;         

    Real length;      /**< Length as defined by ASTM D 4791 */
    Real width;       /**< Width as defined by ASTM D 4791 */
    Real thickness;   /**< Thickness as defined by ASTM D 4791 */
    Real nlength;     /**< Length divided by thickness */
    Real nwidth;      /**< Width divided by thickness */
} Lineitem;
