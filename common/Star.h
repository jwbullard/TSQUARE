/**
@file
@brief Specifies the Star class for star-shaped particles.
*/
#ifndef STARH
#define STARH

#include "Particle.h"

namespace tsquare
{

/**
@class
@brief The Star class derives from the Particle abstract base class.

A 3D star-shaped object is one for which an interior point can be defined from
which a line drawn to any point on the surface remains entirely within the object.
That is, a single-valued distance function can be defined as a function <i>R</i> of
the polar and azimuthal angle coordinates (\f$\theta\f$,\f$\phi\f$) of a
spherical polar coordinate system.

All class members have been defined in the Particle base class, so
this class just overloads a lot of the methods.
*/

class Star : public Particle {

public:

/**
@brief Default constructor.
*/
Star();

/**
@brief Overloaded constructor.

The overloaded constructor will read in the SH coefficients from a file.

@param &fname is the name of the file containing the spherical harmonic coefficients
@param verbose is true iff verbose output is desired
*/
Star(std::string &fname, bool verbose);

/**
@brief Destructor.
*/
~Star() { std::cout << "Destroying Star"; }

/**
@name Integrated Geometry
@brief Functions for getting and setting integral shape features

Among the integral shape features defining the geometry of the particle,
we want to calculate
    - ASTM D 4791 dimensions (length, width, and thickness)
    - Triaxial dimensions (<i>D</i><sub>1</sub>,<i>D</i><sub>2</sub>,<i>D</i><sub>3</sub>)
    - Volume
    - Surface Area
    - Minimum enclosing and maximum inscribed spheres

The length (<i>L</i>), width (<i>W</i>), and thickness (<i>T</i>) are the
particle dimensions defined by ASTM D 4791.  Briefly, <i>L</i> is the length
of the longest line segment that can be drawn between two points both on the particle
surface; <i>W</i> is the length of the longest line segment perpendicular to <i>L</i> that
can be drawn between two points both on the particle surface; <i>T</i> is the length
of the longest line segment perpendicular to both <i>L</i> and <i>W</i> that can
be drawn between two points both on the particle surface.

Related to these are the ratios <i>L</i>/<i>T</i> and <i>W</i>/<i>T</i>

Krumbein (1941) defined the _triaxial_ dimensions of a particle, according to
which the triaxial length, <i>D</i><sub>1</sub> is the length of the longest
line segment that can be drawn between any two points both on the particle surface.
Therefore, <i>D</i><sub>1</sub> = <i>L</i>.  The triaxial width, <i>D</i><sub>2</sub>,
is the length of the longest line segment perpendicular to <i>D</i><sub>1</sub>
that can be drawn between two points on the perimeter of the maximum projected
area of the particle.  <i>D</i><sub>3</sub> is the length of the longest line
segment perpendicular to both <i>D</i><sub>1</sub> and <i>D</i><sub>2</sub>
between two points both on the particle surface.

The procedure for finding triaxial dimensions is:

    -# Rotate the length axis to be parallel to the <i>x</i>-axis
    -# Set up a fine mesh of Gaussian quadrature points on the particle
    -# Let \f$\eta\f$ be an angle of rotation about the <i>x</i>-axis.
    -# Scan from \f$\eta\f$ = 0 to \f$\pi\f$:
    -# At a given \f$\eta\f$, for each \f$\phi\f$:
        + Scan from \f$\theta\f$ = 0 to \f$\pi\f$, compute \f$\rho = x^2+y^2\f$
        + Keep track of the maximum value of \f$\rho\f$ and store it as \f$\rho(\phi)\f$.
    -# Do quadrature on the closed curve \f$\rho(\phi)\f$ and record the area, <i>A</i>(\f$\eta\f$).
    -# Keep track of the maximum value of <i>A</i> and the value of \f$\eta\f$ for it.
    -# Record the maximum dimension of the projected area perpendicular to the <i>x</i>-axis.
    -# When the maximum area is found, the maximum dimension from the previous step is
        <i>D</i><sub>2</sub>.

To rotate the particle about the <i>x</i>-axis using standard Euler angles, we first rotate about
the original <i>z</i>-axis by \f$\pi\f$/2, rotate by the desired angle \f$\beta\f$ about
the new <i>y</i>-axis, and then rotate about the new <i>z</i>-axis by -\f$\pi\f$/2.

It would be nice---and much faster---to just rotate the points themselves, rather
than the SH coefficients and then reconstructing the surface points from scratch.  But we do
it the slower way because it guarantees that the perimeter of the projected plane will be
parametrized as Gaussian quadrature points for easier integration.

To find the projected plane, we conduct a scan at each \f$\phi = \phi_i\f$, looking for the maximum
projected length onto the <i>xy</i>-plane of <i>r</i>(\f$\theta,\phi_i\f$) for
all possible values of \f$\theta\f$.

We are now done with finding the triaxial width, which was the hard part.  Now all
we need to do is compute the maximum triaxial thickness by rotating the particle back
to where its projected area was maximum, and then finding the maximum dimension
perpendicular to the new $yz$-plane.  We do two operations in sequence:

    -# Rotate about <i>x</i>-axis to `angle_at_max_projected_area` to compute the
       triaxial thickness, <i>D</i><sub>3</sub>.
       thickness, <i>T</i>.
*/

/**@{*/
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
virtual Real computeLength(std::vector<std::vector<Real> > sc, std::vector<Real> xc, std::vector<Real> yc,
                   std::vector<Real> zc, std::vector<Real> &lengthvec, Real &t1, Real &t2,
                   Real &p1, Real &p2);

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
virtual Real computeWidth(std::vector<std::vector<Real> > sc, std::vector<Real> xc,
                          std::vector<Real> yc, std::vector<Real> zc, std::vector<Real> lvec,
                          std::vector<Real> &wvec, Real &t1, Real &t2, Real &pp1, Real &pp2);

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
virtual Real computeThickness(std::vector<std::vector<Real> > sc, std::vector<Real> xc,
                          std::vector<Real> yc, std::vector<Real> zc, std::vector<Real> rc,
                          std::vector<Real> lvec, std::vector<Real> wvec, std::vector<Real> &tvec,
                          Real &t1, Real &t2, Real &pp1, Real &pp2);

/**
@brief Computes Cartesian coordinates for a set of surface points.

@param sc is the vector of \f$(\theta,\phi)\f$ pairs to consider
@param &xc will hold the <i>x</i> coordinate of each point in sc
@param &yc will hold the <i>y</i> coordinate of each point in sc
@param &zc will hold the <i>z</i> coordinate of each point in sc
@param &rc will hold the vector of distance from the center to each point
*/
void computeCoords(std::vector<std::vector<Real> > sc, std::vector<Real> &xc, std::vector<Real> &yc,
                   std::vector<Real> &zc, std::vector<Real> &lengthvec);

/**
@brief Computes the dimensions of the object.

This is an interface to other functions for computing the length, width
and thickness, or alternatively the triaxial dimensions <i>D<sub>1</sub></i>,
<i>D<sub>2</sub></i>, and <i>D<sub>3</sub></i>.  See the functions such
as `Star::computeLength` for details.

@param triaxialcalc is true iff the triaxial dimensions, instead of length,
       width, and thickness, should be calculated
*/
virtual void computeDimensions(bool triaxialcalc);

/**
@brief Computes (does not return) the object volume.
*/
virtual void computeVolume(void);

/**
@brief Computes (does not return) the object surface area.
*/
virtual void computeArea(void);

/**
@brief Computes and returns in the argument the object's centroid.
*/
virtual void computeCentroid(void);

/**
@brief Computes (does not return) the object's bounding spheres.

The two "bounding spheres" are
    - the minimum enclosing sphere
    - the maximum inscribed sphere

*/
virtual void computeBoundingspheres(void);

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
virtual Real computeTriaxialwidth(const int numpoints, const Real rstep, Real &rotate_angle,
                          Real &mpa, Real &ampa, Real &amaxwidth, Real &aamw);

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
virtual Real computeTriaxialthickness(const int numpoints, const Real init_angle,
                              const Real d2_angle, const Real w_angle,
                              Real &absmax_thickness);
/**@}*/

/**
@brief Reads the spherical harmonic coefficients from a file.

@param &fname is the name of the file containing the spherical harmonic coefficients
*/
void readSHcoeffs(std::string &fname);

/**
@brief Calculates the object's spherical harmonic coefficients.

This function calculates the SH coefficients of an object based on an
input voxel image.  The input file name is passed to the function as fname.
If the boolean flag surfaceMesh is true, then the file is immediately converted
to a CGAL gray level image, which can then be automatically surface meshed.  We
then use the surface mesh vertices as the reference points to calculate the centroid.
If surfaceMesh is false, we use only the voxel centers to get the reference points.

@param &fname is the name of the file containing the voxel image of the object
@param &outfilename is the file for writing the SH coefficients
@param numGausspts is the number of Gaussian quadrature points to use
@param maxDegree is the maximum degree of the SH expansion
@param surfaceMesh is true iff the input voxel image shall be immediately
    converted to a CGAL gray image for processing
*/
void calculateSHcoeffs(std::string &fname, std::string &outfilename,
                                 int numGausspts, int maxDegree, bool surfaceMesh);

/**
@brief Sets a mesh of surface points.

@param num is the number of points in each angular coordinate
@param nmax is the maximum degree of the SH expansion
@param gaussian is true iff Gaussian quadrature will be used for integration
*/
void setSurface(const int num, const int nmax, bool gaussian);

/**
@brief Populates the surface mesh points according to the Gaussian quadrature
points.

@param n is the number of Gaussian quadrature points in each angular coordinate
@return 0 if successful, non-zero otherwise
*/
int getGausspoints(int n);

/**
@brief Computes and sets length from the center to a point on the surface.

@param theta is the polar angle to use (radians)
@param phi is the azimuthal angle to use (radians)
*/
virtual void setR(const Real theta, const Real phi);

/**
@brief Computes and sets length from the center to a surface mesh point.

@param surfaceindex is the polar angle to use
*/
virtual void setR(const unsigned int surfaceindex);
    
/**
@brief Rotates the particle by Euler angles \f$(\alpha,\beta,\gamma)\f$

If we approximate a function by a spherical harmonic expansion on the unit
circle to degree $p$:
\f[
u(\theta,\phi) = \sum_{n=0}^p \sum_{m=-n}^n a_{nm} Y_n^m (\theta,\phi)
\f]
where \f$(\theta,\phi)\f$ is the position of a point <i>P</i> on the unit
circle, <i>u</i> is the value of the function at <i>P</i>, on the unit circle
in spherical polar coordinates, then a rotation of the coordinate axes is
represented as
\f[
u(\theta',\phi') = \sum_{n=0}^p \sum_{m'=-n}^n a_{nm'} Y_n^{m'} (\theta',\phi')
\f]
where \f$(\theta',\phi')\f$ is the position of <i>P</i> in the rotated
coordinate system.

The rotation is assumed to be implemented in the following manner.
Let <i>xyz</i> be the original coordinate system, <i>x'y'z'</i> be
the coordinate system after the first rotation, <i>x''y''z''</i> be
the coordinate system after the second rotation, and <i>XYZ</i> be the
final coordinate system.  We first rotate by an angle
\f$\alpha\f$ about the <i>z</i>-axis, then by an angle \f$\beta\f$ about
the <i>y'</i>-axis, and then by an angle \f$\gamma\f$ about the
<i>z''</i>-axis, which also happens to be the <i>Z</i>-axis.

The rotated spherical harmonic coefficients are given by:
\f[
a_{nm'} = \sum_{m=-n}^n D_n^{m',m} a_{nm},
\f]
where
\f[
D_n^{m',m} = e^{i m \gamma} d_n^{m',m}(\beta) e^{i m \alpha}
\f]

This method for rotating the spherical coefficients uses
three recursion relations advocated by Gimbutas and Greengard @cite Gimbutas09
\f{eqnarray*}
d_n^{m',m}(\beta) &= \cos^2 \left( \frac{\beta}{2} \right) \sqrt{\frac{(n+m)(n+m-1)}{(n+m')(n+m'-1)}}
d_{n-1}^{m'-1,m-1}(\beta) - 2 \sin \left( \frac{\beta}{2} \right) \cos \left( \frac{\beta}{2} \right)
\sqrt{\frac{(n+m)(n-m)}{(n+m')(n+m'-1)}}
d_{n-1}^{m'-1,m}(\beta) \\
&+ \sin^2 \left( \frac{\beta}{2} \right)
\sqrt{\frac{(n-m)(n-m-1)}{(n+m')(n+m'-1)}} d_{n-1}^{m'-1,m+1}(\beta) \\
d_n^{m',m}(\beta) &= \sin^2 \left( \frac{\beta}{2} \right) \sqrt{\frac{(n+m)(n+m-1)}{(n-m')(n-m'-1)}}
d_{n-1}^{m'+1,m-1}(\beta) + 2 \sin \left( \frac{\beta}{2} \right) \cos \left( \frac{\beta}{2} \right)
\sqrt{\frac{(n+m)(n-m)}{(n-m')(n-m'-1)}}
d_{n-1}^{m'+1,m}(\beta) \\
&+ \cos^2 \left( \frac{\beta}{2} \right)
\sqrt{\frac{(n-m)(n-m-1)}{(n-m')(n-m'-1)}} d_{n-1}^{m'+1,m+1}(\beta) \\
d_n^{m',m}(\beta) &= \sin \left( \frac{\beta}{2} \right) \cos \left( \frac{\beta}{2} \right)
\sqrt{\frac{(n+m)(n+m-1)}{(n+m')(n-m')}}
d_{n-1}^{m',m-1}(\beta)
+ \left( \cos^2 \left( \frac{\beta}{2} \right) - \sin^2 \left( \frac{\beta}{2} \right) \right)
\sqrt{\frac{(n-m)(n+m)}{(n-m')(n+m')}}
d_{n-1}^{m',m}(\beta) \\
&- \sin \left( \frac{\beta}{2} \right) \cos \left( \frac{\beta}{2} \right)
\sqrt{\frac{(n-m)(n-m+1)}{(n-m')(n+m')}} d_{n-1}^{m',m+1}(\beta)
\f}
This recursion procedure requires <i>O</i>(<i>p</i><sup>3</sup>) operations
instead of <i>O</i>(<i>p</i><sup>4</sup>) and is numerically stable
with minimal numerical error up to <i>p</i> = 40.

@param alpha is the angle \f$\alpha\f$ (radians) to rotate counterclockwise
       from above about the original <i>z</i> axis
@param beta is the Euler angle \f$\beta\f$ (radians) to rotate
       counterclockwise about the <i>x'</i> axis,
       which is the line of nodes
@param gamma is the Euler angle \f$\gamma\f$ (radians) to rotate
       about the new <i>z''</i> axis, which is also the <i>Z</i> axis
*/
void doRotate(Real alpha, Real beta, Real gamma);

/**
@brief Rotates the surface mesh decorating the object by
Euler angles \f$(\alpha,\beta,\gamma)\f$

Let <i>xyz</i> be the original coordinate system, <i>x'y'z'</i> be
the coordinate system after the first rotation, <i>x''y''z''</i> be
the coordinate system after the second rotation, and <i>XYZ</i> be the
final coordinate system.

@param alpha is the angle \f$\alpha\f$ (radians) to rotate counterclockwise
       from above about the original <i>z</i> axis
@param beta is the Euler angle \f$\beta\f$ (radians) to rotate
       counterclockwise about the <i>x'</i> axis,
       which is the line of nodes
@param gamma is the Euler angle \f$\gamma\f$ (radians) to rotate
       about the new <i>z''</i> axis, which is also the <i>Z</i> axis
*/
void doRotatesurfacepoints(Real alpha, Real beta, Real gamma);

/**
@name Deringing Filters
@brief Smooth surface artifacts due to truncating the SH expansion
*/
/**@{*/
/**
@brief Applies the Lanczos sigma factor to de-ring the object

Surface artifacts are introduced when using a truncated SH series to
approximate an object. These are called the Gibbs phenomena or
just "ringing".  We wish to reduce ringing on the surface while maintaining
its overall shape as much as possible.  This is accomplished with a windowing
technique that decreases the magnitude of the SH coefficients proportionally
as they get closer and closer to the cutoff degree <i>n<sub>max</sub></i>,
which in our code has the variable name `nmax_`.  One way to do
this windowing is to use Lanczos sigma factors @cite Sloan02.  Given a
SH coefficient <i>a<sub>n</sub><sup>m</sup></i>, its magnitude is multiplied by
\f[
\frac{\sin \left( \frac{\pi (n - n_0}{n_{max} - n_0} \right)}{ \left( \frac{\pi (n - n_0)}{n_{max} - n_0} \right)}
\f]
where <i>n<sub>0</sub></i> is the degree in the SH expansion at which we
should start the de-ringing process.

@param power is the exponent (how many times to apply the sigma factor)
@param nstart is the SH degree at which to start the de-ringing.
*/
void doLanczos(const int power, int nstart);

/**
@brief Applies the inverse of the Lanczos sigma factor.

If the parameters are identical to those used when `Star::doLanczos`
was applied, this function undoes the effect and restores the object to
its original state.

@param power is the exponent (how many times to apply the inverse sigma factor)
@param nstart is the SH degree at which to start.
*/
void undoLanczos(const int power, int nstart);

/**
@brief Applies the Hanning function to de-ring the object
 
An alternative to the Lanczos sigma factor is the Hanning filter, given by
\f[
\frac{1 + \cos \left( \frac{\pi (n - n_0)}{n_{max} - n_0} \right)}{2}
\f]
The Hanning filter decays faster than the Lanczos filter, but not quconst unsigned int surfindexite as fast
as the square of the Lanczos filter.

@param nstart is the SH degree at which to start the de-ringing.
*/
void doHanning(int nstart);

/**
@brief Applies the inverse of the Hanning function.

If the parameter is identical to that used when `Star::doHanning`
was applied, this function undoes the effect and restores the object to
its original state.

@param nstart is the SH degree at which to start.
*/
void undoHanning(int nstart);
/**@}*/

/**
@name Roundness Measures
@brief Functions to calculate various measures of roundness (opposite of angularity).

This section deals with measurements related to the _roundness_ of the particle, in
the spirit originally defined by Wadell @cite Wadell32 and subsequently elaborated by others.
Wadell's definition of roundness of a 2D silhouette or cross section was

\f[
\Omega \equiv \frac{\sum r_i}{N r_{\text{in}}}
\f]
where <i>r<sub>i</sub></i> is the radius of curvature of the
<i>i</i>-th corner, <i>r</i><sub>in</sub> is the radius of the
maximum circle that can be inscribed within the cross section, and <i>N</i>
is the number of corners identified in the cross section.

We would like to maintain the spirit of Wadell's definition for full 3D shapes,
yet characterizing all aspects of the 3D shape rather than just random cross sections.
For example, in a 2D cross section, the only "angular" feature is a so-called corner,
but in 3D there are corners and edges.  There appears to be no accepted definition for
3D roundness, so we are going to make up several:
\f[
\Omega \equiv \frac{1}{SA} \int_0^{2 \pi} \int_0^{\pi}
|\mathbf{X}(\theta,\phi) \cdot \mathbf{n}(\theta,\phi)|\, S\, d\theta\, d\phi
\f]
where <i>SA</i> is the particle surface area, <i><b>X</b></i> is the unit position
vector of the surface, <i><b>n</b></i> is the unit normal vector at any point on the
surface, and \f$S d\theta\, d\phi\f$ is a differential surface patch.
*/

/**@{*/
/**
@brief Computes the roundness of the object, using one of several shape factors.

    - `WADELL' (= 1) uses the 3D analog of Wadell's definition of roundness
    - `DOTPRODUCT') (= 0) uses the surface-weighted integral of the dot product
       of the position vector and the unit normal vector.
    - `ALLROUNDNESSES (= 2) computes both shape factors

@param method is the shape factor to use when computing the roundness
*/
virtual void computeRoundness(const unsigned int method);

/**
@brief Computes all roundness measures of the object

This function is simply an interface to the function of the same name
that computes a given measure of the roundness.  The present function
calls that function with the parameter set to compute all measures
of roundness.
*/
virtual void computeRoundness(void) { computeRoundness(ALLROUNDNESSES); }
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
@brief Gets the distance from the center to a surface mesh point.

@param surfaceindex is the index of the surface mesh point
@return the distance to the surface mesh point
*/
virtual Real getR(const unsigned int surfaceindex);

/**
@brief Gets the distance from the center to any point on the surface.

@param theta is the polar angle \f$\theta\f$ (radians)
@param phi is the azimuthal angle \f$\phi\f$ (radians)
@return the distance from the center to the surface in the direction
         \f$(\theta,\phi)\f$
*/
virtual Real getR(const Real theta, const Real phi);

/**
@brief Gets the distances from the center to each surface point.
@return the vector storing the distances from the center to each point
         in the surface mesh, relative to the center
*/
virtual std::vector<Real> getR(void) const { return r_; }

/**
@brief Gets the derivative \f$R_{\phi}\f$.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\phi}\f$
*/
virtual Real getRp(Real theta, const Real phi);

/**
@brief Gets the derivative \f$R_{\phi}\f$.

@param surfindex is the index of the point in the surface mesh
@return the derivative \f$R_{\phi}\f$
*/
virtual Real getRp(const unsigned int surfindex);

/**
@brief Gets the derivative \f$R_{\theta}\f$.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\theta}\f$
*/
virtual Real getRt(Real theta, const Real phi);

/**
@brief Gets the derivative \f$R_{\theta}\f$.

@param surfindex is the index of the point in the surface mesh
@return the derivative \f$R_{\theta}\f$ at the point
*/
virtual Real getRt(const unsigned int surfindex);

/**
@brief Gets the second derivative \f$R_{\phi \phi}\f$.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\phi \phi}\f$
*/
virtual Real getRpp(Real theta, const Real phi);

/**
@brief Gets the second derivative \f$R_{\phi \phi}\f$.

@param surfindex is the index of the point in the surface mesh
@return the derivative \f$R_{\phi \phi}\f$ at the point
*/
virtual Real getRpp(const unsigned int surfindex);

/**
@brief Gets the second derivative \f$R_{\theta \theta}\f$.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\theta \theta}\f$
*/
virtual Real getRtt(Real theta, const Real phi);

/**
@brief Gets the second derivative \f$R_{\theta \theta}\f$.

@param surfindex is the index of the point in the surface mesh
@return the derivative \f$R_{\theta \theta}\f$ at the point
*/
virtual Real getRtt(const unsigned int surfindex);

/**
@brief Gets the second derivative \f$R_{\theta \phi}\f$.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\theta \phi}\fconst unsigned int surfindex$
*/
virtual Real getRtp(Real theta, const Real phi);

/**
@brief Gets the second derivative \f$R_{\theta \phi}\f$.

@param surfindex is the index of the point in the surface mesh
@return the derivative \f$R_{\theta \phi}\f$ at the point
*/
virtual Real getRtp(const unsigned int surfindex);
/**@}*/

/**
@name X vector
@brief \f$\vec{X}(\theta,\phi)\f$ is a vector pointing from the center
to a point on the surface with angular coordinates \f$(\theta,\phi)\f$.

In differential geometry we view the polar and azimuthal angles
as parameters that specify the three components of the Cartesian vector
\f$\mathbf{X}\f$, which is a vector in \f$\Re^3\f$:
\f[
\mathbf{X}(\theta,\phi) = \left(x(\theta,\phi),y(\theta,\phi),z(\theta,\phi) \right)
\f]
where
\f{eqnarray*}
x(\theta,\phi) &=& R(\theta,\phi) \sin \theta \cos \phi \\
y(\theta,\phi) &=& R(\theta,\phi) \sin \theta \sin \phi \\
z(\theta,\phi) &=& R(\theta,\phi) \cos \theta
\f}

Therefore, we now provide functions that will calculate the
$\mathbf{X}$ with respect to the parameters $\theta$ and $\phi$.
For convenience, we specify functions that return the entire vector at
once, as well as functions that return only the individual components.
*/

/**@{*/
/**
@brief Gets a component of the vector pointing from the center
       a point on the surface.

@param index specifies which coordinate to get
    - 0 is the <i>x</i>-coordinate
    - 1 is the <i>y</i>-coordinate
    - 2 is the <i>z</i>-coordinate
@param theta is the polar angle \f$\theta\f$ at the surface point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the surface point (radians)
@return the requested vector component
*/
Real getX(const unsigned int index, const Real theta, const Real phi);

/**
@brief Gets the vector \f$\vec{X}(\theta,\phi)\f$.

@param theta is the polar angle \f$\theta\f$ at the surface point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the surface point (radians)
@return the requested vector \f$\vec{X}(\theta,\phi)\f$
*/
virtual std::vector<Real> getX(const Real theta, const Real phi);

/**
@brief Gets the vector \f$\vec{X}\f$ for a surface mesh point.

@param surfindex is the index of the surface mesh point
@return the requested vector \f$\vec{X}\f$
*/
virtual std::vector<Real> getX(const unsigned int surfindex);

/**
@brief Gets the <i>x</i>-coordinate of every surface point
@return the vector storing the <i>x</i>-coordinate of each surface point
         in the surface mesh, relative to the center
*/
virtual std::vector<Real> getX(void) const { return x_; }

/**
@brief Gets the derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to azimuthal angle, \f$X_{i,\phi}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\phi}\f$
*/
Real getXp(const unsigned int index, const Real theta, const Real phi);

/**
@brief Gets the derivative of \f$\vec{X}\f$ with respect to
azimuthal angle, \f$\vec{X}_{\phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\phi}\f$
*/
std::vector<Real> getXp(const Real theta, const Real phi);

/**
@brief Gets the derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to azimuthal angle, \f$X_{i,\phi}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\phi}\f$
*/
Real getXp(const unsigned int index, const unsigned int surfaceindex);

/**
@brief Gets the derivative of \f$\vec{X}\f$ with respect to
azimuthal angle, \f$\vec{X}_{\phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\phi}\f$
*/
std::vector<Real> getXp(const unsigned int surfaceindex);

/**
@brief Gets the derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to polar angle, \f$X_{i,\theta}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\theta}\f$
*/
Real getXt(const unsigned int index, const Real theta, const Real phi);

/**
@brief Gets the derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\theta}\f$
*/
std::vector<Real> getXt(const Real theta, const Real phi);

/**
@brief Gets the derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to polar angle, \f$X_{i,\theta}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\theta}\f$
*/
Real getXt(const unsigned int index, const unsigned int surfaceindex);

/**
@brief Gets the derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\theta}\f$
*/
std::vector<Real> getXt(const unsigned int surfaceindex);

/**
@brief Gets the second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to azimuthal angle, \f$X_{i,\phi \phi}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\phi \phi}\f$
*/
Real getXpp(const unsigned int index, const Real theta, const Real phi);

/**
@brief Gets the second derivative of \f$\vec{X}\f$ with respect to
azimuthal angle, \f$\vec{X}_{\phi \phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\phi \phi}\f$
*/
std::vector<Real> getXpp(const Real theta, const Real phi);

/**
@brief Gets the second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to azimuthal angle, \f$X_{i,\phi \phi}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\phi \phi}\f$
*/
Real getXpp(const unsigned int index, const unsigned int surfaceindex);

/**
@brief Gets the second derivative of \f$\vec{X}\f$ with respect to
azimuthal angle, \f$\vec{X}_{\phi \phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\phi \phi}\f$
*/
std::vector<Real> getXpp(const unsigned int surfaceindex);

/**
@brief Gets the second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to polar angle, \f$X_{i,\theta \theta}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\theta \theta}\f$
*/
Real getXtt(const unsigned int index, const Real theta, const Real phi);

/**
@brief Gets the second derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta \theta}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\theta \theta}\f$
*/
std::vector<Real> getXtt(const Real theta, const Real phi);

/**
@brief Gets the second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to polar angle, \f$X_{i,\theta \theta}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\theta \theta}\f$
*/
Real getXtt(const unsigned int index, const unsigned int surfaceindex);

/**
@brief Gets the second derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta \theta}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\theta \theta}\f$
*/
std::vector<Real> getXtt(const unsigned int surfaceindex);

/**
@brief Gets the mixed second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$, \f$X_{i,\theta \phi}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\theta \phi}\f$
*/
Real getXtp(const unsigned int index, const Real theta, const Real phi);

/**
@brief Gets the mixed second derivative of \f$\vec{X}\f$,
\f$\vec{X}_{\theta \phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\theta \phi}\f$
*/
std::vector<Real> getXtp(const Real theta, const Real phi);

/**
@brief Gets the mixed second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$, \f$X_{i,\theta \phi}\f$.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\theta \phi}\f$
*/
Real getXtp(const unsigned int index, const unsigned int surfaceindex);

/**
@brief Gets the mixed second derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta \phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\theta \phi}\f$
*/
std::vector<Real> getXtp(const unsigned int surfaceindex);
/**@}*/

/**
@name First Fundamental Form
@brief Compute and return components of the first fundamental
       form of differential geometry.

The first fundamental form deals with the ability to calculate
the coordinates of a vector in the tangent space of the surface at a particular point
parameterized by the \f$(\theta,\phi)\f$.  The crux of the first fundamental
form is the matrix
\f[
\left( \begin{array}{cc} E & F \\ F & G \end{array} \right)
\f]
where \f$E = \langle \mathbf{X_{\theta}},\mathbf{X_{\theta}} \rangle\f$,
\f$F = \langle \mathbf{X_{\theta}},\mathbf{X_{\phi}} \rangle\f$,
and \f$G = \langle \mathbf{X_{\phi}},\mathbf{X_{\phi}} \rangle\f$.
*/
/**@{*/
/**
@brief Gets the first fundamental form <i>E</i> component of a point
       on the surface with angular coordinates \f$(\theta,\phi)\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the first fundamental form component <i>E</i>
*/
Real getE(const Real theta, const Real phi);

/**
@brief Gets the first fundamental form <i>E</i> component of a point
       in the surface mesh.

@param si is the index of the point in the surface mesh
@return the first fundamental form component <i>E</i>
*/
Real getE(const unsigned int si);

/**
@brief Gets the first fundamental form <i>F</i> component of a point
       on the surface with angular coordinates \f$(\theta,\phi)\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the first fundamental form component <i>F</i>
*/
Real getF(const Real theta, const Real phi);

/**
@brief Gets the first fundamental form <i>F</i> component of a point
       in the surface mesh.

@param si is the index of the point in the surface mesh
@return the first fundamental form component <i>F</i>
*/
Real getF(const unsigned int si);

/**
@brief Gets the first fundamental form <i>G</i> component of a point
       on the surface with angular coordinates \f$(\theta,\phi)\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the first fundamental form component <i>G</i>
*/
Real getG(const Real theta, const Real phi);

/**
@brief Gets the first fundamental form <i>G</i> component of a point
       in the surface mesh.

@param si is the index of the point in the surface mesh
@return the first fundamental form component <i>G</i>
*/
virtual Real getG(const unsigned int si);
/**@}*/

/**
@name Second Fundamental Form
@brief Compute and return components of the second fundamental
       form of differential geometry.

The second fundamental form is a quadratic form that
relates to the normal curvature at a point on a surface.  It is
calculated as the inner product of the differential of normal vector with
a vector in the tangent space.  The crux of the second fundamental form
is the matrix
\f[
\left( \begin{array}{cc} L & M \\ M & N \end{array} \right)
\f]
where \f$L = -\langle \mathbf{\hat{n}},\mathbf{X_{\theta \theta}} \rangle\f$,
\f$M = \langle \mathbf{\hat{n}},\mathbf{X_{\theta \phi}} \rangle\f$,
and \f$N = \langle \mathbf{\hat{n}},\mathbf{X_{\phi \phi}} \rangle\f$.
*/

/**
@name Normal Vector
@brief Functions for computing the normal vector at a point on the surface.
Provides methods to calculate the components of the unit normal vector at
a surface point parameterized by \f$(\theta,\phi)\$:
\f[
\mathbf{\hat{n}} =
     \frac{\mathbf{X_{\theta}} \times \mathbf{X_{\phi}}}{\sqrt{EG - F^2}}
\f]
where <i>E</i>, <i>F</i>, and <i>G</i> are the components of the first
fundamental form.
*/

/**@{*/
/**
@brief Gets the <i>i</i>-th component of the normal vector at a point
       given by its index in the surface mesh.

@param index is the vector component to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param si is the index of the point in the surface mesh
@return the requested component of the normal vector
*/
Real getNormal(const unsigned int index, const unsigned int si);

/**
@brief Gets the <i>i</i>-th component of the normal vector at a point
       given by angular coordinates \f$(\theta,\phi)\f$.

@param index is the vector component to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the requested component of the normal vector
*/
Real getNormal(const unsigned int index, const Real theta, const Real phi);

/**
@brief Gets the normal vector at a point
       given by its index in the surface mesh.

@param si is the index of the point in the surface mesh
@return the normal vector at the point
*/
std::vector<Real> getNormal(const unsigned int si);

/**
@brief Gets the normal vector at a point
       given by angular coordinates \f$(\theta,\phi)\f$.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the normal vector at the point
*/
std::vector<Real> getNormal(const Real theta, const Real phi);
/**@}*/

/**
@name Curvature
@brief Functions for the integrated Gaussian and mean curvatures.

The Gaussian curvature, <i>K</i> and mean curvature, <i>H</i>, are defined
in terms of the principal curvatures \f$\kappa_1\f$ and \f$\kappa_2\f$ as
\f{eqnarray*}
K &\equiv& \kappa_1 \kappa_2 \\
H &\equiv& \frac{1}{2} \left( \kappa_1 + \kappa_2 \right)
\f}
These two curvatures are calculated using a combination of the first and
second fundamental forms:
\f{eqnarray*}
K &=& \frac{LN - M^2}{EG - F^2} \\
H &=& \frac{1}{2} \left( \frac{LG - 2 M F + NE}{EG - F^2} \right)
\f}

Even more importantly, though, we can compute directly the (unit) principal
directions and the principal curvatures using the first and second
fundamental forms.  The two unit principal direction vectors are
\f[
    \left( \frac{\lambda_{+}}{\sqrt{\lambda_{+}^2 + 1}},\frac{1}{\sqrt{\lambda_{+}^2 + 1}} \right)
\qquad \text{and} \qquad 
    \left( \frac{\lambda_{-}}{\sqrt{\lambda_{-}^2 + 1}},\frac{1}{\sqrt{\lambda_{-}^2 + 1}} \right)
\f]
where
\f[
\lambda_{\pm} = \frac{-b \pm \sqrt{ b^2 - 4ac}}{2a}
\f]
with the definitions <i>a</i> = <i>FN</i> - <i>GM</i>, <i>b</i> =
<i>EN</i> - <i>GL</i>, and <i>c</i> = <i>EM</i> - <i>FL</i>.
Furthermore, the principal curvatures \f$\kappa_{+}\f$ and \f$\kappa_{-}\f$ are
\f[
\kappa_{\pm} = \frac{-\beta \pm \sqrt{ \beta^2 - 4\alpha \gamma}}{2\alpha}
\f]
with the definitions \f$\alpha\f$ = <i>EG</i> - <i>F</i><sup>2</sup>,
\f$\beta\f$ = 2<i>FM</i> - <i>EN</i> - <i>GL</i>, and
\f$\gamma\f$ = <i>LN</i> - <i>M</i><sup>2</sup>.
Note, however, that the subscripts for the principal directions and
principal curvatures are not necessarily related to each other.
*/

/**@{*/
/**
@brief Computes and returns the mean curvature integrated over the surface.

@return the mean curvature, <i>H</i>, integrated over the surface
*/
Real getIntegratedH(void);

/**
@brief Computes and returns the Gaussian curvature integrated over the surface.

@return the Gaussian curvature, <i>K</i>, integrated over the surface
*/
Real getIntegratedK(void);
/**@}*/

/**
@brief Computes and returns the moment of inertia tensor

@param itensor is the vector of upper diagonal elements of moment of inertia tensor
@param rg2itensor is the vector of upper diagonal elements of gyration tensor
@param isphere is the moment of inertia of sphere with same volume as particle
@param rg2sphere is the square of the radius of gyratio of sphere with same volume as particle
@return 0 upon successful completion
*/
int getI(std::vector<Real> &itensor, std::vector<Real> &rg2tensor, Real &isphere, Real &rg2sphere);

/**
@name Output
@brief Functions for producing visualizations or formatted output of properties.
*/

/**@{*/
/**
@brief Outputs a VRML file for viewing the object.

@param vrmlname is the name of the VRML file to write
@param surface_poly is the CGAL polyhedron approximation of the surface
@param render_inscribed_sphere is true iff the inscribed sphere should be drawn
@param render_enclosing_sphere is true iff the enclosing sphere should be drawn
@param render_hull is true iff the convex hull should be drawn
@param surface_color_scheme specifies how to color the surface (NOT USED)
@param rotation is the angle to rotate about the y axis in degrees
@return 0 if file created successfully, non-zero otherwise
*/
virtual int createVRML(const std::string vrmlname, bool render_inscribed_sphere,
             bool render_enclosing_sphere, bool render_hull,
             int surface_color_scheme, Real rotation);

/**
@brief Creates a file of the SH coefficients for the particle in its current state.

@param shname is the name of the file to write
@return 0 if successful
*/
virtual int createSHFile(const std::string shname);

/**
@brief Makes a digitized representation of the particle in a bounding box

@param digitizedname is the name of the image file to write
@param resolution is the voxel edge length relative to the unit of measurement
@return 0 if successful
*/
virtual int createDigitizedImage(const std::string digitizedname, Real resolution);

/**
@brief Returns a string description of the Star class.

@return reference to string describing the class type
*/
std::string &getType(void) const { return (std::string &)("STAR-SHAPED PARTICLE"); }

/**
@brief Produce formatted output of the particle properties.
*/
virtual void printProperties();
/**@}*/

/**
@name Operators.
@brief Overloaded operators for the Particle class
*/

/**@{*/
/**
@brief Overloads the assignment operator.

@param &that is a reference to the object on the right side of the assignment
@return a reference to the Star object assigned
*/
Star& operator=(const Star& that);
/**@}*/

};

} // end of tsquare namespace

#endif
