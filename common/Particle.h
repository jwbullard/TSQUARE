/**
@file
@brief Specifies the Particle class.
*/
#ifndef PARTICLEH
#define PARTICLEH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <list>
#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include "ExtendedVector.h"
#include "Sphere.h"
#include "global.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/ImageIO.h>
#include <CGAL/IO/io.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
// #include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include "imagemanipulation.h"
// #include "AffineTransformation.h"

/**
@struct
@brief Calculates unit normal vector of a CGAL plane facet.
*/
struct Facet_normal {
    template <class Facet>
    void operator()( Facet& f) {
        typename Facet::Halfedge_handle h = f.halfedge();
        typename Facet::Normal_3 normal = CGAL::cross_product(
          h->next()->vertex()->point() - h->vertex()->point(),
          h->next()->next()->vertex()->point() - h->next()->vertex()->point());
        f.normal() = normal / std::sqrt( normal * normal);
    }
};

/**
@struct
@brief Calculates unit normal vector associated with a CGAL polyhedral vertex.

This function examines all the planar facets that share the vertex, evaluates
their normal vectors and then averages the collection.
*/
struct Vertex_normal {
    template <class Vertex>
    void operator()( Vertex& v) {
        typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
        typedef typename Vertex::Halfedge_around_vertex_const_circulator Circ;
        Circ c = v.vertex_begin();
        Circ d = c;
        CGAL_For_all( c, d) {
            if ( ! c->is_border())
                normal = normal + c->facet()->normal();
        }
        v.normal() = normal / std::sqrt( normal * normal);
    }
};

/**
@class
@brief Builds more functionality into CGAL's halfedge vertex class.

This ultimately will be used to define a customized polyhedral mesh representation
for the surface of the object.
*/
template <class Refs, class T, class P, class Norm>
class My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P> {
    Norm  norm;
public:
    My_vertex() {} // repeat mandatory constructors
    My_vertex( const P& pt) : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt) {}
    typedef Norm Normal_3;
    Normal_3&       normal()       { return norm; }
    const Normal_3& normal() const { return norm; }
};

/**
@class
@brief Builds more functionality into CGAL's polyhedral facet class.

This ultimately will be used to define a customized polyhedral mesh representation
for the surface of the object.
*/
template <class Refs, class T, class Norm>
class My_facet : public CGAL::HalfedgeDS_face_base<Refs, T> {
    Norm  norm;
public:
    // no constructors to repeat, since only default constructor mandatory
    typedef Norm Normal_3;
    Normal_3&       normal()       { return norm; }
    const Normal_3& normal() const { return norm; }
};

/**
@struct
@brief A customized wrapper for binding a vertex and associated face of a polyhedron

This ultimately will be used to define a customized polyhedral mesh representation
for the surface of the object.
*/
struct My_items : public CGAL::Polyhedron_items_3 {
    template <class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
        typedef My_vertex<Refs, CGAL::Tag_true, Point, Normal> Vertex;
    };
    template <class Refs, class Traits>
    struct Face_wrapper {
        typedef typename Traits::Vector_3 Normal;
        typedef My_facet<Refs, CGAL::Tag_true, Normal> Face;
    };
};

/**
@name CGAL Typedefs
@brief Shorthand names for specifying CGAL classes relevant to the problem at hand.

See CGAL's documentation for details on these implementations and what they mean.
*/

/**@{*/
/**
@typedef
@brief The CGAL kernel to use.
*/
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

/**
@typedef
@brief A customized CGAL 3D polyhedron
*/
typedef CGAL::Polyhedron_3<K,My_items> Polyhedron;
// typedef CGAL::Polyhedron_3<K> Polyhedron;

/**
@typedef
@brief The CGAL 3D tetrahedron
*/
typedef CGAL::Tetrahedron_3<K> Tetrahedron;

/**
@typedef
@brief The CGAL surface triangulation of a 3D point set.
*/
typedef CGAL::Triangulation_3<K> Triangulation_3;

/**
@typedef
@brief The CGAL surface mesh (Delaunay) triangulation of a 2D surface mesh embedded in 3D.
*/
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

/**
@typedef
@brief Stores a Delaunay triangulation complex.
*/
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

/**
@typedef
@brief Geometric trait class of a Delaunay triangulation.
*/
typedef Tr::Geom_traits GT;

/**
@typedef
@brief CGAL class for handling gray-level images of input data
*/
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;

/**
@typedef
@brief CGAL class for handling an implicit surface constructed from a gray-level image.
*/
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

/**
@typedef
@brief An iterator-like object for traversing cells of a CGAL surface mesh.
*/
typedef Triangulation_3::Cell_handle Cell_handle;

/**
@typedef
@brief A point representation of a vertex in a CGAL surface triangulation.
*/
typedef Triangulation_3::Point Point;

/**
@typedef
@brief CGAL's implementation of a point in 3D Cartesian space.
*/
typedef K::Point_3 CartesianPoint;

/**
@typedef
@brief CGAL's implementation of a plane in 3D Cartesian space.
*/
typedef K::Plane_3 CartesianPlane;

/**
@typedef
@brief CGAL's implementation of a 3D vector in Cartesian space.
*/
typedef K::Vector_3 CartesianVector;

/**
@typedef
@brief CGAL's implementation of a line segment in 3D Cartesian space.
*/
typedef K::Segment_3 CartesianSegment;

/**
@typedef
@brief Supports the AABB tree traversal concept for rapid searches.

This concept provides rapid searches for intersections between triangulated surfaces
and line segments or rays.
*/
// typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;

/**
@typedef
@brief Supports the AABB tree traversal concept for rapid searches.

This concept provides rapid searches for intersections between triangulated surfaces
and line segments or rays.
*/
typedef CGAL::AABB_traits<K,Primitive> Traits;

/**
@typedef
@brief Supports the AABB tree traversal concept for rapid searches.

This concept provides rapid searches for intersections between triangulated surfaces
and line segments or rays.
*/
typedef CGAL::AABB_tree<Traits> Tree;

/**
@typedef
@brief Supports 3D affine transformations

The class Aff_transformation_3 represents three-dimensional affine transformations
including translation, reflection, rotation, and scaling
*/
typedef CGAL::Aff_transformation_3<K> Aff_transformation_3;

/**
@typedef
@brief Supports the AABB tree traversal concept for rapid searches.

This concept provides rapid searches for intersections between triangulated surfaces
and line segments or rays.
*/
typedef Tree::Object_and_primitive_id Object_and_primitive_id;

/**
@typedef
@brief Supports the AABB tree traversal concept for rapid searches.

This concept provides rapid searches for intersections between triangulated surfaces
and line segments or rays.
*/
typedef Tree::Primitive_id Primitive_id;

/**
@typedef
@brief A CGAL iterator over vertices in a surface mesh or polyhedron.
*/
typedef Polyhedron::Vertex_iterator Vertex_iterator;

/**
@typedef
@brief A CGAL iterator over facets in a surface mesh or polyhedron.
*/
typedef Polyhedron::Facet_iterator Facet_iterator;

/**
@typedef
@brief CGAL's implementation of a surface facet in a mesh or triangulation.
*/
typedef Polyhedron::Facet Facet;

/**
@typedef
@brief A type of iterator for scanning facets in a surface triangulation.
*/
typedef Polyhedron::Facet_handle Facet_handle;

/**
@typedef
@brief An CGAL iterator for scanning the half edges around a surface facet of a polyhedron.
*/
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;

/**
@typedef
@brief An CGAL circulator for scanning the half edges around a
surface facet of a polyhedron.

A circulator is like an iterator, except it is used for collections of objects that are
known to cycle back to the beginning, such as the edges around a surface facet.
*/
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

/**
@typedef
@brief An CGAL circulator for scanning the half edges around a
surface facet of a polyhedron.

A circulator is like an iterator, except it is used for collections of objects that are
known to cycle back to the beginning, such as the edges common to a surface vertex.
*/
typedef Polyhedron::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;

/**
@typedef
@brief CGAL's implementation of a point in 3D.
*/
typedef K::Point_3 Point_3;
/**@}*/


template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b)
{
    assert(a.size() == b.size());
    std::vector<T> result;
    result.clear();
    result.reserve(a.size());
    std::transform(a.begin(),a.end(),b.begin(),std::back_inserter(result),std::plus<T>());
    return result;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b)
{
    assert(a.size() == b.size());
    std::vector<T> result;
    result.clear();
    result.reserve(a.size());
    std::transform(a.begin(),a.end(),b.begin(),std::back_inserter(result),std::minus<T>());
    return result;
}

namespace tsquare
{

/**
@class
@brief The Particle base class upon which all other types of shapes are based.

The Particle class provides data and methods for implementing a spherical harmonic (SH)
representation of the surface of a 3D object.  The class provides functions for
    - computing the SH coefficients of a particle contained in a digital image
    - reading pre-computed SH coefficients from a file
    - reconstructing particles and operating on them, including
        + dilation
        + rotation
    - calculating integral geometry characteristics such as volume, surface area, and
      moment of intertia tensor
    - calculation differential geometry properties such as first and second fundamental forms,
      local principal curvatures and principal directions, and local Gaussian and mean curvature
    - computing the convex hull
    - computing the maximum inscribed sphere and minimum enclosing sphere
    - computing measures of sphericity and roundness
    - creating 3D visualizations of a particle using VRML

The class relies quite heavily both on the CGAL library for computational geometry (http://www.cgal.org),
and on the Boost library (http://www.boost.org).
*/
class Particle {

protected:
 
std::string name_;            /**< the particle's name */
Real xlo_;                    /**< the particle's minimum <i>x</i> coordinate */
Real ylo_;                    /**< the particle's minimum <i>y</i> coordinate */
Real zlo_;                    /**< the particle's minimum <i>z</i> coordinate */
Real xhi_;                    /**< the particle's maximum <i>x</i> coordinate */
Real yhi_;                    /**< the particle's maximum <i>y</i> coordinate */
Real zhi_;                    /**< the particle's maximum <i>z</i> coordinate */
Real volume_;                 /**< the particle's volume */
Real area_;                   /**< the particle's surface area */
Real narea_;                  /**< the surface area normalized to that of a sphere of the same volume */
Real diam_;                   /**< the particle's volume equivalent spherical diameter (VESD) */
Real itrace_;                 /**< the trace of the moment-of-inertia tensor */
Real ngc_;                    /**< normalized Gaussian curvature */
std::vector<Real> maxroundness_;   /**< the maximum roundness factor found on the particle's surface */
std::vector<Real> minroundness_;   /**< the minimum roundness factor found on the particle's surface */
std::vector<Real> dim_;            /**< the particle's length, width, and thickness measures */
std::vector<Real> ndim_;           /**< the particle's length/thickness and width/thickness ratios */
std::vector<Real> triaxialdim_;    /**< the particle's triaxial dimensions as defined by Wadell @cite Wadell32 */
bool isgaussian_;             /**< true iff surface points are distributed according to Gaussian weightings */
bool verbose_;                /**< true iff extra output is desired */
Sphere inscribedsphere_;      /**< the particle's maximum inscribed sphere (a `Sphere` object) */
Sphere enclosingsphere_;      /**< the particle's minimum enclosing sphere (a `Sphere` object) */
Polyhedron convexhull_;       /**< the convex hull, represented as a CGAL polyhedron object */
Real hullvolume_;             /**< volume of the particle's convex hull */
Real hullarea_;               /**< surface area of the particle's convex hull */
int nmax_;                    /**< maximum degree of the spherical harmonic expansion */
int nmax_from_anmfile_;       /**< maximum degree of SH expansion given in an input file */
std::vector<std::vector<ExtendedVector> > a_;   /**< the collection of complex spherical harmonic coefficients */
int ntheta_;                          /**< number of surface points in the polar angular coordinate */
int nphi_;                            /**< number of surface points in the azimuthal angular coordinate */
std::vector<std::vector<Real> > surface_;       /**< vector of all points decorating the particle surface */
std::vector<std::vector<Real> > normal_;        /**< vector of normal vectors for each point on the particle surface */
std::vector<Real> r_;                      /**< collection of distances from center to each point on the surface */
std::vector<Real> x_;                      /**< collection of <i>x</i> components of vector from center to surface */
std::vector<Real> y_;                      /**< collection of <i>y</i> components of vector from center to surface */
std::vector<Real> z_;                      /**< collection of <i>z</i> components of vector from center to surface */
std::vector<Real> centroid_;               /**< Cartesian coordinates of particle centroid */

/**
@brief vector of differential surface patches used in integrals over the surface.

The differential surface patch in spherical polar coordinates is
\f[
ds = \sqrt{R_{\phi}^2 + \left( R_{\theta}^2 + R^2 \right) \sin^2 \theta}
\f]
*/
std::vector<Real> ddd2_;
std::vector<Real> cumroundness_;           /**< the integrated, cumulative, roundness measure for the particle */
std::vector<std::vector<Real> > roundness_;     /**< the local roundness at each point on the particle surface */
std::vector<std::vector<Real> > pcurv_;         /**< vector of principal curvatures at each point on the surface */

/**
@brief Particle convexity.

The particle's convexity, <i>C</i>, is defined here as the volume of the particle divided by the volume of its
convex hull:

\f[
C = \frac{V}{V_\text{hull}}
\f]
*/
Real convexity_;

std::vector<Real> xg_;                      /**< vector of Gaussian points for Gaussian quadrature */
std::vector<Real> wg_;                      /**< vector of weightings for Gaussian quadrature */
    
public:
    
/**
@brief Default constructor.
*/
Particle();
    
/**
@brief Virtual destructor.
*/
virtual ~Particle() { std::cout << "Destroying Particle"; }

/**
@brief Sets the particle name.

@param &str is the particle name to assign
*/
void setName(const std::string &str) { name_ = str; }

/**
@brief Sets the minimum <i>x</i> coordinate of the particle.

@param xl is the coordinate value to assign
*/
void setXlo(const Real xl) { xlo_ = xl; }

/**
@brief Sets the minimum <i>y</i> coordinate of the particle.

@param yl is the coordinate value to assign
*/
void setYlo(const Real yl) { ylo_ = yl; }

/**
@brief Sets the minimum <i>z</i> coordinate of the particle.

@param zl is the coordinate value to assign
*/
void setZlo(const Real zl) { zlo_ = zl; }

/**
@brief Sets the maximum <i>x</i> coordinate of the particle.

@param xh is the coordinate value to assign
*/
void setXhi(const Real xh) { xhi_ = xh; }

/**
@brief Sets the maximum <i>y</i> coordinate of the particle.

@param yh is the coordinate value to assign
*/
void setYhi(const Real yh) { yhi_ = yh; }

/**
@brief Sets the maximum <i>z</i> coordinate of the particle.

@param zh is the coordinate value to assign
*/
void setZhi(const Real zh) { zhi_ = zh; }

/**
@brief Sets one of the spherical harmonic coefficients.

@note This setter is specific to star-shaped particles.

@param i is the degree of the coefficient to set
@param j is the order of the coefficient to set
@param realvalue is the real component of the (complex) SH coefficient
@param imaginaryvalue is the imaginary component of the (complex) SH coefficient
*/
void setA(const int i, const int j, const Real realvalue, const Real imaginaryvalue)
{
    std::complex<Real> cmplx(realvalue,imaginaryvalue);
    a_[RDIM][i][j] = cmplx;
}

/**
@brief Sets one of the spherical harmonic coefficients.

@note This setter accommodates spherical harmonic parameterization
      functionality, in which there is a separate SH expansion for
      each coordinate in a Cartesian system.

@param coord is the Cartesian coordinate to use
    - 0 = <i>x</i>
    - 1 = <i>y</i>
    - 2 = <i>z</i>
@param i is the degree of the coefficient to set
@param j is the order of the coefficient to set
@param realvalue is the real component of the (complex) SH coefficient
@param imaginaryvalue is the imaginary component of the (complex) SH coefficient
*/
void setA(const int coord, const int i, const int j,
          const Real realvalue, const Real imaginaryvalue)
{
    std::complex<Real> cmplx(realvalue,imaginaryvalue);
    a_[coord][i][j] = cmplx;
}

/**
@brief Sets one of the spherical harmonic coefficients.

@note This setter is specific to star-shaped particles.

@param i is the degree of the coefficient to set
@param j is the order of the coefficient to set
@param complex value is the complex value to set
*/
void setA(const int i, const int j, std::complex<Real> complexvalue)
{
    a_[RDIM][i][j] = complexvalue;
}

/**
@brief Sets one of the spherical harmonic coefficients.

@note This setter accommodates spherical harmonic parameterization
      functionality, in which there is a separate SH expansion for
      each coordinate in a Cartesian system.

@param coord is the Cartesian coordinate to use
    - 0 = <i>x</i>
    - 1 = <i>y</i>
    - 2 = <i>z</i>
@param i is the degree of the coefficient to set
@param j is the order of the coefficient to set
@param complex value is the complex value to set
*/
void setA(const int coord, const int i, const int j, std::complex<Real> complexvalue)
{
    a_[coord][i][j] = complexvalue;
}

/**
@brief Sets the boolean variable signifying verbose output if true.

@param verbose is true iff verbose output is desired
*/
void setVerbose(bool verbose) { verbose_ = verbose; }

/**
@brief Reads all the SH coefficients from a file.

@note Pure virtual function; must be overridden by derived classes.

@param &fname is the file name holding the SH coefficients
*/
virtual void readSHcoeffs(std::string &fname) = 0;

/**
@brief Calculates the object's spherical harmonic coefficients.

This function calculates the SH coefficients of an object based on an
input voxel image.  The input file name is passed to the function as fname.
If the boolean flag surfaceMesh is true, then the file is immediately converted
to a CGAL gray level image, which can then be automatically surface meshed.  We
then use the surface mesh vertices as the reference points to calculate the centroid.
If surfaceMesh is false, we use only the voxel centers to get the reference points.

@note Pure virtual function; must be overridden by derived classes.

@param &fname is the name of the file containing the voxel image of the object
@param &outfilename is the file for writing the SH coefficients
@param numGausspts is the number of Gaussian quadrature points to use
@param maxDegree is the maximum degree of the SH expansion
@param surfaceMesh is true iff the input voxel image shall be immediately
    converted to a CGAL gray image for processing
*/
virtual void calculateSHcoeffs(std::string &fname, std::string &outfilename,
                                 int numGausspts, int maxDegree, bool surfaceMesh) = 0;

/**
@brief Adds a point to the set of points decorating the surface.

Often times we will want to anchor certain points on the surface of the
particle for reference.  This will be needed when we want to approximate
the surface with overlapping spheres, or when we want to perform quadrature.
The surface_ list stores the ordered pairs \f$(\theta_i,\phi_i)\f$.  We will
require this functionality whether or not the particle is star shaped.

@note The two parameters passed are generic surface parameters.  For star-shaped
      objects, the two parameters can be interpreted as the polar and azimuthal
      angles, (\f$\theta\f$,\f$\phi\f$) in a spherical polar coordinate system.
      For convenience, we keep
      the names `theta` and `phi` for the surface parameters even when the object
      is not star-shaped.

@param theta is the value of the first parameter that parameterizes the surface
@param phi is the value of the second parameter that parameterizes the surface
*/
void addtoSurface(Real theta, const Real phi);
 
/**
@brief Sets a mesh of surface points.

@note Pure virtual function; must be overridden by derived classes.

@param num is the number of points in each angular coordinate
@param nmax is
@param gaussian is true iff Gaussian quadrature will be used for integration
*/
virtual void setSurface(const int num, const int nmax, bool gaussian) = 0;

/**
@brief Populates the surface mesh points according to the Gaussian quadrature
points.

@note Pure virtual function; must be overridden by derived classes.
@param n is the number of Gaussian quadrature points in each angular coordinate
@return 0 if successful, non-zero otherwise
*/
virtual int getGausspoints(int n) = 0;

/**
@brief Gets the Legendre polynomials of degree <i>n-1>, <i>n</i>, and <i>n+1</i>
evaluated at <i>x</i>.

@param n is the degree of the Legendre polynomial to evaluate
@param *x is a pointer to the argument to the Legendre polynomial
@param *pn is the value of degree <i>n</i> of the Legendre polynomial
@param *pnm1 is the value of degree <i>n-1</i> of the Legendre polynomial
@param *pnp1 is the value of degree <i>n+1</i> of the Legendre polynomial
*/
void getLegendre(int n, Real *x, Real *pn, Real *pnm1, Real *pnp1);

/**
@brief Gets the Legendre polynomials of degree <i>n</i> and
order <i>m</i> evaluated at <i>u</i>.

@param n is the degree of the Legendre polynomial to evaluate
@param m is the order of the Legendre polynomial to evaluate
@param u is the argument to the Legendre polynomial
@return the value of the Legendre polynomial <i>P<sub>n</sub><sup>m</sup></i>(<i>u</i>)
*/
Real getLegendre(const unsigned int n, const int m, Real u)
{
    return (boost::math::legendre_p(n,m,u));
}

    
/**
@brief Gets the particle name.

@return the particle name
*/
std::string &getName() const { return (std::string &)name_; }

/**
@brief Gets the minimium <i>x</i> coordinate of the particle.

@return the minimum <i>x</i> coordinate
*/
Real getXlo() const { return xlo_; }

/**
@brief Gets the minimium <i>y</i> coordinate of the particle.

@return the minimum <i>y</i> coordinate
*/
Real getYlo() const { return ylo_; }

/**
@brief Gets the minimium <i>z</i> coordinate of the particle.

@return the minimum <i>z</i> coordinate
*/
Real getZlo() const { return zlo_; }

/**
@brief Gets the maximium <i>x</i> coordinate of the particle.

@return the maximum <i>x</i> coordinate
*/
Real getXhi() const { return xhi_; }

/**
@brief Gets the maximium <i>y</i> coordinate of the particle.

@return the maximum <i>y</i> coordinate
*/
Real getYhi() const { return yhi_; }

/**
@brief Gets the maximium <i>z</i> coordinate of the particle.

@return the maximum <i>z</i> coordinate
*/
Real getZhi() const { return zhi_; }

/**
@brief Gets the particle volume.

@return the particle volume
*/
Real getVolume() const { return volume_; }

/**
@brief Gets the particle surface area.

@return the particle surface area
*/
Real getArea() const { return area_; }

/**
@brief Gets the particle surface area normalized to that of a sphere of the same volume.

@return the normalized surface area surface area
*/
Real getNarea() const { return narea_; }

/**
@brief Gets the volume equivalent spherical diameter (VESD).

The VESD is the diameter of a sphere having the same volume as the particle.
@return the volume equivalent spherical diameter
*/
Real getDiam() const { return diam_; }

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
    -# Rotate about <i>x</i>-axis to `angle_at_absmax_width` to compute the usual
       thickness, <i>T</i>.
*/

/**@{*/
/**
@brief Compute the length of the particle.

The length, <i>L</i>, is the length of the longest line segement that
can be drawn entirely within the object. 

@note Pure virtual function; must be overridden.

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
                   Real &p1, Real &p2) = 0;

/**
@brief Compute the refined length of the particle.

The length, <i>L</i>, is the length of the longest line segment that
can be drawn between two points both on the particle surface.
The function `computeLength`
gets a coarse estimate of the length by finding the range of polar and
azimuthal angles that must contain the line segement.  It then calls
this function to find the exact position of the line segement within
those angular ranges.  This function returns the refined length to the
calling function, which then returns it its calling function.

@note This function is common to all particles (star-shaped or not), because
      the surface points are passed in terms of their Cartesian coordinates.

@param sf is the vector of \f$(\theta,\phi)\f$ pairs to consider in the fine set
@param sc is the vector of \f$(\theta,\phi)\f$ pairs to consider in the coarse set
@param xf is the <i>x</i> coordinate of each point in sf (the fine set)
@param yf is the <i>y</i> coordinate of each point in sf (the fine set)
@param zf is the <i>z</i> coordinate of each point in sf (the fine set)
@param xc is the <i>x</i> coordinate of each point in sc (the coarse set)
@param yc is the <i>y</i> coordinate of each point in sc (the coarse set)
@param zc is the <i>z</i> coordinate of each point in sc (the coarse set)
@param rf is the distance from the center to each point in sf (the fine set)
@param rc is the distance from the center to each point in sc (the coarse set)
@param &lengthvec will hold the length vector
@param rm is a real number (WHY WHAT IS THIS?)
@param &tt1 is the estimated lower bound on the polar angle, \f$\theta_1\f$ (radians)
@param &tt2 is the estimated upper bound on the polar angle, \f$\theta_2\f$ (radians)
@param &ppp1 is the estimated lower bound on the azimuthal angle, \f$\phi_1\f$ (radians)
@param &ppp2 is the estimated upper bound on the azimuthal angle, \f$\phi_2\f$ (radians)
@return the refined estimate of the length, <i>L</i>
*/
Real computeRefinedlength(std::vector<std::vector<Real> > sf, std::vector<std::vector<Real> > sc,
                          std::vector<Real> xf, std::vector<Real> yf, std::vector<Real> zf,
                          std::vector<Real> xc, std::vector<Real> yc, std::vector<Real> zc,
                          std::vector<Real> rf, std::vector<Real> rc, std::vector<Real> &lengthvec,
                          Real rm, Real &tt1, Real &tt2, Real &ppp1, Real &ppp2);

/**
@brief Compute the width of the particle.

The width, <i>W</i>, is the length of the longest line segment perpendicular
to <i>L</i> that can be drawn between any two points both on the particle surface.

@note Pure virtual function; must be overridden.

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
                          std::vector<Real> &wvec, Real &t1, Real &t2, Real &pp1, Real &pp2) = 0;
/**
@brief Compute the refined width of the particle.

The width, <i>W</i>, is the length of the longest line segment perpendicular
to <i>L</i> that can be drawn between any two points both on the particle surface.
The function `computeWidth` gets a coarse estimate of the width by finding the
range of polar and azimuthal angles that must contain the line segement.  It
then calls this function to find the exact position of the line segement within
those angular ranges.  This function returns the refined length to the
calling function, which then returns it its calling function.

@note This function is common to all particles (star-shaped or not), because
      the surface points are passed in terms of their Cartesian coordinates.

@param sf is the vector of \f$(\theta,\phi)\f$ pairs to consider in the fine set
@param sc is the vector of \f$(\theta,\phi)\f$ pairs to consider in the coarse set
@param xf is the <i>x</i> coordinate of each point in sf (the fine set)
@param yf is the <i>y</i> coordinate of each point in sf (the fine set)
@param zf is the <i>z</i> coordinate of each point in sf (the fine set)
@param xc is the <i>x</i> coordinate of each point in sc (the coarse set)
@param yc is the <i>y</i> coordinate of each point in sc (the coarse set)
@param zc is the <i>z</i> coordinate of each point in sc (the coarse set)
@param rf is the distance from the center to each point in sf (the fine set)
@param rc is the distance from the center to each point in sc (the coarse set)
@param lvec is the length vector
@param &wvec will hold the refined width vector
@param rm is (WHY WHAT IS THIS?)
@param ddot is the dot product of the length and width vectors
@param &tt1 is the estimated lower bound on the polar angle, \f$\theta_1\f$ (radians)
@param &tt2 is the estimated upper bound on the polar angle, \f$\theta_2\f$ (radians)
@param &ppp1 is the estimated lower bound on the azimuthal angle, \f$\phi_1\f$ (radians)
@param &ppp2 is the estimated upper bound on the azimuthal angle, \f$\phi_2\f$ (radians)
@return the refined estimate of the width, <i>W</i>
*/
Real computeRefinedwidth(std::vector<std::vector<Real> > sf, std::vector<std::vector<Real> > sc,
                          std::vector<Real> xf, std::vector<Real> yf, std::vector<Real> zf,
                          std::vector<Real> xc, std::vector<Real> yc, std::vector<Real> zc,
                          std::vector<Real> rf, std::vector<Real> rc, std::vector<Real> lvec,
                          std::vector<Real> &wvec, Real rm, Real ddot, Real &tt1,
                          Real &tt2, Real &ppp1, Real &ppp2);
/**
@brief Compute the thickness of the particle.

The thickness, <i>T</i>, is the length of the longest line segment perpendicular
to both <i>L</i> and <i>W</i> that can be drawn between any two points both on the
particle surface.

@note Pure virtual function; must be overridden.

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
                          Real &t1, Real &t2, Real &pp1, Real &pp2) = 0;

/**
@brief Compute the refined thickness of the particle.

The thickness, <i>T</i>, is the length of the longest line segment perpendicular
to both <i>L</i> and <i>W</i> that can be drawn between any two points both on the
particle surface.  The function `computeThickness`
gets a coarse estimate of the width by finding the range of polar and
azimuthal angles that must contain the line segement.  It then calls
this function to find the exact position of the line segement within
those angular ranges.  This function returns the refined length to the
calling function, which then returns it its calling function.

@note This function is common to all particles (star-shaped or not), because
      the surface points are passed in terms of their Cartesian coordinates.

@param sf is the vector of \f$(\theta,\phi)\f$ pairs to consider in the fine set
@param sc is the vector of \f$(\theta,\phi)\f$ pairs to consider in the coarse set
@param xf is the <i>x</i> coordinate of each point in sf (the fine set)
@param yf is the <i>y</i> coordinate of each point in sf (the fine set)
@param zf is the <i>z</i> coordinate of each point in sf (the fine set)
@param xc is the <i>x</i> coordinate of each point in sc (the coarse set)
@param yc is the <i>y</i> coordinate of each point in sc (the coarse set)
@param zc is the <i>z</i> coordinate of each point in sc (the coarse set)
@param rf is the distance from the center to each point in sf (the fine set)
@param rc is the distance from the center to each point in sc (the coarse set)
@param lvec is the length vector already found
@param wvec is the refined width vector already found
@param &tvec will hold the refined thickness vector
@param rm is (WHY WHAT IS THIS?)
@param ddotl is the dot product of the length and thickness vectors
@param ddotw is the dot product of the width and thickness vectors
@param &tt1 is the estimated lower bound on the polar angle, \f$\theta_1\f$ (radians)
@param &tt2 is the estimated upper bound on the polar angle, \f$\theta_2\f$ (radians)
@param &ppp1 is the estimated lower bound on the azimuthal angle, \f$\phi_1\f$ (radians)
@param &ppp2 is the estimated upper bound on the azimuthal angle, \f$\phi_2\f$ (radians)
@return the refined estimate of the thickness, <i>T</i>
*/
Real computeRefinedthickness(std::vector<std::vector<Real> > sf, std::vector<std::vector<Real> > sc,
                          std::vector<Real> xf, std::vector<Real> yf, std::vector<Real> zf,
                          std::vector<Real> xc, std::vector<Real> yc, std::vector<Real> zc,
                          std::vector<Real> rf, std::vector<Real> rc, std::vector<Real> lvec,
                          std::vector<Real> wvec, std::vector<Real> &tvec, Real rm,
                          Real ddotl, Real ddotw, Real &tt1, Real &tt2,
                          Real &ppp1, Real &ppp2);

/**
@brief Gets the length (<i>L</i>), width (<i>W</i>), or thickness (<i>T</i>).

@param index specifies the dimension to be returned
    - 0 = <i>L</i>
    - 1 = <i>W</i>
    - 2 = <i>T</i>
@return the requested dimension
*/
Real getDim(const unsigned int index)
{
    try {
        return dim_.at(index);
    }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","getDim","dim_",
                            dim_.size(),index);
        ex.printException();
        exit(1);
    }
}

/**
@brief Gets all three dimensions, (<i>L</i>,<i>W</i>,<i>T</i>).

@return a vector storing <i>L</i>,<i>W</i>, and <i>T</i>.
*/
std::vector<Real> getDim() const { return dim_; }

/**
@brief Gets the ratio <i>L</i>/<i>T</i> or <i>W</i>/<i>T</i>.

@param index specifies the ratio to be returned
    - 0 = <i>L</i>/<i>T</i>
    - 1 = <i>W</i>/<i>T</i>
    - 2 = <i>T</i>/<i>T</i> (= 1)
@return the requested dimension
*/
Real getNdim(const unsigned int index)
{
    try {
        return ndim_.at(index);
    }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","getNdim","ndim_",
                            ndim_.size(),index);
        ex.printException();
        exit(1);
    }
}

/**
@brief Gets all three ratios, (<i>L</i>/<i>T</i>,<i>W</i>/<i>T</i>,<i>T</i>/<i>T</i>).

@return a vector storing <i>L</i>/<i>T</i>,<i>W</i>/<i>T</i>, and <i>T</i>/<i>T</i>.
*/
std::vector<Real> getNdim() const { return ndim_; }

/**
@brief Gets the length, <i>L</i>.

@warning Length must already be computed using `computeLength` for this to be meaningful, unless
         the dimension has already been set by some other means.

@return the length
*/
virtual Real getLength() const
{
    try {
        return dim_.at(0);
    }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","getLength","dim_",
                            dim_.size(),0);
        ex.printException();
        exit(1);
    }
}

/**
@brief Gets the width, <i>W</i>.

@warning Width must already be computed using `computeWidth` for this to be meaningful, unless
         the dimension has already been set by some other means.

@return the width
*/
virtual Real getWidth() const
{
    try {
        return dim_.at(1);
    }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","getWidth","dim_",
                            dim_.size(),1);
        ex.printException();
        exit(1);
    }
}

/**
@brief Gets the thickness, <i>T</i>.

@warning Thickness must already be computed using `computeThickness` for this to be meaningful, unless
         the dimension has already been set by some other means.

@return the thickness
*/
virtual Real getThickness() const
{
    try {
        return dim_.at(2);
    }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","getThickness","dim_",
                            dim_.size(),2);
        ex.printException();
        exit(1);
    }
}

/**
@brief Gets the ratio <i>L</i>/<i>T</i>.

@return the ratio <i>L</i>/<i>T</i>
*/
Real getNlength() const
{
    try {
        return ndim_.at(0);
    }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","getNlength","ndim_",
                            ndim_.size(),0);
        ex.printException();
        exit(1);
    }
}

/**
@brief Gets the ratio <i>W</i>/<i>T</i>.

@return the ratio <i>W</i>/<i>T</i>
*/
Real getNwidth() const
{
    try {
        return ndim_.at(1);
    }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","getNwidth","ndim_",
                            ndim_.size(),1);
        ex.printException();
        exit(1);
    }
}

/**
@brief Gets the ratio <i>T</i>/<i>T</i>.

@return the ratio <i>T</i>/<i>T</i> (= 1)
*/
Real getNthickness() const
{
    try {
        return ndim_.at(2);
    }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","getNthickness","ndim_",
                            ndim_.size(),2);
        ex.printException();
        exit(1);
    }
}

/**
@brief Computes Cartesian coordinates for a set of surface points.

@note Pure virtual function; must be overridden.

@param sc is the vector of \f$(\theta,\phi)\f$ pairs to consider
@param &xc will hold the <i>x</i> coordinate of each point in sc
@param &yc will hold the <i>y</i> coordinate of each point in sc
@param &zc will hold the <i>z</i> coordinate of each point in sc
@param &lengthvec will hold the vector of distance from the center to each point
*/
virtual void computeCoords(std::vector<std::vector<Real> > sc, std::vector<Real> &xc, std::vector<Real> &yc,
                   std::vector<Real> &zc, std::vector<Real> &lengthvec) = 0;
/**
@brief Computes the dimensions of the object.

This is an interface to other functions for computing the length, width
and thickness, or alternatively the triaxial dimensions <i>D<sub>1</sub></i>,
<i>D<sub>2</sub></i>, and <i>D<sub>3</sub></i>.  See the functions such
as `Particle::computeLength` for details.

@note Pure virtual function; must be overridden by derived classes.

@param triaxialcalc is true iff the triaxial dimensions, instead of length,
       width, and thickness, should be calculated
*/
virtual void computeDimensions(bool triaxialcalc) = 0;

/**
@brief Computes (does not return) the object volume.

@note Pure virtual function; must be overridden by derived classes.
*/
virtual void computeVolume(void) = 0;

/**
@brief Computes (does not return) the object surface area.

@note Pure virtual function; must be overridden by derived classes.
*/
virtual void computeArea(void) = 0;

/**
@brief Computes and returns in the argument the object's centroid.

@note Pure virtual function; must be overridden by derived classes.
*/
virtual void computeCentroid(void) = 0;

/**
@brief Computes (does not return) the object's bounding spheres.

The two "bounding spheres" are
    - the minimum enclosing sphere
    - the maximum inscribed sphere

@note Pure virtual function; must be overridden by derived classes.
*/
virtual void computeBoundingspheres(void) = 0;

/**
@brief Gets the radius ratio of the enclosing sphere to the inscribed sphere.

@return the radius ratio of the enclosing sphere to the inscribed sphere.
*/
Real getBoundingsphereratio(void) const
{
    if (enclosingsphere_.getRadius() > 0.0) {
        return (inscribedsphere_.getRadius() / enclosingsphere_.getRadius());
    } else {
        throw FloatException("Particle","getBoundingsphereratio",
                             "Divide by zero enclosing sphere radius");
    }
}

/**
@brief Gets (does not compute) the maximum inscribed sphere.

@return the maximum inscribed sphere (a `Sphere' object)
*/
Sphere getMaxinscribedsphere(void) const { return inscribedsphere_; }

/**
@brief Gets (does not compute) the minimum enclosing sphere.

@return the minimum enclosing sphere (a `Sphere' object)
*/
Sphere getMinenclosingsphere(void) const { return enclosingsphere_; }

/**
@brief Compute the triaxial width of the particle.

The triaxial width, <i>D<sub>2</sub></i>, is the length of the longest
line segment that can be drawn in the plane of the object's maximum cross-sectional
area _and_ perpendicular to the length <i>L</i>.

@note Pure virtual function; must be overridden.

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
                          Real &mpa, Real &ampa, Real &amaxwidth, Real &aamw) = 0;

/**
@brief Compute the triaxial thickness of the particle.

The triaxial thickness, <i>D<sub>3</sub></i>, is the length of the
longest line segment that can be drawn inside the object while simultaneously
being perpendicular to both <i>L</i> _and_ <i>D<sub>2</sub></i>.

@note Pure virtual function; must be overridden.

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
                              Real &absmax_thickness) = 0;

/**
@brief Get the triaxial length, width, or thickness.

@param index specifies which dimension shall be returned
    - 0 = triaxial length (<i>D</i><sub>1</sub> = <i>L</i>)
    - 1 = triaxial width (<i>D</i><sub>2</sub>)
    - 2 = triaxial thickness (<i>D</i><sub>3</sub>)
@return the requested triaxial dimensions
*/
Real getTriaxialdim(const unsigned int index)
{
    try {
        return triaxialdim_.at(index);
    }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","getTriaxialdim","triaxialdim_",
                            triaxialdim_.size(),index);
        ex.printException();
        exit(1);
    }
    return triaxialdim_[index];
}

/**
@brief Get the vector of triaxial dimensions.

@return a vector holding the triaxial dimensions
*/
std::vector<Real> getTriaxialdim(void) const
{
    return triaxialdim_;
}
/**@}*/

/**
@brief Gets the trace of the moment-of-inertia tensor.

@return the trace of the moment-of-inertia tensor
*/
Real getItrace() const { return itrace_; }

/**
@brief Gets the normalized Gaussian curvature, <i>K</i>.

@todo Figure out what "normalized" Gaussian curvature means.
*/
Real getNgc() const { return ngc_; }

/**
@brief Gets the maximum degree of the SH expansion.

@return the maximum degree of the SH expansion.
*/
int getNmax() const { return nmax_; }

/**
@brief Gets the maximum degree of the SH expansion according to the
       traditional geom file.

The geom file specifies an upper limit for the SH expansion degree based
on deviations of the integrated Gaussian curvature from the theoretical
value of 4\f$\pi\f$.  Evidently, above some critical degree, the SH expansion
produces shapes having total curvatures that deviate more and more from the
theoretical value.

@return the maximum degree of the SH expansion.
*/
int getNmax_from_anmfile() const { return nmax_from_anmfile_; }

/**
@brief Gets the entire vector of SH coefficients.

@return the vector of complex SH coefficients
*/
std::vector<std::vector<ExtendedVector> > &getA() const
{
    return (std::vector<std::vector<ExtendedVector> >&) a_;
}

/**
@brief Gets the number of points in the partitioning of the first surface parameter.

@return the number of points in the partitioning of the first surface parameter.
*/
int getNtheta() const { return ntheta_; }

/**
@brief Gets the number of points in the partitioning of the second surface parameter.

@return the number of points in the partitioning of the second surface parameter.
*/
int getNphi() const { return nphi_; }

/**
@brief Gets the entire collection of mesh points decorating the surface.

@return the collection of surface mesh points
*/
std::vector<std::vector<Real> > getSurface() const { return surface_; }

/**
@brief Gets all the surface points in the current partitioning that fall within 
upper and lower bounds of the two surface parameters.

@param thetamin is the lower bound on the first surface parameter
@param thetamax is the upper bound on the first surface parameter
@param phimin is the lower bound on the second surface parameter
@param phimax is the upper bound on the second surface parameter
@param &ssub will hold the collection of surface points in the subset
*/
void getSurfacesubset(Real thetamin, Real thetamax, Real phimin, Real phimax, std::vector<int> &ssub);

/**
@brief Is verbose output set?

@return true iff verbose output has been specified
*/
bool getVerbose() const { return verbose_; }

/**
@brief Gets the value of first surface parameter for a point in the surface mesh.

@param index is the index of the point in the vector of surface mesh points
@return the value of the first surface parameter for the point
*/
Real getTheta(const unsigned int index)
{
     try {
        return surface_.at(index).at(0);
     }
     catch (std::out_of_range &oor) {
         EOBException ex("Particle","getTheta","surface_",
                            surface_.size(),index);
        ex.printException();
        exit(1);
     }
}

/**
@brief Gets the value of second surface parameter for a point in the surface mesh.

@param index is the index of the point in the vector of surface mesh points
@return the value of the second surface parameter for the point
*/
Real getPhi(const unsigned int index)
{
     try {
        return surface_.at(index).at(1);
     }
     catch (std::out_of_range &oor) {
         EOBException ex("Particle","getPhi","surface_",
                            surface_.size(),index);
        ex.printException();
        exit(1);
     }
}
    
/**
@brief Constructs a meshed surface polyhedron of facets based on an image.

Takes a digital image of an object and
creates a 3D polyhedral representation of the surface.  The voxel values
in the original digital image can be binary (0 or 1) or floating point values
on the closed interval [0,1]
 
@param &fname is the name of the file holding the gray level image
@param centroid holds the Cartesian coordinates of the centroid
@param xBB is the <i>x</i> dimensions of a box guaranteed to hold the object
@param yBB is the <i>y</i> dimensions of a box guaranteed to hold the object
@param zBB is the <i>z</i> dimensions of a box guaranteed to hold the object
@return the CGAL Polyhedron object representing the surface
*/
Polyhedron getSurfacePolyhedronFromImage(std::string &fname, std::vector<Real> centroid,
                                         Real xBB, Real yBB, Real zBB);
/**
@brief Dilates the particle.

@param dr is the linear factor by which to multiply each linear dimension of the particle.
*/
virtual void doDilate(const Real dr)
{
    std::complex<Real> cdr(dr,0.0);
    for (int k = 0; k < a_.size(); k++) {
        for (int i = 0; i <= nmax_; i++) {
            for (int j = -i; j <= i; j++) {
               a_[k][i][j] *= cdr;
            }
        }
    }
}

/**
@brief Rotates the particle.

Rotations are specified by the three Euler angles.
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
three recursion relations advocated by Gimbutas and Greengard @cite Gimbutas09:
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

These functions will be slightly different depending on whether we are doing SH for each
Cartesian coordinate (nonstar shapes) or only for the position vector (star shapes).  Therefore,
again we will supply only pure virtual functions here and will defer the details to the individual
classes.  We could define them here, but we can take advantage of numerical efficiency for star
shapes that won't pertain to nonstar shapes.  However, we will keep the function for rotating
the decorated points on the surface because that function is independent of particle shape.

@note Pure virtual function; must be overridden by derived classes.

@param alpha is the first Euler angle, \f$\alpha\f$ (radians)
@param alpha is the second Euler angle, \f$\beta\f$ (radians)
@param alpha is the third Euler angle, \f$\gamma\f$ (radians)
*/
virtual void doRotate(Real alpha, Real beta, Real gamma) = 0;

/**
@brief Rotates all the points in the surface mesh.

If points have been decorated on the surface as reference points to particular
locations, we rotate them when the particle is
rotated, just for the convenience to the user of being able to fix the index of
a surface site with a spot on the particle regardless of its rotation.  However,
the details of the rotation will depend on whether we are using a spherical
polar coordinate system (for star-shaped particles) or a more general two-parameter
surface parametrization.  Therefore, the function itself has been declared as a
pure virtual function.

Rotations are specified by the three Euler angles.
The rotation is assumed to be implemented in the following manner.
Let <i>xyz</i> be the original coordinate system, <i>x'y'z'</i> be
the coordinate system after the first rotation, <i>x''y''z''</i> be
the coordinate system after the second rotation, and <i>XYZ</i> be the
final coordinate system.  We first rotate by an angle
\f$\alpha\f$ about the <i>z</i>-axis, then by an angle \f$\beta\f$ about
the <i>y'</i>-axis, and then by an angle \f$\gamma\f$ about the
<i>z''</i>-axis, which also happens to be the <i>Z</i>-axis.

@note Pure virtual function currently; must be overridden by derived classes.
@todo This should be common to all classes; pure virtual function is not needed here.
@param alpha is the first Euler angle, \f$\alpha\f$ (radians)
@param alpha is the second Euler angle, \f$\beta\f$ (radians)
@param alpha is the third Euler angle, \f$\gamma\f$ (radians)
*/
virtual void doRotatesurfacepoints(Real alpha, Real beta, Real gamma) = 0;

/**
@brief Gets the index of the SH coefficient.

@todo Determine what this function does and why it is needed.

@param n is the degree of the SH coefficient
@param m is the order of the SH coefficient
@return the index
*/
int getIndex(int n, int m) { return (m + n); }

/**
@brief Calculates and returns the factorial of the parameter.

This is simply an interface to the boost library method for calculating
the factorial of a number.

@param j is the number whose factorial is desired
@return the factorial <i>j</i>!
*/
Real fac(int j)
{
    return (boost::math::factorial<Real>((unsigned)j));
}

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
this windowing is to use Lanczos sigma factors [Sloan 2002].  Given a
SH coefficient <i>a<sub>n</sub><sup>m</sup></i>, its magnitude is multiplied by
\f[
\frac{\sin \left( \frac{\pi (n - n_0}{n_{max} - n_0} \right)}{ \left( \frac{\pi (n - n_0)}{n_{max} - n_0} \right)}
\f]
where <i>n<sub>0</sub></i> is the degree in the SH expansion at which we
should start the de-ringing process.

@note Pure virtual function; must be overridden by derived classes.

@param power is the exponent (how many times to apply the sigma factor)
@param nstart is the SH degree at which to start the de-ringing.
*/
virtual void doLanczos(const int power, int nstart) = 0;

/**
@brief Applies the inverse of the Lanczos sigma factor.

If the parameters are identical to those used when `GenusZero::doLanczos`
was applied, this function undoes the effect and restores the object to
its original state.

@note Pure virtual function; must be overridden by derived classes.

@param power is the exponent (how many times to apply the inverse sigma factor)
@param nstart is the SH degree at which to start.
*/
virtual void undoLanczos(const int power, int nstart) = 0;

/**
@brief Applies the Hanning function to de-ring the object
 
An alternative to the Lanczos sigma factor is the Hanning filter, given by
\f[
\frac{1 + \cos \left( \frac{\pi (n - n_0)}{n_{max} - n_0} \right)}{2}
\f]
The Hanning filter decays faster than the Lanczos filter, but not quite as fast
as the square of the Lanczos filter.

@note Pure virtual function; must be overridden by derived classes.

@param nstart is the SH degree at which to start the de-ringing.
*/
virtual void doHanning(int nstart) = 0;

/**
@brief Applies the inverse of the Hanning function.

If the parameter is identical to that used when `GenusZero::doHanning`
was applied, this function undoes the effect and restores the object to
its original state.

@note Pure virtual function; must be overridden by derived classes.

@param nstart is the SH degree at which to start.
*/
virtual void undoHanning(int nstart) = 0;
/**@}*/

/**
@brief Gets the convexity of the object.

This functions generates a polyhedral approximation of the convex hull of
the object <i>S</i>, which is denoted as conv(<i>S</i>).  The
convex hull is generated from algorithms described in @cite ORourke98 and
implemented by the CGAL library.

Once the hull has been found, the function measures the
convexity <i>C</i> using the ratio of the volume enclosed
by <i>S</i> to the volume enclosed by conv(<i>S</i>):
\f[
	C = \frac{V(S)}{V(\mbox{conv}(S)}
\f]
This measure of <i>C</i> is probably the most common one @cite Rosin06
However, it is worthwhile noting that <i>C</i> is an asymmetric measure of convexity
in the sense that intrusions and protrusions have different influences on the
convexity; long and narrow intrusions tend to have a smaller influences on convexity
than long narrow protrusions.  A symmetric convexity measure was proposed by
Rosin and Mumford @cite Rosin06 and may be investigated at a later time.

The computation of the hull volume is done by decomposing the polyhedral representation
of the hull into disjoint polyhedra.  The algorithm is:
    - Let <i>F</i> be the set of facets on the hull, {<i>F</i><sub>1</sub>, <i>F</i><sub>2</sub>, ...}
    - \f$\forall F_i \in F$, create a tetrahedron, <i>T<sub>i</sub></i> at the origin with
      <i>F<sub>i</sub></i> as its base
    - <i>V</i>(conv(<i>S</i>)) = \f$\sum_i V(T_i)\f$

@return the convexity measure
*/
Real getConvexity(void);
    
/**
@brief Gets the resting configuration of the object.

This functions generates a polyhedral approximation of the convex hull of
the object <i>S</i>, which is denoted as conv(<i>S</i>).  The
convex hull is generated from algorithms described in @cite ORourke98 and
implemented by the CGAL library.

Once the hull has been found, its triangular facets are each interrogated to
see if they form a stable tripod on which the particle could rest under
a gravitational field.  It does this by computing the centroid, which is
not necessarily the origin of the coordinate system, and then determining
if the centroid projects within the interior of the tripod.  Finally,
for each tripod into which the base projects, the one with the smallest
distance between the centroid and the plane of the tripod is selected
as the most stable configuration.
*/
void getRestingConfiguration(void);
    
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
surface, and \f$S d\theta\, d\phi\f$ is a differential surface patch, the computation
of which will depend on the type of particle.  For this reason alone, we defer the
definition of functions for computing roundness to the actual particle class, making
them pure virtual in this base class.
*/

/**@{*/
/**
@brief Computes a roundness measure of the object

@note Pure virtual function; must be overridden by derived classes.

@param method is the type of roundness to compute.
    - 0 = Dot product of position vector with unit normal vector
    - 1 = 3D generalization of Wadell's definition @cite Wadell32
*/
virtual void computeRoundness(const unsigned int method) = 0;

/**
@brief Computes all roundness measures of the object

This function is simply an interface to the function of the same name
that computes a given measure of the roundness.  The present function
calls that function with the parameter set to compute all measures
of roundness.
*/
virtual void computeRoundness(void) { computeRoundness(ALLROUNDNESSES); }

/**
@brief Sets the total roundness of the object.

@note Not used.

@param idx is the type of roundness to set.
    - 0 = Dot product of position vector with unit normal vector
    - 1 = 3D generalization of Wadell's definition @cite Wadell32
@param ca is the roundness value
*/
void setRoundness(const unsigned int idx, const Real ca)
{
    try { cumroundness_[idx] = ca; }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","setRoundness","cumroundness_",
                            cumroundness_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Gets every measure of the total roundness of the object.

@note Not used.

@return the vector storing each roundness measure.
    - 0 = Dot product of position vector with unit normal vector
    - 1 = 3D generalization of Wadell's definition @cite Wadell32
*/
std::vector<Real> getRoundness(void) const { return cumroundness_; }

/**
@brief Gets one of the total roundness measures of the object.

@param method is the type of roundness to get.
    - 0 = Dot product of position vector with unit normal vector
    - 1 = 3D generalization of Wadell's definition @cite Wadell32
@return the roundness measure requested
*/
Real getRoundness(const unsigned int method) const
{
    try { return cumroundness_[method]; }
    catch (std::out_of_range &oor) {
         EOBException ex("Particle","setRoundness","cumroundness_",
                            cumroundness_.size(),method);
        ex.printException();
        std::cout << "That roundness measure is not yet implemented." << std::endl;
        return(0.0);
    }
}
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
@brief Sets the distance from the center to a surface mesh point.

@note Pure virtual function; must be overridden by derived class.

@param theta is the first surface parameter specifying the surface point
@param phi is the second surface parameter specifying the surface point
*/
virtual void setR(const Real theta, const Real phi) = 0;

/**
@brief Sets the distance from the center to a surface mesh point.

@note Pure virtual function; must be overridden by derived class.

@param surfaceindex is the index of the point in the vector of surface mesh points
*/
virtual void setR(const unsigned int surfaceindex) = 0;

/**
@brief Gets the distance from the center to a surface mesh point.

@param surfaceindex is the index of the surface mesh point
@return the distance to the surface mesh point
*/
virtual Real getR(const unsigned int surfaceindex);

/**
@brief Gets the distance from the center to any point on the surface.

@note Pure virtual function; must be overridden by derived class.

@param theta is the polar angle \f$\theta\f$ (radians)
@param phi is the azimuthal angle \f$\phi\f$ (radians)
@return the distance from the center to the surface in the direction
         \f$(\theta,\phi)\f$
*/
virtual Real getR(const Real theta, const Real phi) = 0;

/**
@brief Gets the distances from the center to each surface point.
@return the vector storing the distances from the center to each point
         in the surface mesh, relative to the center
*/
virtual std::vector<Real> getR(void) const { return r_; }

/**
@brief Gets the derivative \f$R_{\phi}\f$.

@note Pure virtual function; must be overridden by derived class.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\phi}\f$
*/
virtual Real getRp(Real theta, const Real phi) = 0;

/**
@brief Gets the derivative \f$R_{\phi}\f$.

@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$R_{\phi}\f$
*/
virtual Real getRp(const unsigned int surfaceindex);

/**
@brief Gets the derivative \f$R_{\theta}\f$.

@note Pure virtual function; must be overridden by derived class.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\theta}\f$
*/
virtual Real getRt(Real theta, const Real phi) = 0;

/**
@brief Gets the derivative \f$R_{\theta}\f$.

@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$R_{\theta}\f$
*/
virtual Real getRt(const unsigned int surfaceindex);

/**
@brief Gets the second derivative \f$R_{\phi \phi}\f$.

@note Pure virtual function; must be overridden by derived class.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\phi \phi}\f$
*/
virtual Real getRpp(Real theta, const Real phi) = 0;

/**
@brief Gets the second derivative \f$R_{\phi \phi}\f$.

@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$R_{\phi \phi}\f$
*/
virtual Real getRpp(const unsigned int surfaceindex);

/**
@brief Gets the second derivative \f$R_{\theta \theta}\f$.

@note Pure virtual function; must be overridden by derived class.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\theta \theta}\f$
*/
virtual Real getRtt(Real theta, const Real phi) = 0;

/**
@brief Gets the second derivative \f$R_{\theta \theta}\f$.

@param surfindex is the index of the point in the surface mesh
@return the derivative \f$R_{\theta \theta}\f$
*/
virtual Real getRtt(const unsigned int surfindex);

/**
@brief Gets the second derivative \f$R_{\theta \phi}\f$.

@note Pure virtual function; must be overridden by derived class.

@param theta is the polar angle \f$\theta\f$ at the evaluation point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the evaluation point (radians)
@return the derivative \f$R_{\theta \phi}\f$
*/
virtual Real getRtp(Real theta, const Real phi) = 0;

/**
@brief Gets the second derivative \f$R_{\theta \phi}\f$.

@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$R_{\theta \phi}\f$
*/
virtual Real getRtp(const unsigned int surfaceindex);
/**@}*/

/**
@brief Normalizing prefactor for SH expansions.

For any term in a SH expansion, there is an associated multiplicative prefactor on
its magnitude, given by
\f[
f_{nm} = sqrt{\frac{2n + 1}{4 \pi}\frac{(n-m)!}{(n+m)!}}
\f]
where <i>n</i> is the degree and <i>m</i> is the order.

@param n is the degree of the SH  term
@param m is the order of the SH  term
@return the normalizing prefactor
*/
Real getFnm(const Real n, const Real m);

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

@note Pure virtual function; must be overridden by derived classes.

@param index specifies which coordinate to get
    - 0 is the <i>x</i>-coordinate
    - 1 is the <i>y</i>-coordinate
    - 2 is the <i>z</i>-coordinate
@param theta is the polar angle \f$\theta\f$ at the surface point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the surface point (radians)
@return the requested vector component
*/
virtual Real getX(const unsigned int index, const Real theta, const Real phi) = 0;

/**
@brief Gets the vector \f$\vec{X}(\theta,\phi)\f$.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle \f$\theta\f$ at the surface point (radians)
@param phi is the azimuthal angle \f$\phi\f$ at the surface point (radians)
@return the requested vector \f$\vec{X}(\theta,\phi)\f$
*/
virtual std::vector<Real> getX(const Real theta, const Real phi) = 0;

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
@brief Gets the <i>y</i>-coordinate of every surface point
@return the vector storing the <i>y</i>-coordinate of each surface point
         in the surface mesh, relative to the center
*/
virtual std::vector<Real> getY(void) const { return y_; }

/**
@brief Gets the <i>z</i>-coordinate of every surface point
@return the vector storing the <i>z</i>-coordinate of each surface point
         in the surface mesh, relative to the center
*/
virtual std::vector<Real> getZ(void) const { return z_; }

/**
@brief Gets the derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to azimuthal angle, \f$X_{i,\phi}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\phi}\f$
*/
virtual Real getXp(const unsigned int index, const Real theta, const Real phi) = 0;

/**
@brief Gets the derivative of \f$\vec{X}\f$ with respect to
azimuthal angle, \f$\vec{X}_{\phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\phi}\f$
*/
virtual std::vector<Real> getXp(const Real theta, const Real phi) = 0;

/**
@brief Gets the derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to azimuthal angle, \f$X_{i,\phi}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\phi}\f$
*/
virtual Real getXp(const unsigned int index, const unsigned int surfaceindex) = 0;

/**
@brief Gets the derivative of \f$\vec{X}\f$ with respect to
azimuthal angle, \f$\vec{X}_{\phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\phi}\f$
*/
virtual std::vector<Real> getXp(const unsigned int surfaceindex) = 0;

/**
@brief Gets the derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to polar angle, \f$X_{i,\theta}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\theta}\f$
*/
virtual Real getXt(const unsigned int index, const Real theta, const Real phi) = 0;

/**
@brief Gets the derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\theta}\f$
*/
virtual std::vector<Real> getXt(const Real theta, const Real phi) = 0;

/**
@brief Gets the derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to polar angle, \f$X_{i,\theta}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\theta}\f$
*/
virtual Real getXt(const unsigned int index, const unsigned int surfaceindex) = 0;

/**
@brief Gets the derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\theta}\f$
*/
virtual std::vector<Real> getXt(const unsigned int surfaceindex) = 0;

/**
@brief Gets the second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to azimuthal angle, \f$X_{i,\phi \phi}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\phi \phi}\f$
*/
virtual Real getXpp(const unsigned int index, const Real theta, const Real phi) = 0;

/**
@brief Gets the second derivative of \f$\vec{X}\f$ with respect to
azimuthal angle, \f$\vec{X}_{\phi \phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\phi \phi}\f$
*/
virtual std::vector<Real> getXpp(const Real theta, const Real phi) = 0;

/**
@brief Gets the second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to azimuthal angle, \f$X_{i,\phi \phi}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\phi \phi}\f$
*/
virtual Real getXpp(const unsigned int index, const unsigned int surfaceindex) = 0;

/**
@brief Gets the second derivative of \f$\vec{X}\f$ with respect to
azimuthal angle, \f$\vec{X}_{\phi \phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\phi \phi}\f$
*/
virtual std::vector<Real> getXpp(const unsigned int surfaceindex) = 0;

/**
@brief Gets the second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to polar angle, \f$X_{i,\theta \theta}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\theta \theta}\f$
*/
virtual Real getXtt(const unsigned int index, const Real theta, const Real phi) = 0;

/**
@brief Gets the second derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta \theta}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\theta \theta}\f$
*/
virtual std::vector<Real> getXtt(const Real theta, const Real phi) = 0;

/**
@brief Gets the second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$ with respect to polar angle, \f$X_{i,\theta \theta}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\theta \theta}\f$
*/
virtual Real getXtt(const unsigned int index, const unsigned int surfaceindex) = 0;

/**
@brief Gets the second derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta \theta}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\theta \theta}\f$
*/
virtual std::vector<Real> getXtt(const unsigned int surfaceindex) = 0;

/**
@brief Gets the mixed second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$, \f$X_{i,\theta \phi}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the derivative \f$X_{i,\theta \phi}\f$
*/
virtual Real getXtp(const unsigned int index, const Real theta, const Real phi) = 0;

/**
@brief Gets the mixed second derivative of \f$\vec{X}\f$,
\f$\vec{X}_{\theta \phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle, \f$\theta\f$, at the surface point (radians)
@param phi is the azimuthal angle, \f$\phi\f$, at the surface point (radians)
@return the vector of component derivatives \f$\vec{X}_{\theta \phi}\f$
*/
virtual std::vector<Real> getXtp(const Real theta, const Real phi) = 0;

/**
@brief Gets the mixed second derivative of the <i>i</i>-the component
of \f$\vec{X}\f$, \f$X_{i,\theta \phi}\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the component derivative to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param surfaceindex is the index of the point in the surface mesh
@return the derivative \f$X_{i,\theta \phi}\f$
*/
virtual Real getXtp(const unsigned int index, const unsigned int surfaceindex) = 0;

/**
@brief Gets the mixed second derivative of \f$\vec{X}\f$ with respect to
polar angle, \f$\vec{X}_{\theta \phi}\f$.

Each vector component is differentiated, and the result is returned
in a vector.

@note Pure virtual function; must be overridden by derived classes.

@param surfaceindex is the index of the point in the surface mesh
@return the vector of component derivatives \f$\vec{X}_{\theta \phi}\f$
*/
virtual std::vector<Real> getXtp(const unsigned int surfaceindex) = 0;
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
       in the surface mesh.

@param si is the index of the point in the surface mesh
@return the first fundamental form component <i>E</i>
*/
virtual Real getE(const unsigned int si);

/**
@brief Gets the first fundamental form <i>E</i> component of a point
       on the surface with angular coordinates \f$(\theta,\phi)\f$.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the first fundamental form component <i>E</i>
*/
virtual Real getE(const Real theta, const Real phi) = 0;

/**
@brief Gets the first fundamental form <i>F</i> component of a point
       in the surface mesh.

@param si is the index of the point in the surface mesh
@return the first fundamental form component <i>F</i>
*/
virtual Real getF(const unsigned int si);

/**
@brief Gets the first fundamental form <i>F</i> component of a point
       on the surface with angular coordinates \f$(\theta,\phi)\f$.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the first fundamental form component <i>F</i>
*/
virtual Real getF(const Real theta, const Real phi) = 0;

/**
@brief Gets the first fundamental form <i>G</i> component of a point
       in the surface mesh.

@param si is the index of the point in the surface mesh
@return the first fundamental form component <i>G</i>
*/
virtual Real getG(const unsigned int si);

/**
@brief Gets the first fundamental form <i>G</i> component of a point
       on the surface with angular coordinates \f$(\theta,\phi)\f$.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle (radians)
@param phi is the azimuthal angle (radians)
@return the first fundamental form component <i>G</i>
*/
virtual Real getG(const Real theta, const Real phi) = 0;
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

/**@{*/
/**
@brief Gets the second fundamental form <i>L</i> component of a point
       in the surface mesh.

@param si is the index of the point in the surface mesh
@return the second fundamental form component <i>L</i>
*/
virtual Real getL(const unsigned int si);

/**
@brief Gets the second fundamental form <i>L</i> component of a point
       on the surface with angular coordinates \f$(\theta,\phi)\f$.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return the second fundamental form component <i>L</i>
*/
virtual Real getL(const Real theta, const Real phi);

/**
@brief Gets the second fundamental form <i>M</i> component of a point
       in the surface mesh.

@param si is the index of the point in the surface mesh
@return the second fundamental form component <i>M</i>
*/
virtual Real getM(const unsigned int si);

/**
@brief Gets the second fundamental form <i>M</i> component of a point
       on the surface with angular coordinates \f$(\theta,\phi)\f$.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return the second fundamental form component <i>M</i>
*/
virtual Real getM(const Real theta, const Real phi);

/**
@brief Gets the second fundamental form <i>N</i> component of a point
       in the surface mesh.

@param si is the index of the point in the surface mesh
@return the second fundamental form component <i>N</i>
*/
virtual Real getN(const unsigned int si);

/**
@brief Gets the second fundamental form <i>N</i> component of a point
       on the surface with angular coordinates \f$(\theta,\phi)\f$.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return the second fundamental form component <i>N</i>
*/
virtual Real getN(const Real theta, const Real phi);
/**@}*/

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

@note Pure virtual function; must be overridden by derived classes.

@param index is the vector component to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param si is the index of the point in the surface mesh
@return the requested component of the normal vector
*/
virtual Real getNormal(const unsigned int index, const unsigned int si) = 0;

/**
@brief Gets the <i>i</i>-th component of the normal vector at a point
       given by angular coordinates \f$(\theta,\phi)\f$.

@note Pure virtual function; must be overridden by derived classes.

@param index is the vector component to return
    - 0 = <i>x</i> component
    - 1 = <i>y</i> component
    - 2 = <i>z</i> component
@param theta is the polar angle of the point (radians)
@param phi is the azimuthal angle of the point (radians)
@return the requested component of the normal vector
*/
virtual Real getNormal(const unsigned int index, const Real theta, const Real phi) = 0;

/**
@brief Gets the normal vector at a point
       given by its index in the surface mesh.

@param si is the index of the point in the surface mesh
@return the normal vector at the point
*/
virtual std::vector<Real> getNormal(const unsigned int si);

/**
@brief Gets the normal vector at a point
       given by angular coordinates \f$(\theta,\phi)\f$.

@note Pure virtual function; must be overridden by derived classes.

@param theta is the polar angle of the point (radians)
@param phi is the azimuthal angle of the point (radians)
@return the normal vector at the point
*/
virtual std::vector<Real> getNormal(const Real theta, const Real phi) = 0;
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
@brief Gets the Gaussian curvature at a pre-defined point in the surface mesh.

@param si is the index of the point in the surface mesh
@return the Gaussian curvature
*/
Real getK(const unsigned int si);

/**
@brief Gets the Gaussian curvature at an arbitrary surface point.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return the Gaussian curvature
*/
Real getK(const Real theta, const Real phi);

/**
@brief Gets the mean curvature at a pre-defined point in the surface mesh.

@param si is the index of the point in the surface mesh
@return the mean curvature
*/
Real getH(const unsigned int si);

/**
@brief Gets the mean curvature at an arbitrary surface point.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return the mean curvature
*/
Real getH(const Real theta, const Real phi);

/**
@brief Gets one of the two principal directions at a pre-defined point in the surface mesh.

@param index is 0 for the first principal direction, 1 for the second
@param si is the index of the point in the surface mesh
@return the vector components of the principal direction
*/
std::vector<Real> getPrincipaldirection(const unsigned int index, const unsigned int si);

/**
@brief Gets one of the two principal directions at an arbitrary surface point.

@param index is 0 for the first principal direction, 1 for the second
@param theta is the first surface parameter
@param phi is the second surface parameter
@return the vector components of the principal direction
*/
std::vector<Real> getPrincipaldirection(const int unsigned index, const Real theta, const Real phi);

/**
@brief Gets one of the two principal curvatures at a pre-defined point in the surface mesh.

@param index is 0 for the first principal curvature, 1 for the second
@param si is the index of the point in the surface mesh
@return the principal curvature
*/
Real getPrincipalcurvature(const unsigned int index, const unsigned int si);

/**
@brief Gets one of the two principal curvatures at an arbitrary surface point.

@param index is 0 for the first principal curvature, 1 for the second
@param theta is the first surface parameter
@param phi is the second surface parameter
@return the principal curvature
*/
Real getPrincipalcurvature(const unsigned int index, const Real theta, const Real phi);

/**
@brief Gets one of the two principal curvatures at an arbitrary surface point.

@param index is 0 for the first principal direction, 1 for the second
@param theta is the first surface parameter
@param phi is the second surface parameter
@return the principal curvature
*/
std::vector<Real> getPrincipalcurvatures(const unsigned int si);

/**
@brief Computes both principal curvatures at an arbitrary surface point.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return a vector storing both principal curvatures
*/
std::vector<Real> computePrincipalcurvatures(const Real theta, const Real phi);

/**
@brief Computes one of the principal curvatures at an arbitrary surface point.

@param index is 0 for the first principal direction, 1 for the second
@param theta is the first surface parameter
@param phi is the second surface parameter
@return the principal curvature
*/
Real computePrincipalcurvature(const unsigned int index, const Real theta, const Real phi);

/**
@brief Computes (does not return) both principal curvatures at every point
       in the pre-defined surface mesh.
*/
void computePrincipalcurvatures(void);

/**
@brief Gets one of the principal curvatures at an elliptic umbilic point.

@param index is 0 for the first principal curvature, 1 for the second
@param si is the index of the point in the surface mesh
@return the principal curvature
*/
Real getEllipticumbiliccurvature(const unsigned int index, const unsigned int si);

/**
@brief Gets one of the principal curvatures at an elliptic umbilic point.

@param index is 0 for the first principal curvature, 1 for the second
@param theta is the first surface parameter for the point
@param phi is the second surface parameter for the point
@return the principal curvature
*/
Real getEllipticumbiliccurvature(const unsigned int index,
                                 const Real theta, const Real phi);

/**
@brief Is the given point an umbilic point?

@param si is the index of the point in the surface mesh
@return true iff the point is an umbilic point
*/
bool isUmbilic(const unsigned int si);

/**
@brief Is the given point an umbilic point?

@param theta is the first surface parameter
@param phi is the second surface parameter
@return true iff the point is an umbilic point
*/
bool isUmbilic(const Real theta, const Real phi);

/**
@brief Is the given point an elliptic umbilic point?

@param si is the index of the point in the surface mesh
@return true iff the point is an elliptic umbilic point
*/
bool isEllipticumbilic(const unsigned int si);

/**
@brief Is the given point an elliptic umbilic point?

@param theta is the first surface parameter
@param phi is the second surface parameter
@return true iff the point is an elliptic umbilic point
*/
bool isEllipticumbilic(const Real theta, const Real phi);

/**
@brief Is the given point a parabolic umbilic point?

@param si is the index of the point in the surface mesh
@return true iff the point is a parabolic umbilic point
*/
bool isParabolicumbilic(const unsigned int si);

/**
@brief Is the given point a parabolic umbilic point?

@param theta is the first surface parameter
@param phi is the second surface parameter
@return true iff the point is a parabolic umbilic point
*/
bool isParabolicumbilic(const Real theta, const Real phi);

/**
@brief Gets the minimum of the two principal curvatures at a pre-defined
       surface point

@param si is the index of the point in the surface mesh
@return the minimum principal curvature
*/
Real getMincurvature(const unsigned int si);

/**
@brief Gets the minimum of the principal curvatures at an
       arbitrary surface point.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return the minimum principal curvature
*/
Real getMincurvature(const Real theta, const Real phi);

/**
@brief Gets the maximum of the two principal curvatures at a pre-defined
       surface point

@param si is the index of the point in the surface mesh
@return the maximum principal curvature
*/
Real getMaxcurvature(const unsigned int si);

/**
@brief Gets the maximum of the principal curvatures at an
       arbitrary surface point.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return the maximum principal curvature
*/
Real getMaxcurvature(const Real theta, const Real phi);

/**
@brief Gets the minimum magnitude of the two principal curvatures
       at a pre-defined surface point

@param si is the index of the point in the surface mesh
@return the minimum magnitude of the principal curvatures
*/
Real getMinabscurvature(const unsigned int si);

/**
@brief Gets the minimum magnitude of the principal curvatures at an
       arbitrary surface point.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return the minimum magnitude of the principal curvatures
*/
Real getMinabscurvature(const Real theta, const Real phi);

/**
@brief Gets the maximum magnitude of the two principal curvatures
       at a pre-defined surface point

@param si is the index of the point in the surface mesh
@return the maximum magnitude of the principal curvatures
*/
Real getMaxabscurvature(const unsigned int si);

/**
@brief Gets the maximum magnitude of the principal curvatures at an
       arbitrary surface point.

@param theta is the first surface parameter
@param phi is the second surface parameter
@return the maximum magnitude of the principal curvatures
*/
Real getMaxabscurvature(const Real theta, const Real phi);

/**
@brief Computes and returns the mean curvature integrated over the surface.

@Pure virtual function; must be overridden by derived classes.

@return the mean curvature, <i>H</i>, integrated over the surface
*/
virtual Real getIntegratedH(void) = 0;

/**
@brief Computes and returns the Gaussian curvature integrated over the surface.

@Pure virtual function; must be overridden by derived classes.

@return the Gaussian curvature, <i>K</i>, integrated over the surface
*/
virtual Real getIntegratedK(void) = 0;
/**@}*/

/**
@brief Computes and returns the moment of inertia tensor

@Pure virtual function; must be overridden by derived classes.

@param itensor is the vector of upper diagonal elements of moment of inertia tensor
@param rg2itensor is the vector of upper diagonal elements of gyration tensor
@param isphere is the moment of inertia of sphere with same volume as particle
@param rg2sphere is the square of the radius of gyratio of sphere with same volume as particle

@return 0 upon successful completion
*/
virtual int getI(std::vector<Real> &itensor, std::vector<Real> &rg2tensor, Real &isphere, Real &rg2sphere) = 0;

/**
@name Sphericity
@brief Functions for calculating and returning various measures of the object's sphericity.

Most of these functions rely on the triaxial dimensions of the object to compute
various measures of the object's anisometry.  As such, these functions are generic
to any class of particle and can therefore be defined here.
*/

/**@{*/
/**
@brief Computes the "true" sphericity as defined by Wadell @cite Wadell32.

Wadell defined the true sphericity of an object by comparing its surface
area to that of a sphere having the same volume:
\f[
    S_W = \frac{A_{\text{sphere}}}{A}
\f]

@return the true (Wadell) sphericity, <i>S<sub>W</sub></i>
*/
Real computeSphericity(void) const
{
    if (area_ > 0.0) {
        return ((4.0 * Pi * pow((3.0 * volume_ / (4.0 * Pi)),(2.0/3.0))) / area_);
    } else {
        throw FloatException("Particle","computeSphericity",
                             "Divide by zero area_");
    }
}

/**
@brief Computes the "true" sphericity as defined by Wadell @cite Wadell32.

This is just an interface to the function `Particle::computeSphericity()`.

@return the true (Wadell) sphericity, <i>S<sub>W</sub></i>
*/
Real getWadellsphericity(void) const
{
    try {
        return (computeSphericity());
    }
    catch (FloatException flex) { throw flex; }
}

/**
@brief Computes the Corey shape factor (CSF).

The Corey shape factor @cite Corey49 @cite LeRoux97 is
\f[
\text{CSF} = D_3 / \sqrt{D_1 D_2}
\]

@return the Corey shape factor
*/
Real getCSF(void) const
{
   if (triaxialdim_[0] * triaxialdim_[1] > 0.0) {
       return (triaxialdim_[2] / (sqrt(triaxialdim_[0] * triaxialdim_[1])));
   } else {
       throw FloatException("Particle","getCSF","CSF denominator is zero");
   }
}

/**
@brief Computes the Aschenbrenner shape factor (ASF).

The Aschenbrenner shape factor @cite Aschenbrenner56 @cite LeRoux97 is
\f[
\text{ASF} = \frac{13.392 \left( \frac{D_3^2}{D_1 D_2} \right)^{1/3}}{1 + \frac{D_3}{D_2}
+ \frac{D_3}{D_1} + 6 \sqrt{1 + (D_3/D_2)^2 + (D_3/D_1)^2}
\]

@return the Aschenbrenner shape factor
*/
Real getASF(void);

/**
@brief Computes the maximum projection sphericity (MPS).

The maximum projected sphericity, defined by Sneed and Folk @cite Sneed58 and summarized
in @cite LeRoux97
\f[
\text{MPS} = \left( \frac{D_3^2}{D_1 D_2} \right)^{1/3}
\]

@return the maximum projection sphericity
*/
Real getMPS(void) const
{
    
    Real d1,d2,d3;
    d1 = triaxialdim_[0];
    d2 = triaxialdim_[1];
    d3 = triaxialdim_[2];
    if (d1 > 0.0 && d2 > 0.0) {
        return (pow((d3 * d3 / d1 / d2),(1.0/3.0)));
    } else {
        throw FloatException("Particle","getMPS","MPS denominator is zero");
   }
}

/**
@brief Computes the E shape factor.

The E shape factor, defined by Janke @cite Janke66 and summarized
in @cite LeRoux97 is
\f[
\text{ESF} = D_3 \sqrt{\frac{3}{D_1^2 + D_2^2 + D_3^2}}
\]

@return the E shape factor
*/
Real getESF(void) const
{
    Real d1,d2,d3;
    d1 = triaxialdim_[0];
    d2 = triaxialdim_[1];
    d3 = triaxialdim_[2];
    Real denom = (d3 * d3) + (d2 * d2) + (d1 * d1);
    if (denom > 0.0) {
        return (d3 * sqrt(3.0/denom));
    } else {
        throw FloatException("Particle","getESF","ESF denominator is zero");
   }
}

/**
@brief Computes the Hofmann shape entropy.

The shape entropy, defined by Hofmann @cite Hofmann94 and summarized
in @cite LeRoux97 is
\f[
\S_H = \frac{p_1 \ln p_1 + p_2 \ln p_2 + p_3 \ln p_3}{\ln (1/3)}
\]
where, for example,
\f[
p_1 = \frac{D_1}{D_1 + D_2 + D_3}
\]

@return the Hofmann shape entropy
*/
Real getShapeentropy(void);
/**@}*/
    
/**
@name Output
@brief Functions for producing visualizations or formatted output of properties.
*/

/**@{*/
/**
@brief Outputs a VRML file for viewing the object.

@note Pure virtual function; must be overridden by derived classes.

@param vrmlname is the name of the VRML file to write
@param render_inscribed_sphere is true iff the inscribed sphere should be drawn
@param render_enclosing_sphere is true iff the enclosing sphere should be drawn
@param render_hull is true iff the convex hull should be drawn
@param surface_color_scheme specifies how to color the surface (NOT USED)
@param rotation is the rotation angle about the y axis in degrees
@return 0 if file created successfully, non-zero otherwise
*/
virtual int createVRML(const std::string vrmlname, bool render_inscribed_sphere,
             bool render_enclosing_sphere, bool render_hull,
             int surface_color_scheme, Real rotation) = 0;

/**
@brief Outputs the SH coefficients for the particle in its current state

@param shname is the name of the file to write
@return 0 if successful
*/
virtual int createSHFile(const std::string shname) = 0;

/**
@brief Makes a digitized representation of the particle in a bounding box

@param digitizedname is the name of the image file to write
@param resolution is the voxel edge length relative to the unit of measurement
@return 0 if successful
*/
virtual int createDigitizedImage(const std::string digitizedname, Real resolution) = 0;

/**
@brief Produce formatted output of the particle properties.
*/
virtual void printProperties();
/**@}*/

/**
@brief Get a string describing the particle class.

@return reference to string describing the class type
*/
virtual std::string &getType() const { return (std::string &)("GENERIC PARTICLE"); }

/**
@name Simple operations on 3D vectors
@brief A collection of functions for doing vector operations
*/

/**@{*/
/**
@brief The inner product of two STL vectors

Note that this functionality already exists in the STL using the numeric functions, but this is
a simpler implementation for two 3D vectors.  All bounds need to be right or else this won't work.

@todo Make this more robust by doing bounds checking, type checking, etc.  Perhaps even make this
a template function.
*/
Real getInnerProduct(std::vector<Real> v1, std::vector<Real> v2)
{
    return ((v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]));
}

/**
@brief The Cartesian product of two STL vectors

@todo Make this more robust by doing bounds checking, type checking, etc.  Perhaps even make this
a template function.
*/
std::vector<Real> getCartesianProduct(std::vector<Real> u, std::vector<Real> v)
{
    std::vector<Real> result;
    result.clear();
    result.resize(3,0.0);
    result[0] = u[1]*v[2] - u[2]*v[1];
    result[1] = u[3]*v[0] - u[0]*v[3];
    result[2] = u[0]*v[1] - u[1]*v[0];
    return result;
}

/**
@brief The magnitude of a STL vector

@todo Make this more robust by doing bounds checking, type checking, etc.  Perhaps even make this
a template function.
*/
Real getMagnitude(std::vector<Real> u)
{
    return (std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
}

/**@}*/

/**
@name Operators.
@brief Overloaded operators for the Particle class
*/

/**@{*/
/**
@brief Overloads the assignment operator.

@param &that is a reference to the object on the right side of the assignment
@return a reference to the Particle object assigned
*/
Particle& operator=(const Particle& that);
/**@}*/

};

} // End of tsquare namespace
#endif
