/**
@file
@brief Method definitions for the Particle base class.
*/
#include "Particle.h"

tsquare::Particle::Particle()
{
    xlo_ = ylo_ = zlo_ = 0.0;
    xhi_ = yhi_ = zhi_ = 0.0;
    volume_ = area_ = 0.0;
    narea_ = 0.0;
    ngc_ = diam_ = itrace_ = 0.0;
    nmax_ = nmax_from_anmfile_ = 0;
    dim_.clear();
    dim_.resize(3,0.0);
    ndim_ = dim_;
    triaxialdim_ = dim_;
    a_.clear();
    tsquare::ExtendedVector ev;
    ev.clear();
    r_.clear();
    x_.clear();
    y_.clear();
    z_.clear();
    centroid_.clear();
    centroid_.resize(3,0.0);
    surface_.clear();
    xg_.clear();
    wg_.clear();
    hullvolume_ = hullarea_ = 0.0;
    verbose_ = false;

    inscribedsphere_.setXr(0.0);
    inscribedsphere_.setYr(0.0);
    inscribedsphere_.setZr(0.0);
    inscribedsphere_.setRadius(0.0);

    enclosingsphere_.setXr(0.0);
    enclosingsphere_.setYr(0.0);
    enclosingsphere_.setZr(0.0);
    enclosingsphere_.setRadius(0.0);
}

void tsquare::Particle::addtoSurface(Real theta, const Real phi)
{
    using std::sin;

    std::vector<Real> orderedpair;
    orderedpair.clear();
    if (sin(theta) == 0.0) theta = 0.001 * Pi;
    orderedpair.push_back(theta);
    orderedpair.push_back(phi);
    surface_.push_back(orderedpair);
    return;
}

Real tsquare::Particle::computeRefinedlength(std::vector<std::vector<Real> > sf,
                                    std::vector<std::vector<Real> > sc, std::vector<Real> xf,
                                    std::vector<Real> yf, std::vector<Real> zf,
                                    std::vector<Real> xc, std::vector<Real> yc,
                                    std::vector<Real> zc, std::vector<Real> rf,
                                    std::vector<Real> rc, std::vector<Real> &lengthvec,
                                    Real rm, Real &t1, Real &t2, Real &p1, Real &p2)
{

    using std::sqrt;
    using std::pow;

    Real x1,x2,y1,y2,z1,z2;
    Real xn,yn,zn,rcomp,rmax;

    rmax = rm;
    if ((xf.size() != rf.size()) || (xf.size() != yf.size())
         || (xf.size() != zf.size())) {
        throw tsquare::EOBException("Particle","computeRefinedlength","fine coordinates",
                           rf.size(),xf.size());
    }
    if ((xc.size() != rc.size()) || (xc.size() != yc.size())
         || (xc.size() != zc.size())) {
        throw tsquare::EOBException("Particle","computeRefinedlength","coarse coordinates",
                           rc.size(),xc.size());
    }

    for (register int i = 0; i < xf.size(); i++) {
      x1 = xf[i];
      y1 = yf[i];
      z1 = zf[i];
      for (register int j = 0; j < xc.size(); j++) {
        x2 = xc[j];
        y2 = yc[j];
        z2 = zc[j];

        rcomp = sqrt(pow(x1-x2,2.0)+pow(y1-y2,2.0)+pow(z1-z2,2.0));
    
        if (rmax < rcomp) {
          rmax = rcomp;
          t1 = sf[i][0];
          p1 = sf[i][1];
          t2 = sc[j][0];
          p2 = sc[j][1]; 

          xn = (x2 - x1) / rmax;
          yn = (y2 - y1) / rmax;
          zn = (z2 - z1) / rmax;
        }
      }
    }

    lengthvec.clear();
    lengthvec.resize(3,0.0);
    lengthvec[0] = xn;
    lengthvec[1] = yn;
    lengthvec[2] = zn;
    return rmax;
}

Real tsquare::Particle::computeRefinedwidth(std::vector<std::vector<Real> > sf,
                                   std::vector<std::vector<Real> > sc, std::vector<Real> xf,
                                   std::vector<Real> yf, std::vector<Real> zf,
                                   std::vector<Real> xc, std::vector<Real> yc,
                                   std::vector<Real> zc, std::vector<Real> rf,
                                   std::vector<Real> rc, std::vector<Real> lvec,
                                   std::vector<Real> &wvec, Real rm, Real ddot,
                                   Real &tt1, Real &tt2, Real &pp1, Real &pp2)
{
    using std::sqrt;
    using std::pow;
    using std::fabs;

    Real theta1,theta2,phi1,phi2;
    Real x1,x2,y1,y2,z1,z2;
    Real wx,wy,wz,wwx,wwy,wwz;
    Real t1,t2,p1,p2;
    Real lx,ly,lz,dot;
    Real rmax;

    rmax = rm;
    lx = lvec[0];
    ly = lvec[1];
    lz = lvec[2];
    wwx = wvec[0];
    wwy = wvec[1];
    wwz = wvec[2];

    t1 = tt1;
    t2 = tt2;
    p1 = pp1;
    p2 = pp2;

    if ((rf.size() != xf.size()) || (rf.size() != yf.size())
         || (rf.size() != zf.size())) {
        throw tsquare::EOBException("Particle","computeRefinedlength","fine coordinates",
                           rf.size(),xf.size());
    }
    if ((rc.size() != xc.size()) || (rc.size() != yc.size())
         || (rc.size() != zc.size())) {
        throw tsquare::EOBException("Particle","computeRefinedlength","coarse coordinates",
                           rc.size(),xc.size());
    }
    for (register int i = 0; i < rf.size(); i++) {
      theta1 = sf[i][0];
      phi1 = sf[i][1];
      x1 = xf[i];
      y1 = yf[i];
      z1 = zf[i];
      for (register int j = 0; j < rc.size(); j++) {
        theta2 = sc[j][0];
        phi2 = sc[j][1];
        x2 = xc[j];
        y2 = yc[j];
        z2 = zc[j];

        rm = sqrt(pow(x1-x2,2.0)+pow(y1-y2,2.0)+pow(z1-z2,2.0));

        if (rm > 0.0) {
          wx = (x2 - x1) / rm;
          wy = (y2 - y1) / rm;
          wz = (z2 - z1) / rm;

          dot = (wx * lx) + (wy * ly) + (wz * lz);

          if ((rmax < rm) && (fabs(dot) <= fabs(ddot))) {
              rmax = rm;
              t1 = theta1;
              p1 = phi1;
              t2 = theta2;
              p2 = phi2;
              wwx = wx;
              wwy = wy;
              wwz = wz;
          }
        }
      }
    }

    tt1 = t1;
    tt2 = t2;
    pp1 = p1;
    pp2 = p2;

    wvec[0] = wwx;
    wvec[1] = wwy;
    wvec[2] = wwz;

    return rmax;
}

Real tsquare::Particle::computeRefinedthickness(std::vector<std::vector<Real> > sf,
                                       std::vector<std::vector<Real> > sc, std::vector<Real> xf,
                                       std::vector<Real> yf, std::vector<Real> zf,
                                       std::vector<Real> xc, std::vector<Real> yc,
                                       std::vector<Real> zc, std::vector<Real> rf,
                                       std::vector<Real> rc, std::vector<Real> lvec,
                                       std::vector<Real> wvec, std::vector<Real> &tvec,
                                       Real rm, Real ddotl, Real ddotw, Real &tt1,
                                       Real &tt2, Real &ppp1, Real &ppp2)
{
    using std::cout;
    using std::endl;
    using std::pow;
    using std::sqrt;
    using std::fabs;

    Real theta1,theta2,phi1,phi2;
    Real x1,y1,z1,x2,y2,z2;
    Real tx,ty,tz,ttx,tty,ttz;
    Real lx,ly,lz,wx,wy,wz;
    Real dotl,dotw;
    Real rmax = rm;

    lx = lvec[0];
    ly = lvec[1];
    lz = lvec[2];
    wx = wvec[0];
    wy = wvec[1];
    wz = wvec[2];
    
    // Refine thickness

    if (verbose_) {
      cout << "Computing refined thickness" << endl;
      cout.flush();
    }

    // Now do the actual pairwise search for refined thickness

    if ((rf.size() != xf.size()) || (rf.size() != yf.size())
         || (rf.size() != zf.size())) {
        throw tsquare::EOBException("Particle","computeRefinedthickness","fine coordinates",
                           rf.size(),xf.size());
    }
    if ((rc.size() != xc.size()) || (rc.size() != yc.size())
         || (rc.size() != zc.size())) {
        throw tsquare::EOBException("Particle","computeRefinedthickness","coarse coordinates",
                           rc.size(),xc.size());
    }

    for (register int i = 0; i < rf.size(); i++) {
      theta1 = sf.at(i)[0];
      phi1 = sf.at(i)[1];
      x1 = xf.at(i);
      y1 = yf.at(i);
      z1 = zf.at(i);
      for (register int j = 0; j < rc.size(); j++) {
        theta2 = sc.at(j)[0];
        phi2 = sc.at(j)[1];
        x2 = xc.at(j);
        y2 = yc.at(j);
        z2 = zc.at(j);

        rm = sqrt(pow(x1-x2,2.0)+pow(y1-y2,2.0)+pow(z1-z2,2.0));

        // Compute unit vector along trial width

        if (rm > 0.0) {
          tx = (x2 - x1) / rm;
          ty = (y2 - y1) / rm;
          tz = (z2 - z1) / rm;

          // Next, compute the dot product between this unit vector
          // and the unit vector along the width.  This must be zero
          // (i.e., perpendicular direction to length) as a condition
          // for this to be the width

          dotl = (tx * lx) + (ty * ly) + (tz * lz);
          dotw = (tx * wx) + (ty * wy) + (tz * wz);


          if ((rmax < rm) && (fabs(dotl) <= fabs(ddotl))
                          && (fabs(dotw) <= fabs(ddotw))) {
            rmax = rm;
            tt1 = theta1;
            ppp1 = phi1;
            tt2 = theta2;
            ppp2 = phi2;
            ttx = tx;
            tty = ty;
            ttz = tz;
          }
        }
      }
    }

    if (verbose_) {
        cout << "Refined thickness = " << rmax << endl;
        cout.flush();
    }

    tvec[0] = ttx;
    tvec[1] = tty;
    tvec[2] = ttz;

    return rmax;
}

void tsquare::Particle::getLegendre(int n, Real *x, Real *pn, Real *pnm1, Real *pnp1)
{
    int k;
    std::vector<Real> P;
    P.clear();
    P.push_back(boost::math::legendre_p(0,*x));
    P.push_back(boost::math::legendre_p(1,*x));
    for (k = 1; k <= n; k++) {
        P.push_back(boost::math::legendre_next(k,*x,P[k],P[k-1]));
    }

    *pn = P[n];
    *pnm1 = P[n-1];
    *pnp1 = P[n+1];
    return;
}

void tsquare::Particle::printProperties()
{
    using std::cout;
    using std::endl;

    cout << getType() << ":" << endl;
    cout << "xlow,xhi = " << getXlo() << "," << getXhi() << endl;
    cout << "ylow,yhi = " << getYlo() << "," << getYhi() << endl;
    cout << "zlow,zhi = " << getZlo() << "," << getZhi() << endl;
    cout << "volume = " << getVolume() << endl;
    cout << "surface area = " << getArea() << endl;
    cout << "normalized surface area = " << getNarea() << endl;
    cout << "diameter = " << getDiam() << endl;
    cout << "length, width, thickness = " << getLength() << ", " << getWidth() << ", " << getThickness() << endl;
    cout << "normalized length = " << getNlength() << endl;
    cout << "normalized width = " << getNwidth() << endl;
    cout << "Itrace = " << getItrace() << endl;
    cout << "Normalized Gaussian curvature = " << getNgc() << endl;
    cout << "Max number of spherical harmonics = " << getNmax() << endl;
    cout.flush();
    cout << endl << endl;
    return;
}

void tsquare::Particle::getSurfacesubset(Real tmin, Real tmax, Real pmin, Real pmax, std::vector<int> &ssub)
{

    // First just get the order of the angles correct and make sure they are in range

    Real ttmp,ptmp;
    if (tmin > tmax) {
        ttmp = tmax;
        tmax = tmin;
        tmin = ttmp;
    }
    if (tmax < 0.0) {
        tmin += Pi;
        tmax += Pi;
    }
    if (pmin > pmax) {
        ptmp = pmax;
        pmax = pmin;
        pmin = ptmp;
    }
    if (pmax < 0.0) {
        pmin += (2.0 * Pi);
        pmax += (2.0 * Pi);
    }

    // Next loop over all surface points and add to vector if they fall within range
    // The special case of theta < 0 is handled by taking all phi values for -theta
    // The ssub vector just holds the indices of all the eligible surface points

    ssub.clear();

    for (int i = 0; i < surface_.size(); i++) {
        ttmp = surface_[i][0];
        ptmp = surface_[i][1];
        if ((tmin < 0.0) && (ttmp < -tmin)) {
            ssub.push_back(i);
        } else if ((tmax > Pi) && (ttmp > (2.0 * Pi - tmax))) {
            ssub.push_back(i);
        } else if ((ttmp >= tmin) && (ttmp <= tmax) && (ptmp >= pmin) && (ptmp <= pmax)) {
            ssub.push_back(i);
        }
    }

    return;
}

Polyhedron tsquare::Particle::getSurfacePolyhedronFromImage(std::string &fname, std::vector<Real> centroid,
                                                   Real xBB, Real yBB, Real zBB)
{
    using std::cout;
    using std::endl;

    Polyhedron surface_poly;

    cout << "Using a surface mesh on gray level image" << endl;
    cout.flush();

    Tr tr;  // 3D Delaunay triangulation
    C2t3 c2t3 (tr);  // 2D-complex in 3D-Delaunay triangulation

    cout << "Going into binary2Gray function to create gray level image" << endl;
    cout.flush();

    Gray_level_image image = binary2Gray(fname,0.5);
    
    // Now we must carefully choose the bounding sphere.  The center must be inside
    // the surface defined by image and the radius must be large enough so that the
    // sphere actually bounds the whole image

    cout << "Creating bounding sphere" << endl;
    cout.flush();

    GT::Point_3 bounding_sphere_center(centroid[0],centroid[1],centroid[2]);
    GT::FT bounding_sphere_squared_radius = 2.0 * ((xBB*xBB) + (yBB*yBB) + (zBB*zBB));
    GT::Sphere_3 bounding_sphere(bounding_sphere_center,bounding_sphere_squared_radius);
    
    // Definition of the surface, with 10^-5 as relative precision

    cout << "Creating implicit surface from image" << endl;
    cout.flush();

    Surface_3 surface(image,bounding_sphere,1.0e-5);
    
    // To make the mesh, we first define the meshing criteria with the function
    // criteria below.  The arguments of this function are:
    //     - lower bound on the interior angle of a mesh triangular facet, in degrees
    //     - upper bound on the radii of surface Delaunay balls, in voxel units
    //     - upper bound on distance between the circumcenter of a mesh facet and the
    //         center of a surface Delaunay ball on this facet
    //
    // Generally, lower values of the last two parameters will produce finer meshes,
    // and higher values of the first argument (up to 30 degrees) will produce more
    // regular triangular facets

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.0,0.1,0.1);

    // Meshing surface, with the "manifold without boundary" algorithm

    cout << "Making surface mesh" << endl;
    cout.flush();

    // This function will make the surface triangulation stored in c2t3
    CGAL::make_surface_mesh(c2t3,surface,criteria,CGAL::Manifold_tag());
    
    // We need to output c2t3 as a 3D polyhedron

    cout << "Converting surface mesh to polyhedron" << endl;
    cout.flush();

    // The polyhedron is stored in the surface_poly object, which is
    // created by the next function

    CGAL:: output_surface_facets_to_polyhedron(c2t3,surface_poly);

    return surface_poly;
}

Real tsquare::Particle::getConvexity(void)
{
    using std::cout;
    using std::endl;

    try {
      register int i;
      std::vector<Point_3> points;
      std::vector<Real> coords;
      Point_3 p(0.0,0.0,0.0);
        
      cout << "Computing particle volume now... " << endl;
      cout.flush();
      computeVolume();
      cout << "Done!" << endl;
        
      // We link to the CGAL library for these computational geometry calculations.
      // The hull itself is stored as a Polyhedron_3 object (typedef'd as Polyhedron here).
      // The CGAL library has consistent conventions for specifying the facets, edges, and
      // vertices of a polyhedron, so we just need to follow those conventions.  The first
      // thing to do is create a vector of all the points defining the surface.

      cout << "Preparing vertices for CGAL convex hull computation..." << endl;
      for (i = 0; i < r_.size(); i++) {
        coords = getX(i);
        p = Point_3(coords[XDIM],coords[YDIM],coords[ZDIM]);
        points.push_back(p);
      }

      if (verbose_) {
        cout << "Done preparing vector of points for CGAL convex hull computation... " << endl;
      }
        
      // Calculate and retrieve the convex hull object, which will be stored in the convexhull_ member

      CGAL::convex_hull_3(points.begin(),points.end(),convexhull_);
        
      // We can now extract the facets and vertices using standard calls to the hull's methods.
      // Each orientable triangular facet is contained in a Facet object belonging to the
      // Polyhedron object.  The polyhedron stores all the facets in an order that can be accessed
      // through a facet iterator.  This can be done whenever we want to create a visualization
      // of the convex hull, too.

      // Next, we determine the volume of the convex hull.  We do this in a roundabout way:
      //    - Let F be the set of facets on the hull, (F_1, F_2, etc.)
      //    - For each facet, create a tetrahedron with vertex at the origin and the facet as its base.
      //    - Sum the volumes of all the tetrahedra.

      Facet_iterator fi;
      Halfedge_around_facet_circulator hc;

      CGAL::set_ascii_mode(cout);

      Point_3 O(0.0,0.0,0.0);  // This is the origin of the spherical polar coordinate system
                               // Guaranteed to be inside the object (and the hull)
      int facecount = 1;
      int edgecount = 1;
      std::vector<Point_3> tpoints;
      tpoints.resize(4);
      tpoints[0] = O;
      hullvolume_ = 0.0;
      for (fi = convexhull_.facets_begin(); fi != convexhull_.facets_end(); ++fi) {
        hc = fi->facet_begin();
        if (CGAL::circulator_size(hc) == 3) {
          edgecount = 1;
          do {
            tpoints[edgecount] = hc->vertex()->point();
            edgecount++;
          } while ( ++hc != fi->facet_begin() );
          Tetrahedron tet(tpoints[0],tpoints[1],tpoints[2],tpoints[3]);
          // cout << "Tet " << facecount << ": (" << tpoints[0].x() << "," << tpoints[0].y() << "," << tpoints[0].z() << ") ("
          //                                      << tpoints[1].x() << "," << tpoints[1].y() << "," << tpoints[1].z() << ") ("
          //                                      << tpoints[2].x() << "," << tpoints[2].y() << "," << tpoints[2].z() << ") ("
          //                                      << tpoints[3].x() << "," << tpoints[3].y() << "," << tpoints[3].z() << ")" << endl;
          hullvolume_ += (double)tet.volume();
        } else {
          cout << "  Dangit!  Facet has too many edges!" << endl;
        }
        facecount++;
      }

      cout << endl;
      cout << "Convex Hull Volume = " << hullvolume_ << endl;
      cout << endl;

      // Having calculated the hull volume, we now compute the convexity measure and report
      // it to stdout

      if (hullvolume_ > 0.0) {
        convexity_ = volume_ / hullvolume_;
      } else {
        cout << "WARNING:  Calculated convex hull volume of zero!" << endl;
        convexity_ = 0.0;
      }
        
      cout << endl << "FINAL RESULTS OF CONVEXITY CALCULATION:" << endl;
      cout << "    Particle volume = " << volume_ << endl;
      cout << "        Hull volume = " << hullvolume_ << endl;
      cout << "          Convexity = " << convexity_ << endl << endl;
        
      return convexity_;
    }

    catch (tsquare::FileException fex) { throw fex; }
    catch (tsquare::FloatException flex) { throw flex; }
    catch (tsquare::EOBException ex) { throw ex; }
    catch (tsquare::DataException dex) { throw dex; }
}

void tsquare::Particle::getRestingConfiguration(void)
{
    using std::cout;
    using std::endl;
    using std::acos;
    using std::atan2;

    try {
      Real x1,x2,x3,y1,y2,y3,z1,z2,z3;
      int count;
      register int i;
      std::vector<Point_3> points;
      std::vector<Real> coords;
      Point_3 p(0.0,0.0,0.0);
        
      // We link to the CGAL library for these computational geometry calculations.
      // The hull itself is stored as a Polyhedron_3 object (typedef'd as Polyhedron here).
      // The CGAL library has consistent conventions for specifying the facets, edges, and
      // vertices of a polyhedron, so we just need to follow those conventions.  The first
      // thing to do is create a vector of all the points defining the surface.

      cout << "Preparing vertices for CGAL convex hull computation..." << endl;
      for (i = 0; i < r_.size(); i++) {
        coords = getX(i);
        if (i == 2) {
            cout << "Point 2 coords = (" << coords[XDIM] << "," << coords[YDIM] << "," << coords[ZDIM] << ")" << endl;
            cout.flush();
        }
        p = Point_3(coords[XDIM],coords[YDIM],coords[ZDIM]);
        points.push_back(p);
      }

      if (verbose_) {
        cout << "Done preparing vector of points for CGAL convex hull computation... " << endl;
      }
        
      // Calculate and retrieve the convex hull object, which will be stored in the convexhull_ member

      CGAL::convex_hull_3(points.begin(),points.end(),convexhull_);
        
      // Compute the particle centroid

      computeCentroid();
      cout << "Particle centroid is (" << centroid_[0] << "," << centroid_[1] << "," <<  centroid_[2] << ")" << endl;
      cout.flush();

      // We can now extract the facets and vertices using standard calls to the hull's methods.
      // Each orientable triangular facet is contained in a Facet object belonging to the
      // Polyhedron object.  The polyhedron stores all the facets in an order that can be accessed
      // through a facet iterator.  This can be done whenever we want to create a visualization
      // of the convex hull, too.

      // Next, we determine the volume of the convex hull.  We do this in a roundabout way:
      //    - Let F be the set of facets on the hull, (F_1, F_2, etc.)
      //    - For each facet, create a tetrahedron with vertex at the origin and the facet as its base.
      //    - Sum the volumes of all the tetrahedra.
      
      // Get unit normal vector to each facet on the convex hull

      for_each(convexhull_.facets_begin(),convexhull_.facets_end(), Facet_normal());

      // Now the centroid vector contains the Cartesian coordinates of the centroid

      Facet_iterator fi,save_fi;
      Halfedge_around_facet_circulator hc;

      CGAL::set_ascii_mode(cout);

      int facecount = 1;
      int edgecount = 1;
      std::vector<Point_3> tpoints;
      tpoints.resize(3);

      Real scalefactor = 0.0;
      Real minProjlength = 1.0e10;
      Real maxTripodarea = 0.0;
      Real tripodarea = 0.0;
      Real projlength = 0.0;
      std::vector<Real> proj,norm,norm1,norm2,norm3,diff,diff1,diff2;
      proj.clear();
      proj.resize(3,0.0);
      norm.clear();
      norm.resize(3,0.0);
      norm1.clear();
      norm1.resize(3,0.0);
      norm2.clear();
      norm2.resize(3,0.0);
      norm3.clear();
      norm3.resize(3,0.0);
      diff1.clear();
      diff1.resize(3,0.0);
      diff2.clear();
      diff2.resize(3,0.0);

      save_fi = convexhull_.facets_begin();
      count = 0;
      for (fi = convexhull_.facets_begin(); fi != convexhull_.facets_end(); ++fi) {
        hc = fi->facet_begin();
        if (CGAL::circulator_size(hc) == 3) {
          edgecount = 0;
          do {
            tpoints[edgecount] = hc->vertex()->point();
            edgecount++;
          } while ( ++hc != fi->facet_begin() );

          // cout << "Tet " << facecount << ": (" << tpoints[0].x() << "," << tpoints[0].y() << "," << tpoints[0].z() << ") ("
          //                                      << tpoints[1].x() << "," << tpoints[1].y() << "," << tpoints[1].z() << ") ("
          //                                      << tpoints[2].x() << "," << tpoints[2].y() << "," << tpoints[2].z() << ") ("
          //                                      << tpoints[3].x() << "," << tpoints[3].y() << "," << tpoints[3].z() << ")" << endl;

          // Calculate projection of centroid onto the triangular facet and decide if it lies inside the triangle
          // Unit normal vector of each facet has already been calculated

          // Consider the 3d pyramid with faces ABC, ABP, BCP, CAP. The projection of P onto ABC is inside it iff
          // the dihedral angles between ABC and each of the other 3 triangles are all less than 90 degrees.
          // In turn, these angles are equal to the angle between N and the respective outward-facing triangle
          // normal! So our algorithm is this:

          norm[0] = fi->normal().x();
          norm[1] = fi->normal().y();
          norm[2] = fi->normal().z();

          diff1[0] = tpoints[1].x() - tpoints[0].x();
          diff1[1] = tpoints[1].y() - tpoints[0].y();
          diff1[2] = tpoints[1].z() - tpoints[0].z();

          diff2[0] = centroid_[0] - tpoints[0].x();
          diff2[1] = centroid_[1] - tpoints[0].y();
          diff2[2] = centroid_[2] - tpoints[0].z();
          
          norm1 = getCartesianProduct(diff1,diff2);

          diff1[0] = tpoints[2].x() - tpoints[1].x();
          diff1[1] = tpoints[2].y() - tpoints[1].y();
          diff1[2] = tpoints[2].z() - tpoints[1].z();

          diff2[0] = centroid_[0] - tpoints[1].x();
          diff2[1] = centroid_[1] - tpoints[1].y();
          diff2[2] = centroid_[2] - tpoints[1].z();
          
          norm2 = getCartesianProduct(diff1,diff2);

          diff1[0] = tpoints[0].x() - tpoints[2].x();
          diff1[1] = tpoints[0].y() - tpoints[2].y();
          diff1[2] = tpoints[0].z() - tpoints[2].z();

          diff2[0] = centroid_[0] - tpoints[2].x();
          diff2[1] = centroid_[1] - tpoints[2].y();
          diff2[2] = centroid_[2] - tpoints[2].z();
          
          norm3 = getCartesianProduct(diff1,diff2);

          if (getInnerProduct(norm1,norm) >= 0.0 && getInnerProduct(norm2,norm) >= 0.0
                && getInnerProduct(norm3,norm) >= 0.0) {

              /*  This method used distance from centroid to tripod as criterion
              // The projection of centroid to this plane passes through the facet in question
              // Now find length and compare to minimum length

              // proj is the vector pointing from the centroid to its projection on the plane

              scalefactor = getInnerProduct(diff2,norm);
              proj[0] = centroid_[0] - (scalefactor * norm[0]);
              proj[1] = centroid_[1] - (scalefactor * norm[1]);
              proj[2] = centroid_[2] - (scalefactor * norm[2]);

              // Find length of proj and see if it is the minimum found so far.

              projlength = getMagnitude(proj);

              if (projlength < minProjlength) {
                  minProjlength = projlength;
                  save_fi = fi;
                  x1 = tpoints[0].x();
                  x2 = tpoints[1].x();
                  x3 = tpoints[2].x();
                  y1 = tpoints[0].y();
                  y2 = tpoints[1].y();
                  y3 = tpoints[2].y();
                  z1 = tpoints[0].z();
                  z2 = tpoints[1].z();
                  z3 = tpoints[2].z();
              }
              */

              /* This method uses the area of the tripod as the criterion */
              diff1[0] = tpoints[1].x() - tpoints[0].x();
              diff1[1] = tpoints[1].y() - tpoints[0].y();
              diff1[2] = tpoints[1].z() - tpoints[0].z();

              diff2[0] = tpoints[2].x() - tpoints[0].x();
              diff2[1] = tpoints[2].y() - tpoints[0].y();
              diff2[2] = tpoints[2].z() - tpoints[0].z();

              tripodarea = 0.5 * getMagnitude(getCartesianProduct(diff1,diff2));
              count++;

              if (tripodarea > maxTripodarea) {
                  maxTripodarea = tripodarea;
                  save_fi = fi;
                  x1 = tpoints[0].x();
                  x2 = tpoints[1].x();
                  x3 = tpoints[2].x();
                  y1 = tpoints[0].y();
                  y2 = tpoints[1].y();
                  y3 = tpoints[2].y();
                  z1 = tpoints[0].z();
                  z2 = tpoints[1].z();
                  z3 = tpoints[2].z();
              }
              /*
              if (count == 30) {
                  maxTripodarea = tripodarea;
                  save_fi = fi;
                  x1 = tpoints[0].x();
                  x2 = tpoints[1].x();
                  x3 = tpoints[2].x();
                  y1 = tpoints[0].y();
                  y2 = tpoints[1].y();
                  y3 = tpoints[2].y();
                  z1 = tpoints[0].z();
                  z2 = tpoints[1].z();
                  z3 = tpoints[2].z();
              }
              */

          }

          // proj is the vector pointing from the centroid to its projection on the facet plane

        } else {
          cout << "  Dangit!  Facet has too many edges!" << endl;
        }
        facecount++;
      }

      // At this point, we have the most stable facet saved, if in fact our definition of facet stability
      // is the correct one to use

      // Maybe an alternative definition of stability is the facet with the largest enclosed area?

      // The next thing to do is rotate the particle so that the most stable facet is parallel to the 
      // xy plane and has the centroid above it.

      // Point a facet iterator at the most stable facet on the hull

      norm[0] = save_fi->normal().x();
      norm[1] = save_fi->normal().y();
      norm[2] = save_fi->normal().z();

      cout << "Particle centroid is now (" << centroid_[0] << "," << centroid_[1] << "," <<  centroid_[2] << ")" << endl;
      cout << "Facet points are:" << endl;
      cout << "    Point 1 = (" << x1 << "," << y1 << "," << z1 << ")" << endl;
      cout << "    Point 2 = (" << x2 << "," << y2 << "," << z2 << ")" << endl;
      cout << "    Point 3 = (" << x3 << "," << y3 << "," << z3 << ")" << endl;
      cout << "Tripod area = " << maxTripodarea << endl;
      cout.flush();
      cout << "Resting facet normal is currently (" << norm[0] << "," << norm[1] << "," << norm[2] << ")" << endl;
      cout.flush();
      // Find orientation of the unit normal in spherical polar coordinates

      Real theta = acos(norm[2]);
      Real phi = atan2(norm[1],norm[0]);

      // Determine the rotations needed to get the unit normal pointing in the -z direction

      Real alpha = -phi;
      Real beta = Pi - theta;

      cout << "Need to rotate by alpha = " << alpha * Rad2Deg << ", beta = " << beta * Rad2Deg << endl;
      cout.flush();

      doRotate(-alpha,-beta,0.0);

      // Recalculate and retrieve the convex hull object, which will be stored in the convexhull_ member

      cout << "Preparing vertices for CGAL convex hull computation..." << endl;
      for (i = 0; i < r_.size(); i++) {
        coords = getX(i);
        if (i == 2) {
            cout << "Point 2 coords = (" << coords[XDIM] << "," << coords[YDIM] << "," << coords[ZDIM] << ")" << endl;
            cout.flush();
        }
        p = Point_3(coords[XDIM],coords[YDIM],coords[ZDIM]);
        points.push_back(p);
      }

      if (verbose_) {
        cout << "Done preparing vector of points for CGAL convex hull computation... " << endl;
      }
        
      CGAL::convex_hull_3(points.begin(),points.end(),convexhull_);
        
      // How can we check that it worked, other than qualitative visualization?

      // Recompute the particle centroid

      computeCentroid();

    }

    catch (tsquare::FileException fex) { throw fex; }
    catch (tsquare::FloatException flex) { throw flex; }
    catch (tsquare::EOBException ex) { throw ex; }
    catch (tsquare::DataException dex) { throw dex; }
}

std::vector<Real> tsquare::Particle::getX(const unsigned int surfindex)
{
    Real t,p,st,r;
    std::vector<Real> xp;
    xp.clear();
    xp.resize(3,0.0);
    try {
        xp[XDIM] = x_.at(surfindex);
        xp[YDIM] = y_.at(surfindex);
        xp[ZDIM] = z_.at(surfindex);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getX","surface_",
                        surface_.size(),surfindex);
        ex.printException();
    }
    return xp;
}

Real tsquare::Particle::getFnm(const Real n, const Real m)
{
    Real fnm;
    int in = (int)n;
    int im = (int)m;
    Real factor = fac(in-im)/fac(in+im);
    fnm = std::sqrt((((2.0 * n) + 1.0) * factor)/(4.0 * Pi));
    return fnm;
}

Real tsquare::Particle::getE(const unsigned int si)
{
    using std::cout;
    using std::endl;

    Real theta,phi;
    Real E = 0.0;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        cout << "Particle::getE(" << theta << "," << phi << ")" << endl;
        E = getE(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getE","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return E;
}

Real tsquare::Particle::getF(const unsigned int si)
{
    using std::cout;
    using std::endl;

    Real theta,phi;
    Real F = 0.0;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        cout << "Particle::getF(" << theta << "," << phi << ")" << endl;
        F = getF(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getF","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return F;
}

Real tsquare::Particle::getG(const unsigned int si)
{
    using std::cout;
    using std::endl;

    Real theta,phi;
    Real G = 0.0;
    try {
        theta = surface_[si][0];
        phi = surface_[si][1];
        cout << "Particle::getG(" << theta << "," << phi << ")" << endl;
        G = getG(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getG","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return G;
}

std::vector<Real> tsquare::Particle::getNormal(const unsigned int si)
{
    try { return (normal_.at(si)); }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Particle","getNormal","normal_",
                        normal_.size(),si);
    }
}

Real tsquare::Particle::getL(const unsigned int si)
{
    Real theta,phi;
    Real L = 0.0;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        L = getL(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getNormal","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return L;
}

Real tsquare::Particle::getL(const Real t, const Real p)
{
    // t is theta, p is phi

    Real L = (getNormal(XDIM,t,p) * getXtt(XDIM,t,p))
              + (getNormal(YDIM,t,p) * getXtt(YDIM,t,p))
              + (getNormal(ZDIM,t,p) * getXtt(ZDIM,t,p));
    return L;
}

Real tsquare::Particle::getM(const unsigned int si)
{
    Real theta,phi;
    Real M = 0.0;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        M = getM(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getM","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return M;
}

Real tsquare::Particle::getM(const Real t, const Real p)
{
    // t is theta, p is phi

    Real M = (getNormal(XDIM,t,p) * getXtp(XDIM,t,p))
              + (getNormal(YDIM,t,p) * getXtp(YDIM,t,p))
              + (getNormal(ZDIM,t,p) * getXtp(ZDIM,t,p));
    return M;
}

Real tsquare::Particle::getN(const unsigned int si)
{
    Real theta,phi;
    Real N = 0.0;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        N = getN(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getN","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return N;
}

Real tsquare::Particle::getN(const Real t, const Real p)
{
    // t is theta, p is phi

    Real N = (getNormal(XDIM,t,p) * getXpp(XDIM,t,p))
              + (getNormal(YDIM,t,p) * getXpp(YDIM,t,p))
              + (getNormal(ZDIM,t,p) * getXpp(ZDIM,t,p));
    return N;
}

bool tsquare::Particle::isEllipticumbilic(const unsigned int si)
{
    bool answer = false;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        answer = isEllipticumbilic(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","isEllipticumbilic","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return answer;
}

bool tsquare::Particle::isEllipticumbilic(const Real theta, const Real phi)
{
    using std::fabs;

    Real E = getE(theta,phi);
    Real F = getF(theta,phi);
    Real G = getG(theta,phi);
    Real L = getL(theta,phi);
    Real M = getM(theta,phi);
    Real N = getN(theta,phi);
    
    Real lambda;
    
    if (E != 0.0) {
        lambda = L/E;
    } else {
        tsquare::FloatException flex("Particle","isEllipticumbilic","Divide by E = 0");
        flex.printException();
        return false;
    }
    Real b = (L * N) - (M * M);
    Real g = (E * G) - (F * F);
    if (lambda != 0.0 && g != 0.0) {
        if (fabs(b/g - (lambda * lambda)) < 1.0e-20) return true;
    }
    return false;
}

bool tsquare::Particle::isParabolicumbilic(const unsigned int si)
{
    bool answer = false;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        answer = isParabolicumbilic(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","isParabolicumibilic","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return answer;
}

bool tsquare::Particle::isParabolicumbilic(const Real theta, const Real phi)
{
    Real E = getE(theta,phi);
    Real L = getL(theta,phi);
    
    Real lambda = 1.0e6;
    
    if (E != 0.0) {
        lambda = L/E;
    } else {
        tsquare::FloatException flex("Particle","isEllipticumbilic","Divide by E = 0");
        flex.printException();
        return false;
    }
    if (fabs(lambda) < 1.0e-9) return true;
    return false;
}

bool tsquare::Particle::isUmbilic(const unsigned int si)
{
    bool answer = false;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        answer = isUmbilic(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","isUmbilic","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return answer;
}

bool tsquare::Particle::isUmbilic(const Real theta, const Real phi)
{
    return ((isEllipticumbilic(theta,phi) || isParabolicumbilic(theta,phi)));
}

Real tsquare::Particle::getK(const unsigned int si)
{
    using std::cout;
    using std::endl;

    Real K = 0.0;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        K = getK(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getK","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    catch (tsquare::FloatException flex) {
        flex.printException();
        cout << "Returning K = 0 as default" << endl;
        cout.flush();
    }
    return K;
}

Real tsquare::Particle::getK(const Real theta, const Real phi)
{
    // Gaussian curvature
    Real K = 0.0;
    Real E = getE(theta,phi);
    Real F = getF(theta,phi);
    Real G = getG(theta,phi);
    Real L = getL(theta,phi);
    Real M = getM(theta,phi);
    Real N = getN(theta,phi);
    
    if (((E * G) - (F * F)) != 0.0) {
        K = ((L * N) - (M * M)) / ((E * G) - (F * F));
    } else {
        throw tsquare::FloatException("Particle","getK","Divide by (EG - FF) = 0");
    }
    return K;
}

Real tsquare::Particle::getH(const unsigned int si)
{
    using std::cout;
    using std::endl;

    Real H = 0.0;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        H = getH(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getH","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    catch (tsquare::FloatException flex) {
        flex.printException();
        cout << "Returning H = 0 as default" << endl;
        cout.flush();
    }
    return H;
}

Real tsquare::Particle::getH(const Real theta, const Real phi)
{
    // Mean curvature
    Real H = 0.0;
    Real E = getE(theta,phi);
    Real F = getF(theta,phi);
    Real G = getG(theta,phi);
    Real L = getL(theta,phi);
    Real M = getM(theta,phi);
    Real N = getN(theta,phi);
    
    if (((E * G) - (F * F)) != 0.0) {
        H = ((L * G) - (2.0 * M * F) + (N * E))
              / ((E * G) - (F * F));
    } else {
        throw tsquare::FloatException("Particle","getH","Divide by (EG - FF) = 0");
    }
    return (-0.5 * H);
}

std::vector<Real> tsquare::Particle::getPrincipaldirection(const unsigned int index, const unsigned int si)
{
    std::vector<Real> pd;
    pd.clear();
    pd.resize(2,0.0);
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        pd = getPrincipaldirection(index,theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getPrincipaldirection","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return pd;
}

std::vector<Real> tsquare::Particle::getPrincipaldirection(const unsigned int index,
                                             const Real theta, const Real phi)
{
    using std::sqrt;

    std::vector<Real> pd;
    pd.clear();
    pd.resize(2,0.0);
    Real mag = 1.0;
    Real E = getE(theta,phi);
    Real F = getF(theta,phi);
    Real G = getG(theta,phi);
    Real L = getL(theta,phi);
    Real M = getM(theta,phi);
    Real N = getN(theta,phi);
    
    // Compute principal curvature manually since fundamental forms already computed
    Real k = 1.0;
    Real a = (E * G) - (F * F);
    Real b = (2.0 * F * M) - (E * N) - (G * L);
    Real c = (L * N) - (M * M);
    Real sgn = 1.0;
    if (a != 0.0) {
        k = (-b + (sgn * sqrt((b*b) - (4.0 * a * c))))/(2.0 * a);
    } else if (b != 0.0) {
        k = (-c / b);
    } else {
        tsquare::FloatException flex("Particle","getPrincipaldirection","Divide by zero error");
        flex.printException();
        k = (c > 0.0) ? -1.0e10 : 1.0e10;
    }
    // Done computing principal curvature, k

    Real denom = (N * F) - (M * G);
    if (denom != 0.0) {
        pd[0] = 1.0;
        pd[1] = ((L * G) + k - (M * F)) / denom;
        mag = sqrt(1.0 + (pd[1] * pd[1]));
        pd[0] /= mag;
        pd[1] /= mag;
    } else {
        tsquare::FloatException flex("Particle","getPrincipaldirection","Divide by zero error");
        flex.printException();
        pd[0] = 0.0;
        pd[1] = 1.0;
    }

    if (index == 0) {
        return pd;
    } else {
        Real pd0 = pd[0];
        Real pd1 = pd[1];
        pd[0] = 0.0;
        pd[1] = 0.0;
        if (pd0 != 0.0) {
            pd[1] = 1.0;
            pd[0] = -pd1/pd0;
        } else {
            pd[0] = 1.0;
            pd[1] = -pd0;
        }
        mag = sqrt((pd[0] * pd[0]) + (pd[1] * pd[1]));
        if (mag != 0.0) {
            pd[0] /= mag;
            pd[1] /= mag;
        } else {
            tsquare::FloatException flex("Particle","getPrincipaldirection","Divide by zero error");
            flex.printException();
            pd[0] = 0.0;
            pd[1] = 0.0;
        }
    }
    return pd;
}

Real tsquare::Particle::getPrincipalcurvature(const unsigned int index, const unsigned int si)
{
    Real theta,phi,kappa;
    if (si < pcurv_.size()) {
        try { return (pcurv_.at(si).at(index)); }
        catch (std::out_of_range &oor) {
            throw tsquare::EOBException("Particle","getPrincipalcurvature","pcurv_[si]",
                            pcurv_[si].size(),index);
        }
    } else {
        try {
            theta = surface_.at(si)[0];
            phi = surface_.at(si)[1];
            kappa = computePrincipalcurvature(index,theta,phi);
        }
        catch (std::out_of_range &oor) {
            throw tsquare::EOBException("Particle","getPrincipalcurvature","surface_",
                            surface_.size(),si);
        }
        catch (tsquare::FloatException flex) { throw flex; }
        return (kappa);
    }
}

Real tsquare::Particle::computePrincipalcurvature(const unsigned int index,
                                         const Real theta, const Real phi)
{
    using std::sqrt;

    Real kappa = 1.0;
    Real E = getE(theta,phi);
    Real F = getF(theta,phi);
    Real G = getG(theta,phi);
    Real L = getL(theta,phi);
    Real M = getM(theta,phi);
    Real N = getN(theta,phi);
    
    Real a = (E * G) - (F * F);
    Real b = (2.0 * F * M) - (E * N) - (G * L);
    Real c = (L * N) - (M * M);
    Real sgn = 1.0;
    if (index != 0) sgn = -1.0;
    if (a != 0.0) {
        kappa = (-b + (sgn * sqrt((b*b) - (4.0 * a * c))))/(2.0 * a);
    } else if (b != 0.0) {
        kappa = (-c / b);
    } else {
        throw tsquare::FloatException("Particle","getPrincipaldirection","Divide by zero error");
    }

    return (-kappa);
}


Real tsquare::Particle::getEllipticumbiliccurvature(const unsigned int index, const unsigned int si)
{
    // This function returns the normal curvature at an elliptic umbilic point

    Real kappa = 0.0;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        kappa = getEllipticumbiliccurvature(index,theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getEllipicumbiliccurvature","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return kappa;
}

Real tsquare::Particle::getEllipticumbiliccurvature(const unsigned int index,
                                           const Real theta, const Real phi)
{
    // This function returns the normal curvature at an elliptic umbilic point

    if (!(isEllipticumbilic(theta,phi))) {
        tsquare::DataException dex("Particle","getEllipticumbilicurvature",
                          "Point is not an elliptic umbilic");
        dex.printException();
        return(0.0);
    }

    // We let the parameterization be in the direction (0,1), chosen arbitrarily because
    // all directions should give the same normal curvature at an elliptic umbilic
    // Therefore u' = 0 and v' = 1.

    Real E = getE(theta,phi);
    Real kappa = 0.0;
    if (E != 0.0) {
        kappa = getL(theta,phi)/E;
    } else {
        tsquare::FloatException flex("Particle","getEllipticumbiliccurvature",
                            "Divide by E = 0");
        flex.printException();
        return (0.0);
    }

    return (-kappa);
}

std::vector<Real> tsquare::Particle::getPrincipalcurvatures(const unsigned int si)
{
    Real r = getR(si);
    std::vector<Real> kappa;
    kappa.clear();
    kappa.resize(2,0.0);
    try { kappa = pcurv_[si]; }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getPrincipalcurvatures","pcurv_",
                        pcurv_.size(),si);
        ex.printException();
    }
    return kappa;
}

std::vector<Real> tsquare::Particle::computePrincipalcurvatures(const Real theta, const Real phi)
{
    using std::sqrt;

    Real E = getE(theta,phi);
    Real F = getF(theta,phi);
    Real G = getG(theta,phi);
    Real L = getL(theta,phi);
    Real M = getM(theta,phi);
    Real N = getN(theta,phi);
    
    Real a = (E * G) - (F * F);
    Real b = (2.0 * F * M) - (E * N) - (G * L);
    Real c = (L * N) - (M * M);
    Real mag;
    // Real sgn = 1.0;
    // if (index != 0) sgn = -1.0;

    std::vector<Real> kappa;
    kappa.clear();
    if (a != 0.0) {
        kappa.push_back(-(-b + sqrt((b*b) - (4.0 * a * c)))/(2.0 * a));
        kappa.push_back(-(-b - sqrt((b*b) - (4.0 * a * c)))/(2.0 * a));
    } else if (b != 0.0) {
        kappa.push_back(c / b);
    } else {
        tsquare::FloatException flex("Particle","getPrincipalcurvatures","Divide by zero error");
        flex.printException();
        mag = (c > 0.0) ? -1.0e10 : 1.0e10;
        kappa.push_back(-mag);
    }

    return kappa;
}

void tsquare::Particle::computePrincipalcurvatures(void)
{
    using std::cout;
    using std::endl;
    using std::sqrt;

    Real theta,phi,r,rp,rt,rpp,rtt,rtp,st,ct,sp,cp;
    Real E,F,G,L,M,N,a,b,c,mag;
    int count;
    std::vector<Real> kappa;
    kappa.resize(2,0.0);

    if (pcurv_.size() == 0) {
        kappa.clear();
        kappa.resize(2,0.0);
        pcurv_.clear();
        count = 0;
        for (register int i = 1; i <= ntheta_; i++) {
          theta = surface_[count][0];
          if (verbose_) {
              cout << "    Curvature computation:  Latitude " << i << " out of " << ntheta_ << endl;
              cout.flush();
          }
          for (register int j = 1; j <= nphi_; j++) {
            phi = surface_[count][1];
            E = getE(theta,phi);
            F = getF(theta,phi);
            G = getG(theta,phi);
            L = getL(theta,phi);
            M = getM(theta,phi);
            N = getN(theta,phi);
    
            a = (E * G) - (F * F);
            b = (2.0 * F * M) - (E * N) - (G * L);
            c = (L * N) - (M * M);
            // sgn = 1.0;
            // if (index != 0) sgn = -1.0;

            if (a != 0.0) {
                kappa[0] = (-(-b + sqrt((b*b) - (4.0 * a * c)))/(2.0 * a));
                kappa[1] = (-(-b - sqrt((b*b) - (4.0 * a * c)))/(2.0 * a));
            } else if (b != 0.0) {
                kappa[0] = kappa[1] = (c / b);
            } else {
                tsquare::FloatException flex("Particle","getPrincipalcurvatures","Divide by zero error");
                flex.printException();
                mag = (c > 0.0) ? -1.0e10 : 1.0e10;
                kappa[0] = kappa[1] = (-mag);
            }

            pcurv_.push_back(kappa);
            count++;
          }
        }
    }
    return;
}

Real tsquare::Particle::getMincurvature(const unsigned int si)
{
    Real kappa = 0.0;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        kappa = getMincurvature(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getMincurvature","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return kappa;
}

Real tsquare::Particle::getMincurvature(const Real theta, const Real phi)
{
    std::vector<Real> k = computePrincipalcurvatures(theta,phi);
    if (k.size() == 1) {
        return k[0];
    } else if (k[0] <= k[1]) {
        return k[0];
    } else {
        return k[1];
    }
}

Real tsquare::Particle::getMaxcurvature(const unsigned int si)
{
    Real kappa = 0.0;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        kappa = getMaxcurvature(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getMaxcurvature","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return kappa;
}

Real tsquare::Particle::getMaxcurvature(const Real theta, const Real phi)
{
    std::vector<Real> k = computePrincipalcurvatures(theta,phi);
    if (k.size() == 1) {
        return k[0];
    } else if (k[0] <= k[1]) {
        return k[1];
    } else {
        return k[0];
    }
}

Real tsquare::Particle::getMinabscurvature(const unsigned int si)
{
    Real kappa = 0.0;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        kappa = getMinabscurvature(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getMinabscurvature","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return kappa;
}

Real tsquare::Particle::getMinabscurvature(const Real theta, const Real phi)
{
    std::vector<Real> k = computePrincipalcurvatures(theta,phi);
    if (k.size() == 1) {
        return fabs(k[0]);
    } else if (fabs(k[0]) <= fabs(k[1])) {
        return fabs(k[0]);
    } else {
        return fabs(k[1]);
    }
}

Real tsquare::Particle::getMaxabscurvature(const unsigned int si)
{
    Real kappa = 0.0;
    Real theta,phi;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        kappa = getMaxabscurvature(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getMaxabscurvature","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return kappa;
}

Real tsquare::Particle::getMaxabscurvature(const Real theta, const Real phi)
{
    std::vector<Real> k = computePrincipalcurvatures(theta,phi);
    if (k.size() == 1) {
        return fabs(k[0]);
    } else if (fabs(k[0]) <= fabs(k[1])) {
        return fabs(k[1]);
    } else {
        return fabs(k[0]);
    }
}

Real tsquare::Particle::getASF(void)
{
    using std::pow;
    using std::sqrt;

    Real d1,d2,d3;
    Real term1,term2,numerator;
    d1 = triaxialdim_[0];
    d2 = triaxialdim_[1];
    d3 = triaxialdim_[2];

    if (d1 > 0.0 && d2 > 0.0 && d3 > 0.0) {
        term1 = 1.0 + (d3 / d2) * (1.0 + (d2 / d1));
        term2 = sqrt(1.0 + ((pow((d3 / d2),2.0)) * (1.0 + pow((d2 / d1),2.0))));
        numerator = pow((pow((d3 / d2),2.0)*(d2 / d1)),(1.0/3.0));
        return (13.392 * numerator / (term1 + (6.0 * term2)));
    } else {
        throw tsquare::FloatException("Particle","getASF","ASF denominator is zero");
   }
}

Real tsquare::Particle::getShapeentropy(void)
{
    using std::log;

    Real d1,d2,d3;
    d1 = triaxialdim_[0];
    d2 = triaxialdim_[1];
    d3 = triaxialdim_[2];

    Real denom = d3 + d2 + d1;
    if (denom > 0.0) {
        Real p1 = d1 / denom;
        Real p2 = d2 / denom;
        Real p3 = d3 / denom;
        if (p1 > 0.0 && p2 > 0.0 && p3 > 0.0) {
            return (-((p1 * log(p1))+(p2 * log(p2))+(p3 * log(p3)))/(log(3.0)));
        } else {
            throw tsquare::FloatException("Particle","getShapeentropy",
                                 "Negative argument to log function");
        }
    } else {
        throw tsquare::FloatException("Particle","getShapeentropy","Divide by zero");
    }
}

tsquare::Particle& tsquare::Particle::operator=(const tsquare::Particle& that)
{
    name_ = that.getName();
    xlo_ = that.getXlo();
    ylo_ = that.getYlo();
    zlo_ = that.getZlo();
    xhi_ = that.getXhi();
    yhi_ = that.getYhi();
    zhi_ = that.getZhi();
    volume_ = that.getVolume();
    area_ = that.getArea();
    narea_ = that.getNarea();
    diam_ = that.getDiam();
    dim_ = that.getDim();
    ndim_ = that.getNdim();
    triaxialdim_ = that.getTriaxialdim();
    itrace_ = that.getItrace();
    ngc_ = that.getNgc();
    nmax_ = that.getNmax();
    a_ = that.getA();
    ntheta_ = that.getNtheta();
    nphi_ = that.getNphi();
    surface_ = that.getSurface();
    r_ = that.getR();

    return *this;
}

Real tsquare::Particle::getR(const unsigned int surfaceindex)
{
    try {
        return r_.at(surfaceindex);
    }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Particle","getR","r_",
                            r_.size(),1);
    }
}

Real tsquare::Particle::getRp(const unsigned int surfaceindex)
{
    Real theta,phi,rp;
    try {
        theta = surface_.at(surfaceindex)[0];
        phi = surface_.at(surfaceindex)[1];
        rp = getRp(theta,phi);
        return (rp);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getRp","surface_",
                        surface_.size(),surfaceindex);
        ex.printException();
        exit(1);
    }
}

Real tsquare::Particle::getRt(const unsigned int surfaceindex)
{
    Real theta,phi,rt;
    try {
        theta = surface_.at(surfaceindex)[0];
        phi = surface_.at(surfaceindex)[1];
        rt = getRt(theta,phi);
        return (rt);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getRt","surface_",
                        surface_.size(),surfaceindex);
        ex.printException();
        exit(1);
    }
}

Real tsquare::Particle::getRpp(const unsigned int surfaceindex)
{
    Real theta,phi,rpp;
    try {
        theta = surface_.at(surfaceindex)[0];
        phi = surface_.at(surfaceindex)[1];
        rpp = getRpp(theta,phi);
        return (rpp);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getRpp","surface_",
                        surface_.size(),surfaceindex);
        ex.printException();
        exit(1);
    }
}

Real tsquare::Particle::getRtt(const unsigned int surfaceindex)
{
    Real theta,phi,rtt;
    try {
        theta = surface_.at(surfaceindex)[0];
        phi = surface_.at(surfaceindex)[1];
        rtt = getRtt(theta,phi);
        return (rtt);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getRtt","surface_",
                        surface_.size(),surfaceindex);
        ex.printException();
        exit(1);
    }
}

Real tsquare::Particle::getRtp(const unsigned int surfaceindex)
{
    Real theta,phi,rtp;
    try {
        theta = surface_.at(surfaceindex)[0];
        phi = surface_.at(surfaceindex)[1];
        rtp = getRtp(theta,phi);
        return (rtp);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Particle","getRtp","surface_",
                        surface_.size(),surfaceindex);
        ex.printException();
        exit(1);
    }
}
