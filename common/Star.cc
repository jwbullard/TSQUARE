/**
@file
@brief Method definitions for Star class.
*/
#include "Star.h"

tsquare::Star::Star()
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

tsquare::Star::Star(std::string &fname, bool verbose) {
    std::string buff;
    int nval,mval;
    Real rval,ival;
    std::complex<Real> cnum;
    verbose_ = verbose;

    try {
        std::ifstream in(fname.c_str());
        if (!in) {
            throw tsquare::FileException("Star","Star",fname,
                                "Could not open for reading");
        }

        while (!in.eof()) {
            in >> nval >> mval >> rval >> ival;
        }
        nmax_from_anmfile_ = nval;
        nmax_ = nmax_from_anmfile_;

        in.close();
        in.open(fname.c_str());
        if (!in) {
            throw tsquare::FileException("Star","Star",fname,
                                "Could not open for reading");
        }
        
        std::vector<tsquare::ExtendedVector> ev;
        ev.clear();
        ev.resize(nmax_+1);
        for (int i = 0; i <= nmax_; i++) {
            ev[i].clear();
            ev[i].resize(-i,i);
        }

        a_.clear();
        a_.resize(1,ev);

        for (int n = 0; n <= nmax_; n++) {
            for (int m = -n; m <= n; m++) {
                a_[RDIM][n][m] = std::complex<Real>(0.0,0.0);
            }
        }

        while (!in.eof() && nval <= nmax_) {
            in >> nval >> mval >> rval >> ival;
            if (nval <= nmax_) {
                a_[RDIM][nval][mval] = std::complex<Real>(rval,ival);
            } else {
                throw tsquare::DataException("Star","Star","SH order greater than degree");
            }
        }

        in.close();
    }

    catch (tsquare::FileException fex) { throw fex; }
    catch (tsquare::DataException dex) { throw dex; }

    name_ = fname;
    surface_.clear();
    r_.clear();
    x_.clear();
    y_.clear();
    z_.clear();
    xg_.clear();
    wg_.clear();
    dim_.clear();
    ndim_.clear();
    triaxialdim_.clear();
    dim_.resize(3,0.0);
    ndim_ = dim_;
    triaxialdim_ = dim_;

    inscribedsphere_.setXr(0.0);
    inscribedsphere_.setYr(0.0);
    inscribedsphere_.setZr(0.0);
    inscribedsphere_.setRadius(0.0);

    enclosingsphere_.setXr(0.0);
    enclosingsphere_.setYr(0.0);
    enclosingsphere_.setZr(0.0);
    enclosingsphere_.setRadius(0.0);
}

void tsquare::Star::readSHcoeffs(std::string &fname)
{
    using std::cout;
    using std::endl;

    std::string buff;
    int nval,mval;
    Real rval,ival;
    std::complex<Real> cnum;

    try {
        std::ifstream in(fname.c_str());
        if (!in) {
            throw tsquare::FileException("Star","readSHcoeffs",fname,
                                "Could not open for reading");
        }

        while (!in.eof()) {
            in >> nval >> mval >> rval >> ival;
            nmax_ = nval;
        }

        in.close();
        in.open(fname.c_str());
        if (!in) {
            throw tsquare::FileException("Star","readSHcoeffs",fname,
                                "Could not open for reading");
        }
        
        std::vector<tsquare::ExtendedVector> ev;
        ev.clear();
        ev.resize(nmax_+1);
        for (int i = 0; i <= nmax_; i++) {
            ev[i].clear();
            ev[i].resize(-i,i);
        }

        a_.clear();
        a_.resize(1,ev);

        for (int n = 0; n <= nmax_; n++) {
            for (int m = -n; m <= n; m++) {
                a_[RDIM][n][m] = std::complex<Real>(0.0,0.0);
            }
        }

        while (!in.eof() && nval <= nmax_) {
            in >> nval >> mval >> rval >> ival;

            if (verbose_) cout << "Done: " << nval << ", " << mval << ", "
                               << rval << ", " << ival << endl;
            if (nval <= nmax_) {
                a_[RDIM][nval][mval] = std::complex<Real>(rval,ival);
            } else {
                throw tsquare::DataException("Star","readSHcoeffs","SH order greater than degree");
            }
        }

        in.close();
    }

    catch (tsquare::FileException fex) {
        fex.printException();
        exit(1);
    }
    catch (tsquare::DataException dex) {
        dex.printException();
        exit(1);
    }

    surface_.clear();
    r_.clear();
    x_.clear();
    y_.clear();
    z_.clear();
    xg_.clear();
    wg_.clear();

    inscribedsphere_.setXr(0.0);
    inscribedsphere_.setYr(0.0);
    inscribedsphere_.setZr(0.0);
    inscribedsphere_.setRadius(0.0);

    enclosingsphere_.setXr(0.0);
    enclosingsphere_.setYr(0.0);
    enclosingsphere_.setZr(0.0);
    enclosingsphere_.setRadius(0.0);
}

void tsquare::Star::calculateSHcoeffs(std::string &fname, std::string &outfilename,
                                 int numGausspts, int maxDegree, bool surfaceMesh)
{
    using std::cout;
    using std::endl;
    using std::sin;
    using std::cos;

    std::string buff;
    int nVal,mVal,inputVal,xBB,yBB,zBB,count,id;
    std::vector<Real> centroid;
    Real rval,ival,xsum,ysum,zsum,theta,phi,xval,yval,zval,vol,factor,vscale;
    std::complex<Real> cnum,spharm,cterm,caddend;
    std::vector<std::vector<std::vector<Real> > > voxel;

    std::complex<Real> prefact = std::complex<Real>(0.5 * Pi * Pi,0.0);

    nmax_ = nmax_from_anmfile_ = maxDegree;

    centroid.clear();
    centroid.resize(3,0.0);
    voxel.clear();

    // Clear out the array of SH coefficients and reset all to zero
    
    std::vector<tsquare::ExtendedVector> ev;
    ev.clear();
    ev.resize(nmax_+1);
    for (int i = 0; i <= nmax_; i++) {
        ev[i].clear();
        ev[i].resize(-i,i);
    }

    a_.clear();
    a_.resize(1,ev);

    for (int n = 0; n <= nmax_; n++) {
        for (int m = -n; m <= n; m++) {
            a_[RDIM][n][m] = std::complex<Real>(0.0,0.0);
        }
    }

    // Next get the Gaussian quadrature points.  We will need them regardless
    // of whether we are using a surface mesh or raw voxel data

    cout << "Getting " << numGausspts << " Gaussian quadrature points";
    cout.flush();

    surface_.clear();
    getGausspoints(numGausspts);

    // Open the input voxel files and read the solid voxel coordinates

    voxel = getVoxelImage(fname,xBB,yBB,zBB,vscale);

    // Calculate the centroid

    centroid = getImageCentroid(voxel);

    cout << "Found centroid at (" << centroid[XDIM] << ","
         << centroid[YDIM] << "," << centroid[ZDIM] << ")" << endl;
    cout.flush();

    if (surfaceMesh) {

        // The next chunk creates the 3D polyhedron, storing it in surface_poly

        Polyhedron surface_poly = getSurfacePolyhedronFromImage(fname,centroid,xBB,yBB,zBB);
        
        // First translate its centroid to the origin
        const CartesianVector translation_vector(-centroid[0],-centroid[1],-centroid[2]);

        Aff_transformation_3 transl(CGAL::TRANSLATION, translation_vector);
        transform(surface_poly.points_begin(),surface_poly.points_end(),
                  surface_poly.points_begin(),transl);

        // Now the centroid is the origin

        centroid.resize(3,0.0);
        CartesianPoint origin(0.0,0.0,0.0);
     
        // Construct the AABB tree for quick intersection queries

        cout << "Creating AABB tree from polyhedron" << endl;
        cout.flush();

        // The line below is when using AABB_polyhedron_triangle_primitive.h
        // Tree tree(surface_poly.facets_begin(),surface_poly.facets_end());

        // The line below is when using AABB_face_graph_triangle_primitive.h
        Tree tree(faces(surface_poly).first,faces(surface_poly).second,surface_poly);
    
        // Object intersection will hold the point of intersection with the surface

        boost::optional<Object_and_primitive_id> intersection;

        Real bounding_sphere_squared_radius = 2.0 * ((xBB*xBB) + (yBB*yBB) + (zBB*zBB));

        try {
            count = 0;
            vol = 0.0;  // object volume
            for (unsigned int nt = 0; nt < ntheta_; nt++) {
                theta = 0.5 * Pi * (xg_[nt] + 1.0);
                for (unsigned int np = 0; np < nphi_; np++) {
                    phi = Pi * (xg_[np] + 1.0);
                    // theta = surface_[count][0];
                    // phi = surface_[count][1];

                    // Construct a sample ray from the centroid
                    // in the (theta,phi) direction

                    xval = (bounding_sphere_squared_radius * sin(theta)*cos(phi));
                    yval = (bounding_sphere_squared_radius * sin(theta)*sin(phi));
                    zval = (bounding_sphere_squared_radius * cos(theta));

                    CartesianPoint b(xval,yval,zval);

                    CartesianSegment ray_query(origin,b);

                    // Compute first encountered ray intersection(s) with surface_poly

                    intersection = tree.any_intersection(ray_query);
                    if (intersection) {

                        // Get intersection object
                        Object_and_primitive_id op = *intersection;
                        CGAL::Object object = op.first;
                        CartesianPoint point;
                        if (CGAL::assign(point,object)) {

                            // Compute distance to point from centroid

                            xval = point.x();
                            yval = point.y();
                            zval = point.z();
                            rval = sqrt( (xval * xval) + (yval * yval) + (zval * zval) );
                    
                            factor = sin(theta) * wg_[nt] * wg_[np];
                            cterm = std::complex<Real>(rval*factor,0.0);
                            spharm = boost::math::spherical_harmonic(0,0,theta,phi);
                            caddend = cterm * std::conj(spharm);
                            a_[0][0][0] += (caddend);
                            vol += (factor * rval * rval * rval / 3.0);

                            for (register int n = 1; n <= nmax_; n++) {
                                for (register int m = n; m >= -n; m--) {
                                    spharm = boost::math::spherical_harmonic(n,m,theta,phi);
                                    caddend = cterm * std::conj(spharm);
                                    a_[0][n][m] += (caddend);
                                }
                            }
                        
                        } else {
                            cout << "Dangit! Intersection with surface is NOT a point" << endl;
                            throw tsquare::DataException("Star","calculateSHcoeffs",
                                  "Intersection of ray with meshed surface is not a point");
                        }

                    } else {
                        cout << "Dangit! No intersection with surface was found "
                             << "for orientation (" << theta << "," << phi << ")" << endl;
                        throw tsquare::DataException("Star","calculateSHcoeffs",
                                          "No intersection of ray with meshed surface");
                    }

                    count++;
  
                }  // End of phi loop
            }      // End of theta loop
            
        }

        catch (tsquare::DataException dex) {
            dex.printException();
            return;
        }

    } else {    // We are dealing only with a voxel collection
                // Not a surface mesh

        // For each direction, we crawl along the ray and check when we get to
        // a surface voxel, then interpolate to get the radius

        Real xval2,yval2,zval2,rstart;
        Real rr1,rr2,st,sp,ct,cp;
        int i1,i2,j1,j2,k1,k2,id2;

        Real len = (Real)xBB;       // len is minimum dimension of bounding box
        if (yBB < len) len = yBB;
        if (zBB < len) len = zBB;
        
        Real cutoff = len / 100000.0;
        Real dd = len / 10.0;
        bool done = false;
      
        count = 0;
        vol = 0.0;  // object volume

        cout << "Calculating on binary image directly" << endl;
        cout.flush();
        for (unsigned int nt = 0; nt < ntheta_; nt++) {
            for (unsigned int np = 0; np < nphi_; np++) {
                theta = surface_[count][0];
                phi = surface_[count][1];
                st = sin(theta);
                sp = sin(phi);
                ct = cos(theta);
                cp = cos(phi);

                // First we must find the distance to the surface in this direction

                rstart = 0.0;
                done = false;
                do {
                    for (register int ii = 1; ii <= 10; ii++) {
                        rr1 = rstart + ((ii - 1) * dd);
                        rr2 = rr1 + dd;
                        done = false;

                        // position with respect to centroid

                        xval = rr1 * st * cp;
                        yval = rr1 * st * sp;
                        zval = rr1 * ct;
                        xval2 = rr2 * st * cp;
                        yval2 = rr2 * st * sp;
                        zval2 = rr2 * ct;

                        // actual position

                        xval += centroid[0];
                        yval += centroid[1];
                        zval += centroid[2];
                        xval2 += centroid[0];
                        yval2 += centroid[1];
                        zval2 += centroid[2];

                        // nearest voxel center

                        i1 = (int)(xval + 0.5);
                        j1 = (int)(yval + 0.5);
                        k1 = (int)(zval + 0.5);
                        i2 = (int)(xval2 + 0.5);
                        j2 = (int)(yval2 + 0.5);
                        k2 = (int)(zval2 + 0.5);

                        id = i1 + (j1 * xBB) + (k1 * xBB * yBB);
                        id2 = i2 + (j2 * xBB) + (k2 * xBB * yBB);

                        // check if 2nd point goes beyond image box on the first try
                        if (i2 > xBB || i2 < 0 || j2 > yBB || j2 < 0
                                     || k2 > zBB || k2 < 0) {
                            dd *= 0.5;

                        } else if ((int)(voxel[i1][j1][k1]) != 0 && (int)(voxel[i2][j2][k2]) != 0) {

                            if (ii == 10) rstart = rr2;

                        } else if ((int)(voxel[i1][j1][k1]) != 0 && (int)(voxel[i2][j2][k2]) == 0) {

                            if (fabs(dd) < cutoff) {
                                done = true;
                            } else {
                                rstart = rr1;
                                dd *= 0.1;
                            }

                        }
                    }

                } while (!done);

                rval = 0.5 * (rr1 + rr2);
                factor = sin(theta) * wg_[nt] * wg_[np];
                cout << "Radius(" << theta * 180.0 / Pi << "," << phi * 180.0 / Pi << ") = "
                     << rval << "; factor = " << factor << endl;
                cout.flush();
                spharm = boost::math::spherical_harmonic(0,0,theta,phi);
                a_[RDIM][0][0] += (rval * factor * std::conj(spharm));
                vol += (factor * rval * rval * rval / 3.0);

                for (register int n = 1; n <= nmax_; n++) {
                    for (register int m = n; m >= -n; m--) {
                        spharm = boost::math::spherical_harmonic(n,m,theta,phi);
                        a_[RDIM][n][m] += (rval * factor * std::conj(spharm));
                    }
                }

                count++;

            }  // End of phi loop
        }      // End of theta loop

    }      // Done with block dealing with voxel image input

    // By now, all SH coefficients have been computed one way or the other
    // We need to multiply each one by (1/2) Pi^2

    vol *= real(prefact);
    a_[RDIM][0][0] *= prefact;
    for (register int n = 1; n <= nmax_; n++) {
        for (register int m = n; m >= -n; m--) {
            a_[RDIM][n][m] *= prefact;
        }
    }

    cout << "Calculated particle volume = " << vol << endl;
    cout.flush();
    
    // Now we just need to write the SH coefficients to a file

    try {
        int result = createSHFile(outfilename);
    }
    catch (tsquare::FileException fex) {
        fex.printException();
    }

    return;

}

void tsquare::Star::setSurface(const int num, const int nmax, bool gaussian=false)
{
    using std::cout;
    using std::endl;

    std::string junk;
    Real theta,phi;
    std::vector<Real> rval;

    nmax_ = (nmax < nmax_from_anmfile_) ? nmax : nmax_from_anmfile_;
    maxroundness_.clear();
    minroundness_.clear();
    maxroundness_.resize(ALLROUNDNESSES,0.0);
    minroundness_.resize(ALLROUNDNESSES,1.0e9);
    cumroundness_.resize(ALLROUNDNESSES,0.0);
    xg_.clear();
    wg_.clear();
    if (num > 0) {
        try {
            ntheta_ = nphi_= num;
            if (gaussian) {
                isgaussian_ = true;
                if (xg_.size() != (num + 1)) {
                    roundness_.clear();
                    cumroundness_.clear();
                    surface_.clear();
                    r_.clear();
                    x_.clear();
                    y_.clear();
                    z_.clear();
                    pcurv_.clear();
                    ddd2_.clear();
                    normal_.clear();
                    volume_ = -1.0;
                    area_ = -1.0;
                    getGausspoints(num);
                    for (int i = 0; i < surface_.size(); i++) {
                        if (verbose_) {
                            cout << "Setting the radius vector for point... "
                                 << i << endl;
                            cout.flush();
                        }
                        setR(i);
                    }
                    rval.resize(r_.size(),0.0);
                    roundness_.resize(ALLROUNDNESSES,rval);
                    computeVolume();
                    computeArea();
                }
            } else {
                isgaussian_ = false;
                if (verbose_) cout << "Setting surface to Regular, num = "
                                   << num << ", r size = " << r_.size() << endl;
                if (xg_.size() != 0 || r_.size() != num) {
                    roundness_.clear();
                    cumroundness_.clear();
                    surface_.clear();
                    r_.clear();
                    x_.clear();
                    y_.clear();
                    z_.clear();
                    pcurv_.clear();
                    ddd2_.clear();
                    normal_.clear();
                    volume_ = -1.0;
                    area_ = -1.0;
                    wg_.resize(num+1,1.0);
                    for (int i = 0; i < num; i++) {
                        theta = ((Real)(i) * Pi) / ((Real)(num));
                        if (i == 0) theta = 0.001 * Pi;
                        if (i == num - 1) theta = 0.999 * Pi;
                        for (int j = 0; j < num; j++) {
                            phi = (2.0 * Pi * (Real)(j))/((Real)((num)));
                            if (verbose_) {
                                cout << "Adding (theta,phi) = ("
                                     << theta << "," << phi
                                     << ") to surface... " << endl;
                                cout.flush();
                            }
                            addtoSurface(theta,phi);
                        }
                    }

                    for (int i = 0; i < surface_.size(); i++) {
                        if (verbose_) {
                            cout << "Setting the radius vector for point... " << i << endl;
                            cout.flush();
                        }
                        setR(i);
                    }
                    rval.resize(r_.size(),0.0);
                    roundness_.resize(ALLROUNDNESSES,rval);
                }
            }
        }
        catch (tsquare::EOBException ex) { throw ex; }
    } else {
        xg_.clear();
        wg_.clear();
        r_.clear();
        x_.clear();
        y_.clear();
        z_.clear();
        surface_.clear();
        pcurv_.clear();
        ddd2_.clear();
        normal_.clear();
        roundness_.clear();
        cumroundness_.clear();
        volume_ = -1.0;
        area_ = -1.0;
    }
    return;
}

int tsquare::Star::getGausspoints(int n)
{
    using std::cos;
    using std::pow;

    int m,i,j;
    Real eps = 1.0e-14;
    Real dn,e1,t,x0,den,d1,pnm1,pn,pnp1,dpn,d2pn,u,v,x1;
    
    if ((n % 2) != 0) {
        throw tsquare::DataException("Star","getGausspoints",
                            "Number of Gauss points must be even");
    }

    ntheta_ = nphi_ = n;

    std::vector<Real> x,w;

    xg_.clear();
    wg_.clear();
    x.clear();
    w.clear();
   
    xg_.resize(n + 1,0.0);
    wg_.resize(n + 1,0.0);
    x.resize(n + 1,0.0);   
    w.resize(n + 1,0.0);   

    /* Now populate the quadrature points */

    m = (n + 1)/2;
    dn = (Real)(n);
    e1 = dn * (dn + 1.0);
    for (i = 1; i <= m; i++) {
        t = ((4.0 * (Real)(i)) - 1.0) * Pi / (4.0 * dn + 2.0);
        x0 = (1.0 - (1.0-(1.0/dn)) / (8.0 * dn * dn)) * cos(t);
        do {
            getLegendre(n,&x0,&pn,&pnm1,&pnp1);
            den = 1.0 - (x0 * x0);
            d1 = dn * (pnm1 - (x0 * pn));
            dpn = d1 / den;
            d2pn = ((2.0 * x0 * dpn) - (e1 * pn)) / den;
            u = pn / dpn;
            v = d2pn / dpn;
            x1 = x0 - (u * (1.0 + (0.50 * u * v)));
            x0 = x1;
        } while ((x1 - x0) >= eps);
        x0 = x1;
        getLegendre(n,&x0,&pn,&pnm1,&pnp1);
        x[i] = x0;
        w[i] = (2.0 * (1.0 - (x0 * x0))) / pow((dn * pnm1),2.0);
    }

    if (m + m > n) x[m] = 0.0;
    
    for (i = 1; i <= m; i++) {
        xg_[i] = -x[i];
        xg_[i+m] = x[m+1-i];
        wg_[i] = w[i];
        wg_[i+m] = w[m+1-i];
    }
    
    std::vector<Real> op;
    op.resize(2,0.0);
    surface_.clear();
    for (i = 1; i <= n; i++) {
        op[0] = (0.5 * Pi * xg_[i]) + (0.5 * Pi);
        for (j = 1; j <= n; j++) {
            op[1] = (Pi * xg_[j]) + Pi;
            surface_.push_back(op);
        }
    }

    return(0);
}

void tsquare::Star::setR(const Real theta, const Real phi)
{
    using std::sin;
    using std::cos;

    std::complex<Real> r1 = a_[RDIM][0][0] * boost::math::spherical_harmonic(0,0,theta,phi);
    for (int n = 0; n <= nmax_; n++) {
        for (int m = n; m >= -n; m--) {
            r1 += (a_[RDIM][n][m] * boost::math::spherical_harmonic(n,m,theta,phi));
        }
    }

    // Check to see if this surface site is already listed.  If so,
    // then set the distance for that site. If not, add it now
 
    bool found = false;
    register int i;
    try {
        for (i = 0; (i < surface_.size() && !found); i++) {
            if (theta == surface_[i][0] && phi == surface_[i][1]) {
                found = true;
                r_.at(i) = real(r1);
                x_.at(i) = real(r1) * cos(phi) * sin(theta);
                y_.at(i) = real(r1) * sin(phi) * sin(theta);
                z_.at(i) = real(r1) * cos(theta);
            }
        }
    }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Star","setR","r_",r_.size(),i);
    }

    if (!found) {
        std::vector<Real> angles;
        angles.clear();
        angles.resize(2,0.0);
        angles[0] = theta;
        angles[1] = phi;
        surface_.push_back(angles);
        try {
            r_.push_back(real(r1));
            x_.push_back(real(r1) * cos(phi) * sin(theta));
            y_.push_back(real(r1) * sin(phi) * sin(theta));
            z_.push_back(real(r1) * cos(theta));
        }
        catch (std::out_of_range &oor) {
            throw tsquare::EOBException("Star","setR","r_",
                                  r_.size(),surface_.size()-1);
        }
    }
    return;
}

void tsquare::Star::setR(const unsigned int surfaceindex)
{
    using std::cos;
    using std::sin;

    Real theta,phi;
    try {
        theta = surface_.at(surfaceindex).at(0);
    }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Star","setR","surface_",
                           surface_.size(),surfaceindex);
    }

    try {
        phi = surface_.at(surfaceindex).at(1);
    }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Star","setR","surface_",
                           surface_.size(),surfaceindex);
    }

    std::complex<Real> r1(0.0,0.0);
    std::complex<Real> carg;
    Real rm;

    for (int n = 0; n <= nmax_; n++) {
        for (int m = n; m >= -n; m--) {
            rm = (Real)m;
            carg = std::complex<Real>(0.0,rm*phi);
            r1 += (a_[RDIM][n][m] * boost::math::spherical_harmonic(n,m,theta,phi));
        }
    }
 
    if (surfaceindex < r_.size()) {
        r_[surfaceindex] = real(r1);
        x_.at(surfaceindex) = real(r1) * cos(phi) * sin(theta);
        y_.at(surfaceindex) = real(r1) * sin(phi) * sin(theta);
        z_.at(surfaceindex) = real(r1) * cos(theta);
    } else {
        r_.push_back(real(r1));
        x_.push_back(real(r1) * cos(phi) * sin(theta));
        y_.push_back(real(r1) * sin(phi) * sin(theta));
        z_.push_back(real(r1) * cos(theta));
    }
    return;
}

void tsquare::Star::printProperties()
{
    using std::cout;
    using std::endl;

    cout << getType() << ":" << endl;
    cout << "    xlow,xhi = " << getXlo() << "," << getXhi() << endl;
    cout << "    ylow,yhi = " << getYlo() << "," << getYhi() << endl;
    cout << "    zlow,zhi = " << getZlo() << "," << getZhi() << endl;
    cout << "    volume = " << getVolume() << endl;
    cout << "    surface area = " << getArea() << endl;
    cout << "    normalized surface area = " << getNarea() << endl;
    cout << "    diameter = " << getDiam() << endl;
    cout << "    length, width, thickness = " << getLength() << ", " << getWidth() << ", " << getThickness() << endl;
    cout << "    normalized length = " << getNlength() << endl;
    cout << "    normalized width = " << getNwidth() << endl;
    cout << "    Itrace = " << getItrace() << endl;
    cout << "    Normalized Gaussian curvature = " << getNgc() << endl;
    cout << "    Max number of spherical harmonics = " << getNmax() << endl;
    cout.flush();
    cout << endl << endl;
    return;
}

Real tsquare::Star::getR(const Real theta, const Real phi)
{
//    using std::cos;
//    using std::sin;
//    using std::sqrt;

//  Below is for general SH expansion of particle
    std::complex<Real> r1(0.0,0.0);
    for (register int n = 0; n <= nmax_; n++) {
        for (int m = n; m >= -n; m--) {
            r1 += (a_[RDIM][n][m] * boost::math::spherical_harmonic(n,m,theta,phi));
        }
    }
    return real(r1);

//
//  Below is for a regular tetrahedron with side length sqrt(2)
//
//  Real t = theta;
//  Real p = phi;
//  Real rad = 1.0e9;
//  Real r1 = 1.0 / 2.0 / ((sin(t)*cos(p)) + (sin(t)*sin(p)) + cos(t));
//  if (r1 > 0.0 && r1 < rad) rad = r1;
//  Real r2 = 1.0 / 2.0 / ((sin(t)*cos(p)) - (sin(t)*sin(p)) - cos(t));
//  if (r2 > 0.0 && r2 < rad) rad = r2;
//  Real r3 = 1.0 / 2.0 / ((sin(t)*sin(p)) - (sin(t)*cos(p)) - cos(t));
//  if (r3 > 0.0 && r3 < rad) rad = r3;
//  Real r4 = 1.0 / 2.0 / (cos(t) - (sin(t)*sin(p)) - (sin(t)*cos(p)));
//  if (r4 > 0.0 && r4 < rad) rad = r4;
//  return rad;

//  Below is for a regular octahedron with side length 1/(2 sqrt 2)
//
//  Real t = theta;
//  Real p = phi;
//  Real rad = 1.0e9;
//  Real r1 = 1.0 / ((2.0*sin(t)*cos(p)) + (sqrt(2.0)*cos(t)));
//  Real r2 = 1.0 / ((2.0*sin(t)*sin(p)) + (sqrt(2.0)*cos(t)));
//  Real r3 = 1.0 / ((sqrt(2.0)*cos(t)) - (2.0*sin(t)*cos(p)));
//  Real r4 = 1.0 / ((sqrt(2.0)*cos(t)) - (2.0*sin(t)*sin(p)));
//  Real r5 = -r2;
//  Real r6 = 1.0 / ((2.0*sin(t)*cos(p)) - (sqrt(2.0)*cos(t)));
//  Real r7 = -r1;
//  Real r8 = 1.0 / ((2.0*sin(t)*sin(p)) - (sqrt(2.0)*cos(t)));
//  if (r1 > 0.0 && r1 < rad) rad = r1;
//  if (r2 > 0.0 && r2 < rad) rad = r2;
//  if (r3 > 0.0 && r3 < rad) rad = r3;
//  if (r4 > 0.0 && r4 < rad) rad = r4;
//  if (r5 > 0.0 && r5 < rad) rad = r5;
//  if (r6 > 0.0 && r6 < rad) rad = r6;
//  if (r7 > 0.0 && r7 < rad) rad = r7;
//  if (r8 > 0.0 && r8 < rad) rad = r8;

//  Below is for a regular dodecahedron

//  Real t = theta;
//  Real p = phi;
//  Real cp = cos(p);
//  Real sp = sin(p);
//  Real ct = cos(t);
//  Real st = sin(t);
//  Real rad = 1.0e9;
//  Real sr3 = sqrt(3.0);
//  Real sr5 = sqrt(5.0);
//  Real r1 = (sr3*(2.0 + sr5))/(sqrt(14.0 + 6.0*sr5)*ct - 4.0*cp*st);
//  Real r2 = -((sr3*(2.0 + sr5))/(sqrt(14.0 + 6.0*sr5)*ct - 4.0*cp*st));
//  Real r3 = (6.0*(2.0 + sr5))/
// ((-3.0 + sr5)*sqrt(6.0*(3.0+ sr5))*ct + 
//   2.0*(sqrt(6.0*(3.0 + sr5))*cp - 3.0*(1.0 + sr5)*sp)*st);
//  Real r4 = (sr3*(2.0 + sr5))/
// ((-1.0 + sr5)*ct + 2.0*sqrt(2.0*(3.0 + sr5))*cp*st);
//  Real r5 = (-3.0*(2.0 + sr5))/
// (sr3*(-1.0 + sr5)*ct - 
//   (sqrt(6.0*(3.0 + sr5))*cp + 3.0*(1.0 + sr5)*sp)*st);
//  Real r6 = (-3.0*(2.0 + sr5))/
// (sqrt(6.0*(7.0 + 3.0*sr5))*ct + 2.0*sr3*cp*st - 6.0*sp*st);
//  Real r7 = (-3.0*(2.0 + sr5))/(sr3*(3.0 + sr5)*ct + 2.0*(sr3*cp + 3.0*sp)*st);
//  Real r8 = (3.0*(2.0 + sr5))/
// (sr3*(-1.0 + sr5)*ct + 
//   (-(sqrt(6.0*(3.0 + sr5))*cp) + 3.0*(1.0 + sr5)*sp)*st);
//  Real r9 = ((2.0 + sr5)*sqrt(6.0/(3.0 + sr5)))/((-3.0 + sr5)*ct - 4.0*cp*st);
//  Real r10 = (3.0*(2.0 + sr5))/
// (sr3*(-1.0 + sr5)*ct - 
//   (sqrt(6.0*(3.0 + sr5))*cp + 3.0*(1.0 + sr5)*sp)*st);
//  Real r11 = (3.0*(2.0 + sr5))/(sr3*(3.0 + sr5)*ct + 2.0*(sr3*cp - 3.0*sp)*st);
//  Real r12 = (3.0*(2.0 + sr5))/(sr3*(3.0 + sr5)*ct + 2.0*(sr3*cp + 3.0*sp)*st);
//  if (r1 > 0.0 && r1 < rad) rad = r1;
//  if (r2 > 0.0 && r2 < rad) rad = r2;
//  if (r3 > 0.0 && r3 < rad) rad = r3;
//  if (r4 > 0.0 && r4 < rad) rad = r4;
//  if (r5 > 0.0 && r5 < rad) rad = r5;
//  if (r6 > 0.0 && r6 < rad) rad = r6;
//  if (r7 > 0.0 && r7 < rad) rad = r7;
//  if (r8 > 0.0 && r8 < rad) rad = r8;
//  if (r9 > 0.0 && r9 < rad) rad = r9;
//  if (r10 > 0.0 && r10 < rad) rad = r10;
//  if (r11 > 0.0 && r11 < rad) rad = r11;
//  if (r12 > 0.0 && r12 < rad) rad = r12;

//  Below is for a regular icosahedron

//  Real t = theta;
//  Real p = phi;
//  Real cp = cos(p);
//  Real sp = sin(p);
//  Real ct = cos(t);
//  Real st = sin(t);
//  Real sr5 = sqrt(5.0);
//  Real sr2 = sqrt(2.0);
//  Real rad = 1.0e9;

//  Real r1 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((5.0 - sr5)*ct + ((15.0 - 7.0*sr5)*cp + 
//      sqrt(250.0 - 110.0*sr5)*sp)*st);
//  Real r2 = sqrt((5.0 + sr5)/2.)/(2.*(ct + (-3.0 + sr5)*cp*st));
//  Real r3 = -(sqrt((5.0*(5.0 - sr5))/2.)/
//   ((-5.0 + sr5)*ct + ((-15.0 + 7.0*sr5)*cp + 
//        sqrt(250.0 - 110.0*sr5)*sp)*st));
//  Real r4 = sqrt(2.5)/(sqrt(5.0 - sr5)*ct + 
//   sr2*(sqrt(5.0 - 2.0*sr5)*cp + (-5.0 + 2.0*sr5)*sp)*st);
//  Real r5 = sqrt(2.5)/(sqrt(5.0 - sr5)*ct + 
//   sr2*(sqrt(5.0 - 2.0*sr5)*cp + (5.0 - 2.0*sr5)*sp)*st);
//  Real r6 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((-5.0 + sr5)*ct + ((-15.0 + 7.0*sr5)*cp - 
//      sqrt(250.0 - 110.0*sr5)*sp)*st);
//  Real r7 = -sqrt((5.0 + sr5)/2.)/(2.*(ct + (-3.0 + sr5)*cp*st));
//  Real r8 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((-5.0 + sr5)*ct + ((-15.0 + 7.0*sr5)*cp + 
//      sqrt(250.0 - 110.0*sr5)*sp)*st);
//  Real r9 = -(sqrt(2.5)/(sqrt(5.0 - sr5)*ct + 
//     sr2*(sqrt(5.0 - 2.0*sr5)*cp + (-5.0 + 2.0*sr5)*sp)*st));
//  Real r10 = sqrt(10.0)/(-2.0*sqrt(5.0 - sr5)*ct + 
//   ((-3.0 + sr5)*sqrt(5.0 + sr5)*cp + 
//      2.0*sr2*(-5.0 + 2.0*sr5)*sp)*st);
//  Real r11 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((-15.0 + 7.0*sr5)*ct + 2.0*((-5.0 + 2.0*sr5)*cp + 
//      sqrt(5.0*(5.0 - 2.0*sr5))*sp)*st);
//  Real r12 = sqrt((5.0*(5.0 + sr5))/2.)/
// (2.*((5.0 - 2.0*sr5)*ct + (-5.0 + sr5)*cp*st));
//  Real r13 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((-15.0 + 7.0*sr5)*ct + 2.0*((-5.0 + 2.0*sr5)*cp - 
//      sqrt(5.0*(5.0 - 2.0*sr5))*sp)*st);
//  Real r14 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((-15.0 + 7.0*sr5)*ct - ((-5.0 + sr5)*cp + 
//      sqrt(250.0 - 110.0*sr5)*sp)*st);
//  Real r15 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((-15.0 + 7.0*sr5)*ct + ((5.0 - sr5)*cp + 
//      sqrt(250.0 - 110.0*sr5)*sp)*st);
//  Real r16 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((15.0 - 7.0*sr5)*ct - 2.0*((-5.0 + 2.0*sr5)*cp + 
//      sqrt(5.0*(5.0 - 2.0*sr5))*sp)*st);
//  Real r17 = sqrt((5.0*(5.0 + sr5))/2.)/
// (2.*((-5.0 + 2.0*sr5)*ct - (-5.0 + sr5)*cp*st));
//  Real r18 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((15.0 - 7.0*sr5)*ct + 2.0*((5.0 - 2.0*sr5)*cp + 
//      sqrt(5.0*(5.0 - 2.0*sr5))*sp)*st);
//  Real r19 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((15.0 - 7.0*sr5)*ct + ((-5.0 + sr5)*cp + 
//      sqrt(250.0 - 110.0*sr5)*sp)*st);
//  Real r20 = sqrt((5.0*(5.0 - sr5))/2.)/
// ((15.0 - 7.0*sr5)*ct + ((-5.0 + sr5)*cp - 
//      sqrt(250.0 - 110.0*sr5)*sp)*st);
//
//  if (r1 > 0.0 && r1 < rad) rad = r1;
//  if (r2 > 0.0 && r2 < rad) rad = r2;
//  if (r3 > 0.0 && r3 < rad) rad = r3;
//  if (r4 > 0.0 && r4 < rad) rad = r4;
//  if (r5 > 0.0 && r5 < rad) rad = r5;
//  if (r6 > 0.0 && r6 < rad) rad = r6;
//  if (r7 > 0.0 && r7 < rad) rad = r7;
//  if (r8 > 0.0 && r8 < rad) rad = r8;
//  if (r9 > 0.0 && r9 < rad) rad = r9;
//  if (r10 > 0.0 && r10 < rad) rad = r10;
//  if (r11 > 0.0 && r11 < rad) rad = r11;
//  if (r12 > 0.0 && r12 < rad) rad = r12;
//  if (r13 > 0.0 && r13 < rad) rad = r13;
//  if (r14 > 0.0 && r14 < rad) rad = r14;
//  if (r15 > 0.0 && r15 < rad) rad = r15;
//  if (r16 > 0.0 && r16 < rad) rad = r16;
//  if (r17 > 0.0 && r17 < rad) rad = r17;
//  if (r18 > 0.0 && r18 < rad) rad = r18;
//  if (r19 > 0.0 && r19 < rad) rad = r19;
//  if (r20 > 0.0 && r20 < rad) rad = r20;
//
//  return rad;
} 

Real tsquare::Star::getR(const unsigned int surfaceindex)
{
    try {
        return r_.at(surfaceindex);
    }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Star","getR","r_",
                            r_.size(),1);
    }
}

void tsquare::Star::doRotate(Real alpha, Real beta, Real gamma)
{
    using std::cos;
    using std::sin;
    using std::sqrt;
    using std::pow;
    using std::exp;

    register int i;
    std::string buff;
    int klow,khigh;
    std::complex<Real> ddd,icmplx,aarg,garg;
    std::vector<tsquare::ExtendedVector> aa;
    Real realnum,abc,cosbeta,sinbeta,total,sgn;

    // Get all angles into their standard range

    sgn = 1.0;
    if (alpha < 0.0) sgn = -1.0;
    while (alpha < 0.0 || alpha >= 2.0 * Pi) {
        alpha -= (sgn * 2.0 * Pi);
    }

    sgn = 1.0;
    if (beta < 0.0) sgn = -1.0;
    while (beta < 0.0 || beta >= 2.0 * Pi) {
        beta -= (sgn * 2.0 * Pi);
    }

    sgn = 1.0;
    if (gamma < 0.0) sgn = -1.0;
    while (gamma < 0.0 || gamma >= 2.0 * Pi) {
        gamma -= (sgn * 2.0 * Pi);
    }

    // aa is for temporary storage
    aa.clear();
    aa.resize(nmax_+1);
    for (int i = 0; i <= nmax_; i++) {
       aa[i].clear();
       aa[i].resize(-i,i);
       for (int j = -i; j <= i; j++) {
           aa[i][j] = std::complex<Real>(0.0,0.0);
       }
    }

    
    cosbeta = cos(beta/2.0);
    sinbeta = sin(beta/2.0);
    if (cosbeta == 0.0) {
        beta += 1.0e-10;
        cosbeta = cos(beta/2.0);
    }
    if (sinbeta == 0.0) {
        beta += 1.0e-10;
        sinbeta = sin(beta/2.0);
    }

    std::vector<std::vector<std::vector<Real> > > dterm;
    std::vector<std::vector<Real> > dvec2;
    std::vector<Real> dvec1;

    // Initialize matrix of rotations about beta

    dterm.clear();
    dvec2.clear();
    dvec1.clear();
    dvec1.resize(2*(nmax_ + 1),0.0);
    for (i = 0; i < 2*(nmax_ + 1); i++) {
        dvec2.push_back(dvec1);
    }

    for (i = 0; i <= nmax_; i++) {
       dterm.push_back(dvec2);
    }

    // Now get all the dterm components at once using recursion relations
    // First, seed the recursion with n = 0 and n = 1 degree elements

    for (int n = 0; n <= 1; n++) {
        for (int m = -n; m <= n; m++) {
            for (int mp = -n; mp <= n; mp++) {
                realnum = sqrt(fac(n+mp)*fac(n-mp)/fac(n+m)/fac(n-m));
                ddd = std::complex<Real>(realnum,0.0);
                klow = std::max(0,m-mp);
                khigh = std::min(n-mp,n+m);
                total = 0.0;
                for (int k = klow; k <= khigh; k++) {
                    abc = pow(-1.0,k+mp-m);
                    abc *= (fac(n+m)/fac(k)/fac(n+m-k));
                    abc *= (fac(n-m)/fac(n-mp-k)/fac(mp+k-m));
                    total += abc * (pow(cosbeta,2*n+m-mp-2*k)) * (pow(sinbeta,2*k+mp-m));
                }
                dterm[n][getIndex(n,mp)][getIndex(n,m)] = total * realnum;
            }
        }
    }

    // Now we have seeded the recursion relations, build the rest using recursion

    Real term = 0.0;
    Real a,b,nb,c,d,nd,rn,rm,rmp;
    Real ss,cc,sc,cms;
    ss = sinbeta * sinbeta;
    cc = cosbeta * cosbeta;
    sc = sinbeta * cosbeta;
    cms = (cosbeta * cosbeta) - (sinbeta * sinbeta);
    for (int n = 2; n <= nmax_; n++) {
        rn = (Real)(n);
        for (int m = -n; m <= n; m++) {
            rm = (Real)(m);
            for (int mp = -n; mp <= n; mp++) {
                rmp = (Real)(mp);
                term = 0.0;
                if (mp > -n && mp < n) {
                    a = cms * sqrt((rn+rm)*(rn-rm)/(rn+rmp)/(rn-rmp));
                    b = sc * sqrt((rn+rm)*(rn+rm-1.0)/(rn+rmp)/(rn-rmp));
                    nb = -(sc * sqrt((rn-rm)*(rn-rm-1.0)/(rn+rmp)/(rn-rmp)));
                    term += a * dterm[n-1][getIndex(n-1,mp)][getIndex(n-1,m)];
                    term += b * dterm[n-1][getIndex(n-1,mp)][getIndex(n-1,m-1)];
                    term += nb * dterm[n-1][getIndex(n-1,mp)][getIndex(n-1,m+1)];
                    dterm[n][getIndex(n,mp)][getIndex(n,m)] = term;
                } else if (mp == -n) {
                    c = 2.0 * sc * sqrt((rn+rm)*(rn-rm)/(rn-rmp)/(rn-rmp-1.0));
                    d = ss * sqrt((rn+rm)*(rn+rm-1.0)/(rn-rmp)/(rn-rmp-1.0));
                    nd = cc * sqrt((rn-rm)*(rn-rm-1.0)/(rn-rmp)/(rn-rmp-1.0));
                    term += c * dterm[n-1][getIndex(n-1,mp+1)][getIndex(n-1,m)];
                    term += d * dterm[n-1][getIndex(n-1,mp+1)][getIndex(n-1,m-1)];
                    term += nd * dterm[n-1][getIndex(n-1,mp+1)][getIndex(n-1,m+1)];
                    dterm[n][getIndex(n,mp)][getIndex(n,m)] = term;
                } else {
                    c = -(2.0 * sc * sqrt((rn+rm)*(rn-rm)/(rn+rmp)/(rn+rmp-1.0)));
                    d = cc * sqrt((rn+rm)*(rn+rm-1.0)/(rn+rmp)/(rn+rmp-1.0));
                    nd = ss * sqrt((rn-rm)*(rn-rm-1.0)/(rn+rmp)/(rn+rmp-1.0));
                    term += c * dterm[n-1][getIndex(n-1,mp-1)][getIndex(n-1,m)];
                    term += d * dterm[n-1][getIndex(n-1,mp-1)][getIndex(n-1,m-1)];
                    term += nd * dterm[n-1][getIndex(n-1,mp-1)][getIndex(n-1,m+1)];
                    dterm[n][getIndex(n,mp)][getIndex(n,m)] = term;
                }
            }
        }
    }

    // Now we have all dterm rotation matrix elements stored

    for (int n = 0; n <= nmax_; n++) {
        for (int m = -n; m <= n; m++) {
            for (int mp = -n; mp <= n; mp++) {
                total = dterm[n][getIndex(n,mp)][getIndex(n,m)];
                ddd = std::complex<Real>(total,0.0);
                aarg = std::complex<Real>(0.0,mp*alpha);
                garg = std::complex<Real>(0.0,mp*gamma);
                icmplx = exp(garg) * (ddd * (exp(aarg)*(a_[RDIM][n][mp])));
                aa[n][m] += icmplx;
            }
        }
    }

    for (int n = 0; n <= nmax_; n++) {
        for (int m = -n; m <= n; m++) {
            a_[RDIM][n][m] = aa[n][m];
        }
    }

    doRotatesurfacepoints(alpha,beta,gamma);
}

void tsquare::Star::doRotatesurfacepoints(Real alpha, Real beta, Real gamma)
{
    using std::cos;
    using std::sin;
    using std::atan2;

    Real x,y,z,x_p,y_p,z_p;  // Letter p appended to a variable
                             // stands for transformed value
    Real theta,phi; 
    Real c1,c2,c3,s1,s2,s3; 
    Real a11,a12,a13;
    Real a21,a22,a23;
    Real a31,a32,a33;
    std::vector<Real> coords;

    c1 = cos(alpha);
    s1 = sin(alpha);
    c2 = cos(beta);
    s2 = sin(beta);
    c3 = cos(gamma);
    s3 = sin(gamma);

    a11 = (c3 * c2 * c1) - (s3 * s1);
    a12 = (c3 * c2 * s1) + (s3 * c1);
    a13 = (-c3 * s2);
    a21 = (-s3 * c2 * c1) - (c3 * s1);
    a22 = (c3 * c1) - (s3 * c2 * s1);
    a23 = (s3 * s2);
    a31 = (s2 * c1);
    a32 = (s2 * s1);
    a33 = c2;

    // Rotate the theta and phi surface coordinates
    for (register int i = 0; i < surface_.size(); i++) {
        theta = surface_[i][0];
        phi = surface_[i][1];
        x = cos(phi) * sin(theta);
        y = sin(phi) * sin(theta);
        z = cos(theta);
        x_p = (a11 * x) + (a12 * y) + (a13 * z);
        y_p = (a21 * x) + (a22 * y) + (a23 * z);
        z_p = (a31 * x) + (a32 * y) + (a33 * z);
        surface_[i][1] = atan2(y_p,x_p);
        surface_[i][0] = acos(z_p);
    }

    // Now rotate any actual points decorated on the surface
    for (register int i = 0; i < x_.size(); i++) {
        coords = getX(i);
        x = coords[XDIM];
        y = coords[YDIM];
        z = coords[ZDIM];
        x_p = (a11 * x) + (a12 * y) + (a13 * z);
        y_p = (a21 * x) + (a22 * y) + (a23 * z);
        z_p = (a31 * x) + (a32 * y) + (a33 * z);
        x_[i] = x_p;
        y_[i] = y_p;
        z_[i] = z_p;
    }
    return;
}

void tsquare::Star::doLanczos(const int power, int nstart)
{
    using std::cout;
    using std::endl;
    using std::sin;
    using std::pow;

    Real rval,ival,facbase,factor;
    rval = std::real(a_[RDIM][0][0]);
    ival = std::imag(a_[RDIM][0][0]);
    if (nstart < 0) {
        throw tsquare::DataException("Star","doLanczos","nstart less than zero");
    }
    if (nstart >= nmax_) {
        throw tsquare::DataException("Star","doLanczos","nstart greater than nmax");
    }
    if (power < 0) {
        throw tsquare::DataException("Star","undoLanczos","power < 0");
    }
    for (int n = nstart + 1; n <= nmax_; n++) {
        facbase = sin(Pi * ((Real)(n - nstart)) / ((Real)(nmax_ - nstart)))
		/ (Pi * ((Real)(n - nstart)) / ((Real)(nmax_ - nstart)));
        factor = pow(facbase,((Real)power));
        for (int m = -n; m <= n; m++) {
            if (verbose_) cout << "Reducing ringing... a(" << n << "," << m << ")" << endl;
            rval = std::real(a_[RDIM][n][m]) * factor;
            ival = std::imag(a_[RDIM][n][m]) * factor;
            a_[RDIM][n][m] = std::complex<Real>(rval,ival);
        }
    }
    return;
}

void tsquare::Star::undoLanczos(const int power, int nstart)
{
    using std::sin;
    using std::pow;

    Real rval,ival,facbase,factor;
    if (nstart < 0) {
        throw tsquare::DataException("Star","undoLanczos","nstart less than zero");
    }
    if (nstart >= nmax_) {
        throw tsquare::DataException("Star","undoLanczos","nstart greater than nmax");
    }
    if (power < 0) {
        throw tsquare::DataException("Star","undoLanczos","power < 0");
    }
    for (int n = nstart + 1; n <= nmax_; n++) {
        facbase = sin(Pi * ((Real)(n - nstart)) / ((Real)(nmax_ - nstart)))
		/ (Pi * ((Real)(n - nstart)) / ((Real)(nmax_ - nstart)));
        factor = pow(facbase,((Real)power));
        for (int m = -n; m <= n; m++) {
            rval = std::real(a_[RDIM][n][m]) / factor;
            ival = std::imag(a_[RDIM][n][m]) / factor;
            a_[RDIM][n][m] = std::complex<Real>(rval,ival);
        }
    }
    return;
}

void tsquare::Star::doHanning(int nstart)
{
    using std::cos;
    using std::cout;
    using std::endl;

    Real rval,ival,factor;
    if (nstart < 0) {
        throw tsquare::DataException("Star","doHanning","nstart less than zero");
    }
    if (nstart >= nmax_) {
        throw tsquare::DataException("Star","doHanning","nstart greater than nmax");
    }
    for (int n = nstart + 1; n <= nmax_; n++) {
        factor = 0.5 * (1.0 + cos(Pi * ((Real)(n - nstart)) / ((Real)(nmax_ - nstart))));
        for (int m = -n; m <= n; m++) {
            if (verbose_) cout << "Reducing ringing... a(" << n << "," << m << ")" << endl;
            rval = std::real(a_[RDIM][n][m]) * factor;
            ival = std::imag(a_[RDIM][n][m]) * factor;
            a_[RDIM][n][m] = std::complex<Real>(rval,ival);
        }
    }
    return;
}

void tsquare::Star::undoHanning(int nstart)
{
    using std::cos;

    Real rval,ival,factor;
    if (nstart < 0) nstart = 0;
    if (nstart >= nmax_) nstart = nmax_ - 1;
    for (int n = nstart + 1; n <= nmax_; n++) {
        factor = 0.5 * (1.0 + cos(Pi * ((Real)(n - nstart)) / ((Real)(nmax_ - nstart))));
        for (int m = -n; m <= n; m++) {
            rval = std::real(a_[RDIM][n][m]) / factor;
            ival = std::imag(a_[RDIM][n][m]) / factor;
            a_[RDIM][n][m] = std::complex<Real>(rval,ival);
        }
    }
    return;
}

void tsquare::Star::computeRoundness(const unsigned int method)
{
    using std::cout;
    using std::endl;
    using std::sqrt;
    using std::abs;
    using std::fabs;

    Real round,sa,wsa,rins,kins;
    Real factor = 0.5 * Pi * Pi;
    Real theta,phi;
    Real dotprod,k1,k2,kmax,curvterm,xmag,weight;
    Real f1,tmp,tmp1;
    std::vector<Real> X;
    std::vector<Real> rval;
    
    sa = area_;
    wsa = 0.0;
    round = 0.0;

    unsigned int count = 0;
    switch (method) {
      case DOTPRODUCT:  // Dot product definition
        try {
            for (register int i = 1; i <= ntheta_; i++) {
              for (register int j = 1; j <= nphi_; j++) {
                theta = surface_[count][0];
                phi = surface_[count][1];
                X = getX(count);
                xmag = sqrt(X[XDIM]*X[XDIM] + X[YDIM]*X[YDIM] + X[ZDIM]*X[ZDIM]);
                if (xmag > 0.0) {
                    dotprod = ((X[XDIM]/xmag) * normal_.at(count)[XDIM]);
                    dotprod += ((X[YDIM]/xmag) * normal_.at(count)[YDIM]);
                    dotprod += ((X[ZDIM]/xmag) * normal_.at(count)[ZDIM]);
                    f1 = abs(dotprod) * getR(theta,phi) * ddd2_.at(count);
                    roundness_.at(DOTPRODUCT).at(count) = abs(dotprod);
                    if (abs(dotprod) < minroundness_.at(DOTPRODUCT))
                        minroundness_.at(DOTPRODUCT) = abs(dotprod);
                    if (abs(dotprod) > maxroundness_.at(DOTPRODUCT))
                        maxroundness_.at(DOTPRODUCT) = abs(dotprod);
                    round += (wg_.at(i) * wg_.at(j) * f1);
                } else {
                    throw tsquare::DataException("Star","computeRoundness(DOTPRODUCT)",
                                    "Found zero position vector X.  Skipping.");
                }
                count++;
              }
            }
            if (sa > 0.0) {     
                round *= (factor/sa);
            } else {
                throw tsquare::FloatException("Star","computeRoundness(DOTPRODUCT)",
                                     "Divide by zero surface area");
            }
            cumroundness_.at(DOTPRODUCT) = round;
        }
        catch (std::out_of_range &oor) {
            throw tsquare::EOBException("Star","computeRoundness(DOTPRODUCT)","hard to say",
                            normal_.size(),count);
        }
        catch (tsquare::FloatException flex) { throw flex; }
        catch (tsquare::EOBException ex) { throw ex; }
        catch (tsquare::DataException dex) { throw dex; }
        break;
      case WADELL: // Max curvature --- analog of Wadell
        try {
            computePrincipalcurvatures();
            rins = inscribedsphere_.getRadius();
            if (fabs(rins) > 0.0) {
                kins = 1.0/rins;
            } else {
                throw tsquare::DataException("Star","computeRoundness(WADELL)",
                                  "Inscribed sphere radius of zero");
            }
            for (register int i = 1; i <= ntheta_; i++) {
              for (register int j = 1; j <= nphi_; j++) {
                theta = surface_.at(count)[0];
                phi = surface_.at(count)[1];
                k1 = abs(getPrincipalcurvature(0,count));
                k2 = abs(getPrincipalcurvature(1,count));

                kmax = (k1 > k2) ? k1 : k2;
                // kmax = 0.5 * (k1 + k2);
                kmax = (kmax > 1.0e-9) ? kmax : 1.0e-9;
                weight = (kmax < kins) ? 0.0 : 1.0;
                curvterm = (kins * weight / kmax);

                tmp = getR(theta,phi) * ddd2_.at(count) * curvterm;
                tmp1 = getR(theta,phi) * ddd2_.at(count) * weight;

                roundness_.at(WADELL).at(count) = (weight < 0.5) ? 1.0 : curvterm;

                if (roundness_.at(WADELL).at(count) < minroundness_.at(WADELL))
                    minroundness_.at(WADELL) = roundness_.at(WADELL).at(count);
                if (roundness_.at(WADELL).at(count) > maxroundness_.at(WADELL))
                    maxroundness_.at(WADELL) = roundness_.at(WADELL).at(count);
                round += (wg_.at(i) * wg_.at(j) * tmp);
                wsa += (wg_.at(i) * wg_.at(j) * tmp1);
                count++;
              }
            }
            if (wsa > 0.0) {
                cumroundness_.at(WADELL) = round / wsa;
             
            } else {
                throw tsquare::FloatException("Star","computeRoundness(WADELL)",
                                     "Divide by zero wsa");
            }
        }
        catch (std::out_of_range &oor) {
            throw tsquare::EOBException("Star","computeRoundness(WADELL)","hard to say",
                            surface_.size(),count);
        }
        catch (tsquare::FloatException flex) { throw flex; }
        catch (tsquare::EOBException ex) { throw ex; }
        catch (tsquare::DataException dex) { throw dex; }
        break;
      case ALLROUNDNESSES:  // Calculate all roundnesses at once
        try {
        	rins = 0.0;
            rins = inscribedsphere_.getRadius();
            if (fabs(rins) > 0.0) {
                kins = 1.0/rins;
            } else {
                throw tsquare::DataException("Star","computeRoundness(ALL)",
                                  "Inscribed sphere radius of zero");
            }
            computePrincipalcurvatures();
            for (register int i = 0; i < surface_.size(); i++) {
                theta = surface_[i][0];
                phi = surface_[i][1];
                if (getPrincipalcurvature(0,i) > getPrincipalcurvature(1,i)) {
                    k1 = getPrincipalcurvature(0,i);
                    k2 = getPrincipalcurvature(1,i);
                } else {
                    k1 = getPrincipalcurvature(1,i);
                    k2 = getPrincipalcurvature(0,i);
                }
            }
            rval.clear();
            rval.resize(ALLROUNDNESSES,0.0);
            for (register int i = 1; i <= ntheta_; i++) {
              for (register int j = 1; j <= nphi_; j++) {
                theta = surface_.at(count)[0];
                phi = surface_.at(count)[1];

                // dot product part

                X = getX(count);
                xmag = sqrt((X[XDIM]*X[XDIM])+(X[YDIM]*X[YDIM])+(X[ZDIM]*X[ZDIM]));
                if (xmag > 0.0) {
                    dotprod = ((X[XDIM]/xmag) * normal_.at(count)[XDIM]);
                    dotprod += ((X[YDIM]/xmag) * normal_.at(count)[YDIM]);
                    dotprod += ((X[ZDIM]/xmag) * normal_.at(count)[ZDIM]);
                    f1 = abs(dotprod) * getR(theta,phi) * ddd2_.at(count);
                    roundness_.at(DOTPRODUCT).at(count) = abs(dotprod);
                    if (abs(dotprod) < minroundness_.at(DOTPRODUCT))
                        minroundness_.at(DOTPRODUCT) = abs(dotprod);
                    if (abs(dotprod) > maxroundness_.at(DOTPRODUCT))
                        maxroundness_.at(DOTPRODUCT) = abs(dotprod);
                    rval[DOTPRODUCT] += (wg_.at(i) * wg_.at(j) * f1);
                } else {
                    tsquare::DataException dex("Star","computeRoundness(ALL)",
                                    "Found zero position vector X.  Skipping.");
                    dex.printException();
                }

                // Wadell part

                k1 = abs(getPrincipalcurvature(0,count));
                k2 = abs(getPrincipalcurvature(1,count));

                kmax = (k1 > k2) ? k1 : k2;
                // kmax = 0.5 * (k1 + k2);
                kmax = (kmax > 1.0e-9) ? kmax : 1.0e-9;
                weight = (kmax < kins) ? 0.0 : 1.0;
                curvterm = (kins * weight / kmax);

                tmp = getR(theta,phi) * ddd2_.at(count) * curvterm;
                tmp1 = getR(theta,phi) * ddd2_.at(count) * weight;

                roundness_.at(WADELL).at(count) = (weight < 0.5) ? 1.0 : curvterm;

                if (roundness_.at(WADELL).at(count) < minroundness_.at(WADELL))
                    minroundness_.at(WADELL) = roundness_.at(WADELL).at(count);
                if (roundness_.at(WADELL).at(count) > maxroundness_.at(WADELL))
                    maxroundness_.at(WADELL) = roundness_.at(WADELL).at(count);
                round += (wg_.at(i) * wg_.at(j) * tmp);
                wsa += (wg_.at(i) * wg_.at(j) * tmp1);
                cout << "(" << theta*Rad2Deg << "," << phi*Rad2Deg << ") :  k1 = " << k1
                     << ", k2 = " << k2 << ", curvterm = " << curvterm << ", wsa = " << wsa << endl;
                count++;
              }
            }
            if (verbose_) cout << "Calculating cumulative roundness measures... ";
            if (sa > 0.0) {
                cumroundness_[DOTPRODUCT] = rval[DOTPRODUCT] * factor / sa;
            } else {
                throw tsquare::FloatException("Star","computeRoundness(ALL)",
                                     "Divide by zero sa");
            }
            if (wsa > 0.0) {
                cumroundness_[WADELL] = round / wsa;
            } else {
                throw tsquare::FloatException("Star","computeRoundness(ALL)",
                                     "Divide by zero wsa");
            }
            if (verbose_) cout << "Done!" << endl;
        }
        catch (std::out_of_range &oor) {
            throw tsquare::EOBException("Star","computeRoundness(ALL)","hard to say",
                                surface_.size(),count);
        }
        catch (tsquare::FloatException flex) { throw flex; }
        catch (tsquare::EOBException ex) { throw ex; }
        catch (tsquare::DataException dex) { throw dex; }

        break;
      default:
        throw tsquare::DataException("Star","computeRoundness",
                              "Unrecognized roundness method");
    }

    return;
}

Real tsquare::Star::getRp(Real theta, const Real phi)
{
    std::complex<Real> rp(0.0,0.0);
    std::complex<Real> i(0.0,1.0);
    std::complex<Real> cm;
    for (int n = 0; n <= nmax_; n++) {
        for (int m = n; m >= -n; m--) {
            cm = std::complex<Real>((Real)m,0.0);
            rp += (i * cm * a_[RDIM][n][m] * boost::math::spherical_harmonic(n,m,theta,phi));
        }
    }
    return (real(rp));
}

Real tsquare::Star::getRp(const unsigned int surfindex)
{
    std::complex<Real> rp(0.0,0.0);
    std::complex<Real> i(0.0,1.0);
    std::complex<Real> cm;
    Real theta,phi;
    try {
        theta = surface_.at(surfindex)[0];
        phi = surface_.at(surfindex)[1];
        for (int n = 0; n <= nmax_; n++) {
            for (int m = n; m >= -n; m--) {
                cm = std::complex<Real>((Real)m,0.0);
                rp += (i * cm * a_[RDIM][n][m] * boost::math::spherical_harmonic(n,m,theta,phi));
            }
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getRp","surface_",
                        surface_.size(),surfindex);
        ex.printException();
    }
    return (real(rp));
}

Real tsquare::Star::getRpp(Real theta, const Real phi)
{
    std::complex<Real> rpp(0.0,0.0);
    Real rm;
    for (int n = 0; n <= nmax_; n++) {
        for (int m = n; m >= -n; m--) {
            rm = (Real)m;
            rpp -= (rm * rm * a_[RDIM][n][m] * boost::math::spherical_harmonic(n,m,theta,phi));
        }
    }
    return (real(rpp));
}

Real tsquare::Star::getRpp(const unsigned int surfindex)
{
    std::complex<Real> rpp(0.0,0.0);
    Real rm;
    Real theta,phi;
    try {
        theta = surface_.at(surfindex)[0];
        phi = surface_.at(surfindex)[1];
        for (int n = 0; n <= nmax_; n++) {
            for (int m = n; m >= -n; m--) {
                rm = (Real)m;
                rpp -= (rm * rm * a_[RDIM][n][m] * boost::math::spherical_harmonic(n,m,theta,phi));
            }
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getX","surface_",
                        surface_.size(),surfindex);
        ex.printException();
    }
    return (real(rpp));
}

Real tsquare::Star::getRt(Real theta, const Real phi)
{
    using std::sin;
    using std::cos;
    using std::exp;

    Real rn,rm,term;
    std::complex<Real> rt(0.0,0.0);
    std::complex<Real> carg;
    if (sin(theta) == 0.0) theta += 1.0e-5;
    Real st = sin(theta);
    Real u = cos(theta);
    for (int n = 0; n <= nmax_; n++) {
        rn = (Real)n;
        for (int m = n; m >= -n; m--) {
            rm = (Real)m;
            carg = std::complex<Real>(0.0,rm*phi);
            term = ((rn + 1.0) * u * getLegendre(n,m,u));
            term -= ((rn - rm + 1.0) * getLegendre(n+1,m,u));
            rt -= (a_[RDIM][n][m] * getFnm(rn,rm) * exp(carg) * term / st);
        }
    }
    return (real(rt));
}

Real tsquare::Star::getRt(const unsigned int surfindex)
{
    using std::sin;
    using std::cos;
    using std::exp;

    Real rn,rm,term;
    std::complex<Real> rt(0.0,0.0);
    std::complex<Real> carg;
    Real theta,phi;
    try {
        theta = surface_.at(surfindex)[0];
        phi = surface_.at(surfindex)[1];
        if (sin(theta) == 0.0) theta += 1.0e-5;
        Real st = sin(theta);
        Real u = cos(theta);
        for (int n = 0; n <= nmax_; n++) {
            rn = (Real)n;
            for (int m = n; m >= -n; m--) {
                rm = (Real)m;
                carg = std::complex<Real>(0.0,rm*phi);
                term = ((rn + 1.0) * u * getLegendre(n,m,u));
                term -= ((rn - rm + 1.0) * getLegendre(n+1,m,u));
                rt -= (a_[RDIM][n][m] * getFnm(rn,rm) * exp(carg) * term / st);
            }
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getRt","surface_",
                        surface_.size(),surfindex);
        ex.printException();
    }
    return (real(rt));
}

Real tsquare::Star::getRtt(Real theta, const Real phi)
{
    using std::sin;
    using std::cos;
    using std::exp;

    Real rn,rm,term1,term2,term3;
    std::complex<Real> rtt(0.0,0.0);
    std::complex<Real> carg;
    if (sin(theta) == 0.0) theta += 1.0e-5;
    Real st = sin(theta);
    Real u = cos(theta);

    for (int n = 0; n <= nmax_; n++) {
        rn = (Real)n;
        for (int m = n; m >= -n; m--) {
            rm = (Real)m;
            term1 = (rn + 1.0 + ((rn + 1.0) * (rn + 1) * u * u)) * getLegendre(n,m,u);
            term2 = (2.0 * u * (rn - rm + 1.0) * (rn + 2.0)) * getLegendre(n+1,m,u);
            term3 = ((rn - rm + 1.0) * (rn - rm + 2.0)) * getLegendre(n+2,m,u);
            carg = std::complex<Real>(0.0,rm*phi);
            rtt += (a_[RDIM][n][m] * getFnm(rn,rm) * exp(carg) * (term1 - term2 + term3)/st/st);
        }
    }
    return (real(rtt));
}

Real tsquare::Star::getRtt(const unsigned int surfindex)
{
    using std::sin;
    using std::cos;
    using std::exp;

    Real rn,rm,term1,term2,term3;
    std::complex<Real> rtt(0.0,0.0);
    std::complex<Real> carg;
    Real theta,phi;
    try {
        theta = surface_.at(surfindex)[0];
        phi = surface_.at(surfindex)[1];
        if (sin(theta) == 0.0) theta += 1.0e-5;
        Real st = sin(theta);
        Real u = cos(theta);

        for (int n = 0; n <= nmax_; n++) {
            rn = (Real)n;
            for (int m = n; m >= -n; m--) {
                rm = (Real)m;
                term1 = (rn + 1.0 + ((rn + 1.0) * (rn + 1) * u * u)) * getLegendre(n,m,u);
                term2 = (2.0 * u * (rn - rm + 1.0) * (rn + 2.0)) * getLegendre(n+1,m,u);
                term3 = ((rn - rm + 1.0) * (rn - rm + 2.0)) * getLegendre(n+2,m,u);
                carg = std::complex<Real>(0.0,rm*phi);
                rtt += (a_[RDIM][n][m] * getFnm(rn,rm) * exp(carg) * (term1 - term2 + term3)/st/st);
            }
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getRtt","surface_",
                        surface_.size(),surfindex);
        ex.printException();
    }
    return (real(rtt));
}

Real tsquare::Star::getRtp(Real theta, const Real phi)
{
    using std::sin;
    using std::cos;
    using std::exp;

    Real rn,rm,term;
    std::complex<Real> rtp(0.0,0.0);
    std::complex<Real> carg,im;
    if (sin(theta) == 0.0) theta += 1.0e-5;
    Real st = sin(theta);
    Real u = cos(theta);
    for (int n = 0; n <= nmax_; n++) {
        rn = (Real)n;
        for (int m = n; m >= -n; m--) {
            rm = (Real)m;
            carg = std::complex<Real>(0.0,rm*(phi));
            im = std::complex<Real>(0.0,rm);
            term = ((rn + 1.0) * u * getLegendre(n,m,u));
            term -= ((rn - rm + 1.0) * getLegendre(n+1,m,u));
            rtp -= (im * a_[RDIM][n][m] * getFnm(rn,rm) * exp(carg) * term / st);
        }
    }
    return (real(rtp));
}

Real tsquare::Star::getRtp(const unsigned int surfindex)
{
    using std::sin;
    using std::cos;
    using std::exp;

    Real rn,rm,term;
    std::complex<Real> rtp(0.0,0.0);
    std::complex<Real> carg,im;
    Real theta,phi;
    try {
        theta = surface_.at(surfindex)[0];
        phi = surface_.at(surfindex)[1];
        if (sin(theta) == 0.0) theta += 1.0e-5;
        Real st = sin(theta);
        Real u = cos(theta);
        for (int n = 0; n <= nmax_; n++) {
            rn = (Real)n;
            for (int m = n; m >= -n; m--) {
                rm = (Real)m;
                carg = std::complex<Real>(0.0,rm*(phi));
                im = std::complex<Real>(0.0,rm);
                term = ((rn + 1.0) * u * getLegendre(n,m,u));
                term -= ((rn - rm + 1.0) * getLegendre(n+1,m,u));
                rtp -= (im * a_[RDIM][n][m] * getFnm(rn,rm) * exp(carg) * term / st);
            }
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getRtp","surface_",
                        surface_.size(),surfindex);
        ex.printException();
    }
    return (real(rtp));
}

Real tsquare::Star::getX(const unsigned int index, const Real t, const Real p)
{
    using std::sin;
    using std::cos;

    Real x = 0.0;
    Real st = sin(t);
    switch (index) {
        case XDIM:
            x = (getR(t,p) * cos(p) * st);
            break;
        case YDIM:
            x = (getR(t,p) * sin(p) * st);
            break;
        case ZDIM:
            x = (getR(t,p) * cos(t));
            break;
        default:
            tsquare::DataException dex("Star","getX","Unrecognized index");
            dex.printException();
            break;
    }
    return x;
}

std::vector<Real> tsquare::Star::getX(const Real t, const Real p)
{
    using std::sin;
    using std::cos;

    Real st = sin(t);
    Real r = getR(t,p);
    std::vector<Real> xp;
    xp.clear();
    xp.resize(3,0.0);
    xp[XDIM] = (r * cos(p) * st);
    xp[YDIM] = (r * sin(p) * st);
    xp[ZDIM] = (r * cos(t));
    return xp;
}

std::vector<Real> tsquare::Star::getX(const unsigned int surfindex)
{
    std::vector<Real> xp;
    xp.clear();
    xp.resize(3,0.0);
    try {
        xp[XDIM] = x_.at(surfindex);
        xp[YDIM] = y_.at(surfindex);
        xp[ZDIM] = z_.at(surfindex);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getX","surface_",
                        surface_.size(),surfindex);
        ex.printException();
    }
    return xp;
}

Real tsquare::Star::getXp(const unsigned int index, const Real t, const Real p)
{
    using std::sin;
    using std::cos;

    Real dxdp = 0.0;
    Real st = sin(t);
    switch (index) {
        case XDIM:
            dxdp = (getRp(t,p) * st * cos(p)) - (getR(t,p) * st * sin(p));
            break;
        case YDIM:
            dxdp = (getRp(t,p) * st * sin(p)) + (getR(t,p) * st * cos(p));
            break;
        case ZDIM:
            dxdp = (getRp(t,p) * cos(t));
            break;
        default:
            tsquare::DataException dex("Star","getXp","Unrecognized index");
            dex.printException();
            break;
    }
    return dxdp;
}

std::vector<Real> tsquare::Star::getXp(const Real t, const Real p)
{
    std::vector<Real> xp;
    xp.clear();
    xp.resize(3,0.0);
    xp[XDIM] = getXp(XDIM,t,p);
    xp[YDIM] = getXp(YDIM,t,p);
    xp[ZDIM] = getXp(ZDIM,t,p);
    return xp;
}

Real tsquare::Star::getXp(const unsigned int index, const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,st;
    Real dxdp = 0.0;
    try {
        t = surface_.at(si)[0];
        p = surface_.at(si)[1];
        st = sin(t);
        switch (index) {
            case XDIM:
                dxdp = (getRp(si) * st * cos(p)) - (getR(si) * st * sin(p));
                break;
            case YDIM:
                dxdp = (getRp(si) * st * sin(p)) + (getR(si) * st * cos(p));
                break;
            case ZDIM:
                dxdp = (getRp(si) * cos(t));
                break;
            default:
                tsquare::DataException dex("Star","getXp","Unrecognized index");
                dex.printException();
                break;
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getX","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return dxdp;
}

std::vector<Real> tsquare::Star::getXp(const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,st,sp,cp,ct,rp,r;
    std::vector<Real> xp;
    xp.clear();
    xp.resize(3,0.0);
    try {
        t = surface_.at(si)[0];
        p = surface_.at(si)[1];
        st = sin(t);
        sp = sin(p);
        cp = cos(p);
        ct = cos(t);
        rp = getRp(t,p);
        r = getR(si);
        xp[XDIM] = (rp * st * cp) - (r * st * sp);
        xp[YDIM] = (rp * st * sp) + (r * st * cp);
        xp[ZDIM] = (rp * ct);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getXp","surface_",
                        surface_.size(),si);
        ex.printException();
        xp.clear();
        xp.resize(3,0.0);
    }
    return xp;
}

Real tsquare::Star::getXt(const unsigned int index, const Real t, const Real p)
{
    using std::sin;
    using std::cos;

    Real dxdt = 0.0;
    Real cp,sp;
    switch (index) {
        case XDIM:
            cp = cos(p);
            dxdt = (getRt(t,p) * sin(t) * cp) - (getR(t,p) * cos(t) * cp);
            break;
        case YDIM:
            sp = sin(p);
            dxdt = (getRt(t,p) * sin(t) * sp) + (getR(t,p) * cos(t) * sp);
            break;
        case ZDIM:
            dxdt = (getRt(t,p) * cos(t)) - (getR(t,p) * sin(t));
            break;
        default:
            tsquare::DataException dex("Star","getXt","Unrecognized index");
            dex.printException();
            break;
    }
    return dxdt;
}

std::vector<Real> tsquare::Star::getXt(const Real t, const Real p)
{
    std::vector<Real> xt;
    xt.clear();
    xt.resize(3,0.0);
    xt[XDIM] = getXt(XDIM,t,p);
    xt[YDIM] = getXt(YDIM,t,p);
    xt[ZDIM] = getXt(ZDIM,t,p);
    return xt;
}

Real tsquare::Star::getXt(const unsigned int index, const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,cp,sp;
    Real dxdt = 0.0;
    try {
        t = surface_.at(si)[0];
        p = surface_.at(si)[1];
        switch (index) {
            case XDIM:
                cp = cos(p);
                dxdt = (getRt(si) * sin(t) * cp) - (getR(si) * cos(t) * cp);
                break;
            case YDIM:
                sp = sin(p);
                dxdt = (getRt(si) * sin(t) * sp) + (getR(si) * cos(t) * sp);
                break;
            case ZDIM:
                dxdt = (getRt(si) * cos(t)) - (getR(si) * sin(t));
                break;
            default:
                tsquare::DataException dex("Star","getXt","Unrecognized index");
                dex.printException();
                break;
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getXt","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return dxdt;
}

std::vector<Real> tsquare::Star::getXt(const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,st,sp,cp,ct,rt,r;
    std::vector<Real> xt;
    xt.clear();
    xt.resize(3,0.0);
    try {
        t = surface_.at(si)[0];
        p = surface_.at(si)[1];
        st = sin(t);
        sp = sin(p);
        cp = cos(p);
        ct = cos(t);
        rt = getRt(t,p);
        r = getR(si);
        xt[XDIM] = (rt * st * cp) - (r * ct * cp);
        xt[YDIM] = (rt * st * sp) + (r * ct * sp);
        xt[ZDIM] = (rt * ct) - (r * st);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getXt","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return xt;
}

Real tsquare::Star::getXpp(const unsigned int index, const Real t, const Real p)
{
    using std::sin;
    using std::cos;

    Real dxdp = 0.0;
    Real st = sin(t);
    Real cp = cos(p);
    Real sp = sin(p);
    switch (index) {
        case XDIM:
            dxdp = st*(-(cp*getR(t,p)) - (2.0*sp*getRp(t,p)) + (cp*getRpp(t,p)));
            break;
        case YDIM:
            dxdp = st*(-(getR(t,p)*sp) + (2.0*cp*getRp(t,p)) + (sp*getRpp(t,p)));
            break;
        case ZDIM:
            dxdp = getRpp(t,p) * cos(t);
            break;
        default:
            tsquare::DataException dex("Star","getXpp","Unrecognized index");
            dex.printException();
            break;
    }
    return dxdp;
}

std::vector<Real> tsquare::Star::getXpp(const Real t, const Real p)
{
    std::vector<Real> xpp;
    xpp.clear();
    xpp.push_back(getXpp(XDIM,t,p));
    xpp.push_back(getXpp(YDIM,t,p));
    xpp.push_back(getXpp(ZDIM,t,p));
    return xpp;
}

Real tsquare::Star::getXpp(const unsigned int index, const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,st;
    Real dxdp = 0.0;
    try {
        t = surface_.at(si)[0];
        p = surface_.at(si)[1];
        st = sin(t);
        switch (index) {
            case XDIM:
                dxdp = ((getRpp(si) - getR(si)) * st * cos(p))
                       - (2.0 * getRp(si) * st * sin(p));
                break;
            case YDIM:
                dxdp = ((getRpp(si) - getR(si)) * st * sin(p))
                       + (2.0 * getRp(si) * st * cos(p));
                break;
            case ZDIM:
                dxdp = getRpp(si) * cos(t);
                break;
            default:
                tsquare::DataException dex("Star","getXpp","Unrecognized index");
                dex.printException();
                break;
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getXt","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return dxdp;
}

std::vector<Real> tsquare::Star::getXpp(const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,st,sp,cp,ct,rpp,rp,r;
    std::vector<Real> xpp;
    xpp.clear();
    xpp.resize(3,0.0);
    try {
        t = surface_[si][0];
        p = surface_[si][1];
        st = sin(t);
        sp = sin(p);
        cp = cos(p);
        ct = cos(t);
        rpp = getRpp(t,p);
        rp = getRp(t,p);
        r = getR(si);
        xpp[XDIM] = ((rpp - r) * st * cp) - (2.0 * rp * st * sp);
        xpp[YDIM] = ((rpp - r) * st * sp) + (2.0 * rp * st * cp);
        xpp[ZDIM] = rpp * ct;
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getXpp","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return xpp;
}

Real tsquare::Star::getXtt(const unsigned int index, const Real t, const Real p)
{
    using std::sin;
    using std::cos;

    Real dxdt = 0.0;
    Real cp = cos(p);
    Real st = sin(t);
    Real ct = cos(t);
    Real sp = sin(p);
    switch (index) {
        case XDIM:
            dxdt = cp*(-(getR(t,p)*st) + (2.0*ct*getRt(t,p)) + (st*getRtt(t,p)));
            break;
        case YDIM:
            dxdt = sp*(-(getR(t,p)*st) + (2.0*ct*getRt(t,p)) + (st*getRtt(t,p)));
            break;
        case ZDIM:
            dxdt = -(ct*getR(t,p)) - (2.0*st*getRt(t,p)) + (ct*getRtt(t,p));
            break;
        default:
            tsquare::DataException dex("Star","getXtt","Unrecognized index");
            dex.printException();
            break;
    }
    return dxdt;
}

std::vector<Real> tsquare::Star::getXtt(const Real t, const Real p)
{
    using std::sin;
    using std::cos;

    std::vector<Real> xtt;
    xtt.clear();
    xtt.resize(3,0.0);
    Real rtt = getRtt(t,p);
    Real rt = getRt(t,p);
    Real r = getR(t,p);
    Real st = sin(t);
    Real sp = sin(p);
    Real cp = cos(p);
    Real ct = cos(t);
    xtt[XDIM] = cp*(-(r*st) + (2.0*ct*rt) + (st*rtt));
    xtt[YDIM] = sp*(-(r*st) + (2.0*ct*rt) + (st*rtt));
    xtt[ZDIM] = -(ct*r) - (2.0*st*rt) + (ct*rtt);
    return xtt;
}

Real tsquare::Star::getXtt(const unsigned int index, const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,cp,sp;
    Real dxdt = 0.0;
    try {
        t = surface_.at(si)[0];
        p = surface_.at(si)[1];
        switch (index) {
            case XDIM:
                cp = cos(p);
                dxdt = ((getRtt(si) - getR(si)) * sin(t) * cp)
                       + (2.0 * getRt(si) * cos(t) * cp);
                break;
            case YDIM:
                sp = sin(p);
                dxdt = ((getRtt(si) - getR(si)) * sin(t) * sp)
                       + (2.0 * getRt(si) * cos(t) * sp);
                break;
            case ZDIM:
                dxdt = ((getRtt(si) - getR(si)) * cos(t))
                         - (2.0 * getRt(si) * sin(t));
                break;
            default:
                tsquare::DataException dex("Star","getXtt","Unrecognized index");
                dex.printException();
                break;
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getXtt","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return dxdt;
}

std::vector<Real> tsquare::Star::getXtt(const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,st,sp,cp,ct,rtt,rt,r;
    std::vector<Real> xtt;
    xtt.clear();
    xtt.resize(3,0.0);
    try {
        t = surface_.at(si)[0];
        p = surface_.at(si)[1];
        st = sin(t);
        sp = sin(p);
        cp = cos(p);
        ct = cos(t);
        rtt = getRtt(t,p);
        rt = getRt(t,p);
        r = getR(si);
        xtt[XDIM] = ((rtt - r) * st * cp) + (2.0 * rt * ct * cp);
        xtt[YDIM] = ((rtt - r) * st * sp) + (2.0 * rt * ct * sp);
        xtt[ZDIM] = ((rtt - r) * ct) - (2.0 * rt * st);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getXtt","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return xtt;
}

Real tsquare::Star::getXtp(const unsigned int index, const Real t, const Real p)
{
    using std::sin;
    using std::cos;

    Real dxdtp = 0.0;
    Real st = sin(t);
    Real sp = sin(p);
    Real ct = cos(t);
    Real cp = cos(p);
    switch (index) {
        case XDIM:
            dxdtp = -(ct*getR(t,p)*sp) + cp*ct*getRp(t,p)
			+ (st*(-(sp*getRt(t,p)) + cp*getRtp(t,p)));
            break;
        case YDIM:
            dxdtp = cp*ct*getR(t,p) + ct*sp*getRp(t,p) + (st*(cp*getRt(t,p) + sp*getRtp(t,p)));
            break;
        case ZDIM:
            dxdtp = -(st*getRp(t,p)) + ct*getRtp(t,p);
            break;
        default:
            tsquare::DataException dex("Star","getXtp","Unrecognized index");
            dex.printException();
            break;
    }
    return dxdtp;
}

std::vector<Real> tsquare::Star::getXtp(const Real t, const Real p)
{
    using std::sin;
    using std::cos;

    Real st = sin(t);
    Real sp = sin(p);
    Real cp = cos(p);
    Real ct = cos(t);
    Real rtp = getRtp(t,p);
    Real rt = getRt(t,p);
    Real rp = getRp(t,p);
    Real r = getR(t,p);
    std::vector<Real> xtp;
    xtp.clear();
    xtp.resize(3,0.0);
    xtp[XDIM] = (-(ct*r*sp) + cp*ct*rp + (st*(-(sp*rt) + cp*rtp)));
    xtp[YDIM] = (cp*ct*r + ct*sp*rp + (st*(cp*rt + sp*rtp)));
    xtp[ZDIM] = (-(st*rp) + ct*rtp);
    return xtp;
}

Real tsquare::Star::getXtp(const unsigned int index, const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,st,sp,ct,cp;
    Real dxdtp = 0.0;
    try {
        t = surface_.at(si)[0];
        p = surface_.at(si)[1];
        st = sin(t);
        sp = sin(p);
        ct = cos(t);
        cp = cos(p);
        switch (index) {
            case XDIM:
                dxdtp = (getRtp(si) * st * cp) - (getRt(si) * st * sp)
                        + (getRp(si) * ct * cp) - (getR(si) * ct * sp);
                break;
            case YDIM:
                dxdtp = (getRtp(si) * st * sp) + (getRt(si) * st * cp)
                        + (getRp(si) * ct * sp) + (getR(si) * ct * cp);
                break;
            case ZDIM:
                dxdtp = (getRtp(si) * ct) - (getRp(si) * st);
                break;
            default:
                tsquare::DataException dex("Star","getXtp","Unrecognized index");
                dex.printException();
                break;
        }
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getXtp","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return dxdtp;
}

std::vector<Real> tsquare::Star::getXtp(const unsigned int si)
{
    using std::sin;
    using std::cos;

    Real t,p,st,sp,cp,ct,rtp,rp,rt,r;
    std::vector<Real> xtp;
    xtp.clear();
    xtp.resize(3,0.0);
    try {
        t = surface_.at(si)[0];
        p = surface_.at(si)[1];
        st = sin(t);
        sp = sin(p);
        cp = cos(p);
        ct = cos(t);
        rtp = getRtp(t,p);
        rp = getRp(t,p);
        rt = getRt(t,p);
        r = getR(si);
        xtp[XDIM] = ((rtp * st * cp) - (rt * st * sp) + (rp * ct * cp) - (r * ct * sp));
        xtp[YDIM] = ((rtp * st * sp) + (rt * st * cp) + (rp * ct * sp) + (r * ct * cp));
        xtp[ZDIM] = ((rtp * ct) - (rp * st));
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getXtp","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return xtp;
}

Real tsquare::Star::getE(const Real theta, const Real phi)
{
    // using std::cout;
    // using std::endl;

    Real rt = getRt(theta,phi);
    Real r = getR(theta,phi);
    // cout << "--> Star::getE(" << theta * 180.0/Pi << "," << phi * 180.0/Pi << ")" << endl;
    Real E = (rt * rt) + (r * r);
    return E;
}

Real tsquare::Star::getE(const unsigned int si)
{
    using std::cout;
    using std::endl;

    Real theta,phi;
    Real E = 0.0;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        cout << "Star::getE(" << theta << "," << phi << ")" << endl;
        E = getE(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getE","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return E;
}

Real tsquare::Star::getF(const Real theta, const Real phi)
{
    using std::sin;
    using std::cos;

    Real st = sin(theta);
    Real ct = cos(theta);
    Real cp = cos(phi);
    Real F = (getRt(theta,phi) * getRp(theta,phi))
               * ((2.0 * st * st * cp * cp) + (ct * ct));
    return F;
}

Real tsquare::Star::getF(const unsigned int si)
{
    using std::cout;
    using std::endl;

    Real theta,phi;
    Real F = 0.0;
    try {
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        cout << "Star::getF(" << theta << "," << phi << ")" << endl;
        F = getF(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getF","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return F;
}

Real tsquare::Star::getG(const Real theta, const Real phi)
{
    using std::sin;

    Real st = sin(theta);
    Real rp = getRp(theta,phi);
    Real r = getR(theta,phi);
    Real G = (rp * rp) + (r * r * st * st);
    return G;
}

Real tsquare::Star::getG(const unsigned int si)
{
    using std::cout;
    using std::endl;

    Real theta,phi;
    Real G = 0.0;
    try {
        theta = surface_[si][0];
        phi = surface_[si][1];
        cout << "Star::getG(" << theta << "," << phi << ")" << endl;
        G = getG(theta,phi);
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getG","surface_",
                        surface_.size(),si);
        ex.printException();
    }
    return G;
}

Real tsquare::Star::getNormal(const unsigned int index, const unsigned int si)
{
    using std::sin;
    using std::cos;
    using std::sqrt;

    Real theta,phi,ct,st,cp,sp,val,denom,r,rp,rt;
    std::vector<Real> ncomp;
    ncomp.clear();
    ncomp.resize(3,0.0);
    try{
        theta = surface_.at(si)[0];
        phi = surface_.at(si)[1];
        ct = cos(theta);
        st = sin(theta);
        cp = cos(phi);
        sp = sin(phi);
        r = getR(si);
        rp = getRp(theta,phi);
        rt = getRt(theta,phi);

        ncomp[XDIM] = r * ((cp*r*st*st) + (sp*rp) - (cp*ct*st*rt));
        ncomp[YDIM] = r * ((r*sp*st*st) - (cp*rp) - (ct*sp*st*rt));
        ncomp[ZDIM] = r * st * ((ct*r) + (st*rt));

        denom = sqrt((ncomp[XDIM]*ncomp[XDIM]) +
                     (ncomp[YDIM]*ncomp[YDIM]) +
                     (ncomp[ZDIM]*ncomp[ZDIM]));
        if (denom > 0.0) {
            switch (index) {
                case XDIM:
                    val = ncomp[XDIM]/denom;
                    break;
                case YDIM:
                    val = ncomp[YDIM]/denom;
                    break;
                case ZDIM:
                    val = ncomp[ZDIM]/denom;
                    break;
                default:
                    tsquare::DataException dex("Star","getNormal","Unrecognized index");
                    dex.printException();
                    val = 0.0;
                    break;
            }
        } else {
            tsquare::DataException dex("Star","getNormal","Zero Unit normal vector");
            dex.printException();
            val = 0.0;
        }
    }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Star","getNormal","surface_",
                        surface_.size(),si);
    }
    return (val);
}

Real tsquare::Star::getNormal(const unsigned int index, const Real t, const Real p)
{
    using std::sin;
    using std::cos;
    using std::sqrt;

    // t is theta, and p is phi

    Real ct = cos(t);
    Real st = sin(t);
    Real cp = cos(p);
    Real sp = sin(p);

    Real val,denom;
    std::vector<Real> ncomp;
    ncomp.clear();
    ncomp.resize(3,0.0);
    ncomp[XDIM] = getR(t,p) * ((cos(p)*getR(t,p)*st*st)
		+ (sp*getRp(t,p)) - (cp*ct*st*getRt(t,p)));

    ncomp[YDIM] = getR(t,p) * ((getR(t,p)*sp*st*st)
		- (cp*getRp(t,p)) - (ct*sp*st*getRt(t,p)));

    ncomp[ZDIM] = getR(t,p) * st * ((ct*getR(t,p)) + (st*getRt(t,p)));

    denom = sqrt((ncomp[XDIM]*ncomp[XDIM]) + (ncomp[YDIM]*ncomp[YDIM])
                  + (ncomp[ZDIM]*ncomp[ZDIM]));
    if (denom > 0.0) {
        switch (index) {
            case XDIM:
                val = ncomp[XDIM]/denom;
                break;
            case YDIM:
                val = ncomp[YDIM]/denom;
                break;
            case ZDIM:
                val = ncomp[ZDIM]/denom;
                break;
            default:
                tsquare::DataException dex("Star","getNormal","Unrecognized index");
                dex.printException();
                val = 0.0;
                break;
        }
    } else {
        throw tsquare::DataException("Star","getNormal","Zero Unit normal vector");
    }
    return (val);
}

std::vector<Real> tsquare::Star::getNormal(const Real t, const Real p)
{
    using std::sin;
    using std::cos;
    using std::sqrt;

    std::vector<Real> nvec;
    nvec.clear();
    nvec.resize(3,0.0);

    // t is theta, and p is phi

    Real ct = cos(t);
    Real st = sin(t);
    Real cp = cos(p);
    Real sp = sin(p);
    Real r = getR(t,p);
    Real rp = getRp(t,p);
    Real rt = getRt(t,p);

    Real denom;
    std::vector<Real> ncomp;
    ncomp.clear();
    ncomp.resize(3,0.0);
    ncomp[XDIM] = r * ((cp*r*st*st) + (sp*rp) - (cp*ct*st*rt));
    ncomp[YDIM] = r * ((r*sp*st*st) - (cp*rp) - (ct*sp*st*rt));
    ncomp[ZDIM] = r * st * ((ct*r) + (st*rt));

    denom = sqrt((ncomp[XDIM]*ncomp[XDIM]) +
                 (ncomp[YDIM]*ncomp[YDIM]) +
                 (ncomp[ZDIM]*ncomp[ZDIM]));
    if (denom > 0.0) {
        nvec[XDIM] = ncomp[XDIM]/denom;
        nvec[YDIM] = ncomp[YDIM]/denom;
        nvec[ZDIM] = ncomp[ZDIM]/denom;
    } else {
        throw tsquare::DataException("Star","getNormal","Zero Unit normal vector");
    }
    return (nvec);
}

std::vector<Real> tsquare::Star::getNormal(const unsigned int si)
{
    try { return (normal_.at(si)); }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Star","getNormal","normal_",
                        normal_.size(),si);
    }
}

void tsquare::Star::computeVolume(void)
{
    using std::pow;
    using std::sin;

    Real theta,v1;
    Real factor = 0.5 * Pi * Pi;
    int count;
    volume_ = v1 = 0.0;

    count = 0;
    try {
        for (register int i = 1; i <= ntheta_; i++) {
            for (register int j = 1; j <= nphi_; j++) {
                theta = surface_.at(count)[0];
                v1 = sin(theta)/3.0;
                v1 *= pow(r_.at(count),3.0);
                v1 *= (wg_.at(i) * wg_.at(j));
                volume_ += v1;
                count++;
            }
        }
        volume_ *= factor;
        return;
    }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Star","computeVolume","surface_",
                        surface_.size(),count);
    }
}

void tsquare::Star::computeArea(void)
{
    using std::cout;
    using std::endl;
    using std::cos;
    using std::sin;
    using std::sqrt;

    Real theta,phi,sa,st,ct,sp,cp,r,rp,rt,ddd,f1,denom;
    Real factor = 0.5 * Pi * Pi;
    int count;
    std::vector<Real> ncomp;
    ncomp.resize(3,0.0);
    area_ = sa = 0.0;
    count = 0;

    try {
        for (register int i = 1; i <= ntheta_; i++) {
            if (verbose_) {
              cout << "    Area computation:  Latitude " << i << " out of " << ntheta_ << endl;
              cout.flush();
            }
            for (register int j = 1; j <= nphi_; j++) {
                theta = surface_.at(count)[0];
                phi = surface_.at(count)[1];
                st = sin(theta);
                ct = cos(theta);
                sp = sin(phi);
                cp = cos(phi);
                r = getR(count);
                rp = getRp(theta,phi);
                rt = getRt(theta,phi);
                ddd = (rp * rp) + (rt * rt * st * st) + (r * r * st * st);
                ddd2_.push_back(sqrt(ddd)); // Area of surface patch
                ncomp[XDIM] = r * ((cp*r*st*st) + (sp*rp) - (cp*ct*st*rt));
                ncomp[YDIM] = r * ((r*sp*st*st) - (cp*rp) - (ct*sp*st*rt));
                ncomp[ZDIM] = r * st * ((ct*r) + (st*rt));

                denom = sqrt((ncomp[XDIM]*ncomp[XDIM]) +
                             (ncomp[YDIM]*ncomp[YDIM]) +
                             (ncomp[ZDIM]*ncomp[ZDIM]));
                if (denom > 0.0) {
                    ncomp[XDIM] /= denom;
                    ncomp[YDIM] /= denom;
                    ncomp[ZDIM] /= denom;
                }
                normal_.push_back(ncomp);
                f1 = r * sqrt(ddd);
                sa += (wg_.at(i) * wg_.at(j) * f1);
                count++;
            }
        }
        area_ = sa * factor;
        return;
    }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Star","computeVolume","surface_",
                        surface_.size(),count);
    }
}

void tsquare::Star::computeCentroid(void)
{
    using std::sin;
    using std::cos;
    using std::pow;

    Real theta,phi,xcoord,ycoord,zcoord,x1,x2,x3,v1;
    Real factor = 0.5 * Pi * Pi;
    int count;
    volume_ = xcoord = ycoord = zcoord = v1 = 0.0;

    count = 0;
    try {
        for (register int i = 1; i <= ntheta_; i++) {
            for (register int j = 1; j <= nphi_; j++) {
                theta = surface_.at(count)[0];
                phi = surface_.at(count)[1];
                v1 = sin(theta)/3.0;
                v1 *= pow(r_.at(count),3.0);
                v1 *= (wg_.at(i) * wg_.at(j));
                x1 = v1 * r_.at(count) * sin(theta) * cos(phi);
                x2 = v1 * r_.at(count) * sin(theta) * sin(phi);
                x3 = v1 * r_.at(count) * cos(theta);
                volume_ += v1;
                xcoord += x1;
                ycoord += x2;
                zcoord += x3;
                count++;
            }
        }
        volume_ *= factor;
        xcoord *= (factor/volume_);
        ycoord *= (factor/volume_);
        zcoord *= (factor/volume_);
        centroid_[0] = xcoord;
        centroid_[1] = ycoord;
        centroid_[2] = zcoord;
        return;
    }
    catch (std::out_of_range &oor) {
        throw tsquare::EOBException("Star","computeVolume","surface_",
                        surface_.size(),count);
    }
}

void tsquare::Star::computeCoords(std::vector<std::vector<Real> > sc, std::vector<Real> &xc,
                                   std::vector<Real> &yc,
                                   std::vector<Real> &zc, std::vector<Real> &rc)
{
    using std::sin;
    using std::cos;

    Real theta,phi,r,x,y,z;

    int numpoints = sc.size();
    xc.clear();
    yc.clear();
    zc.clear();
    rc.clear();

    for (register int i = 0; i < numpoints; i++) {
        theta = sc[i][0];
        phi = sc[i][1];
        r = getR(theta,phi);
        rc.push_back(r);
        x = r * sin(theta) * cos(phi);
        y = r * sin(theta) * sin(phi);
        z = r * cos(theta);
        xc.push_back(x);
        yc.push_back(y);
        zc.push_back(z);
    }
}

Real tsquare::Star::computeLength(std::vector<std::vector<Real> > sc, std::vector<Real> xc,
                   std::vector<Real> yc, std::vector<Real> zc, std::vector<Real> &lvec,
                   Real &t1, Real &t2, Real &p1, Real &p2)
{
    using std::cout;
    using std::endl;
    using std::sqrt;
    using std::pow;

    int numpoints = (int)(sqrt((Real)xc.size()));
    Real theta1,phi1,theta2,phi2;
    Real x1,y1,z1,x2,y2,z2,xn,yn,zn,rm;
    Real tt1,tt2,pp1,pp2;
    Real rmax = 0.0;
    std::vector<Real> orderedpair;
    std::vector<std::vector<Real> > surface_fine,surface_comp;
    std::vector<Real> x_fine,y_fine,z_fine,r_fine;
    std::vector<Real> x_comp,y_comp,z_comp,r_comp;

    // Unrefined length first

    if ((xc.size() != yc.size()) || (xc.size() != zc.size())
         || (xc.size() != sc.size())) {
        throw tsquare::EOBException("Star","computeLength","coarse coordinates",
                           sc.size(),xc.size());
    }
    for (register int i = 0; i < xc.size(); i++) {
        theta1 = sc[i][0];
        phi1 = sc[i][1];
        x1 = xc[i];
        y1 = yc[i];
        z1 = zc[i];
        for (register int j = i; j < xc.size(); j++) {
            theta2 = sc[j][0];
            phi2 = sc[j][1];
            x2 = xc[j];
            y2 = yc[j];
            z2 = zc[j];

            rm = sqrt(pow(x1-x2,2.0)+pow(y1-y2,2.0)+pow(z1-z2,2.0));

            if (rmax < rm) {
                rmax = rm;
                t1 = theta1;
                p1 = phi1;
                t2 = theta2;
                p2 = phi2;

                // unit vector for length axis

                xn = (x2 - x1) / rmax;
                yn = (y2 - y1) / rmax;
                zn = (z2 - z1) / rmax;
            }
        }
    }

    lvec.clear();
    lvec.resize(3,0.0);
    lvec[0] = xn;
    lvec[1] = yn;
    lvec[2] = zn;

    if (verbose_) {
        cout << "Done comparing for unrefined lengths" << endl;
        cout.flush();
        cout << "Unrefined length = " << rmax << " at (" << t1 << ","
             << p1 << ") and (" << t2 << "," << p2 << ")" << endl;
        cout.flush();
    }

    // Now do the refined length

    surface_fine.clear();
    orderedpair.clear();
    orderedpair.resize(2,0.0);

    for (register int i = 1; i <= numpoints; i++) {
        theta1 = t1 + (((Real)(i - 10) * Pi) / ((Real)numpoints) / 10.0);
        orderedpair[0] = theta1;
        for (register int j = 1; j <= numpoints; j++) {
            phi1 = p1 + (((Real)(j - 10)) * 2.0 * Pi / ((Real)numpoints) / 10.0);
            orderedpair[1] = phi1;
            surface_fine.push_back(orderedpair);
        }
    }
    computeCoords(surface_fine,x_fine,y_fine,z_fine,r_fine);

    surface_comp.clear();
    for (register int i = 1; i <= numpoints; i++) {
        theta1 = t2 + (((Real)(i - 10) * Pi) / ((Real)numpoints) / 10.0);
        orderedpair[0] = theta1;
        for (register int j = 1; j <= numpoints; j++) {
            phi1 = p2 + (((Real)(j - 10)) * 2.0 * Pi / ((Real)numpoints) / 10.0);
            orderedpair[1] = phi1;
            surface_comp.push_back(orderedpair);
        }
    }
    computeCoords(surface_comp,x_comp,y_comp,z_comp,r_comp);

    tt1 = t1;
    tt2 = t2;
    pp1 = p1;
    pp2 = p2;

    if (verbose_) {
        cout << "Computing refined length... ";
        cout.flush();
    }
    try {
        rmax = computeRefinedlength(surface_fine,surface_comp,x_fine,y_fine,z_fine,
                                x_comp,y_comp,z_comp,r_fine,r_comp,lvec,
                                rm,tt1,tt2,pp1,pp2);
    }
    catch (tsquare::EOBException ex) { throw ex; }

    // cout << "Done!" << endl;
    // cout.flush();

    t1 = tt1;
    t2 = tt2;
    p1 = pp1;
    p2 = pp2;

    return rmax;
}

Real tsquare::Star::computeWidth(std::vector<std::vector<Real> > sc, std::vector<Real> xc,
                          std::vector<Real> yc, std::vector<Real> zc, std::vector<Real> lvec,
                          std::vector<Real> &wvec, Real &t1, Real &t2, Real &p1, Real &p2)
{
    using std::sqrt;
    using std::acos;
    using std::pow;
    using std::fabs;

    int nang = (int)(sqrt((Real)sc.size()));
    Real theta1,theta2,phi1,phi2,rm;
    Real x1,x2,y1,y2,z1,z2,wwx,wwy,wwz,wx,wy,wz,dot,ddot;
    Real tt1,tt2,pp1,pp2;
    Real lx,ly,lz;
    Real rmax = 0.0;

    std::vector<Real> orderedpair;
    std::vector<Real> x_fine, y_fine, z_fine, x_comp, y_comp, z_comp;
    std::vector<Real> r_fine, r_comp;
    std::vector<std::vector<Real> > surface_fine,surface_comp;
    
    
    lx = lvec[0];
    ly = lvec[1];
    lz = lvec[2];

    if ((sc.size() != xc.size()) || (sc.size() != yc.size())
         || (sc.size() != zc.size())) {
        throw tsquare::EOBException("Star","computeRefinedlength","coarse coordinates",
                           xc.size(),sc.size());
    }

    for (register int i = 0; i < sc.size(); i++) {
        theta1 = sc[i][0];
        phi1 = sc[i][1];
        x1 = xc[i];
        y1 = yc[i];
        z1 = zc[i];
        for (register int j = 0; j < sc.size(); j++) {
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

                if ((fabs(dot) < 0.1) && (rmax < rm)) {
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

    ddot = (wwx * lx) + (wwy * ly) + (wwz * lz);
    dot = acos(ddot);

    wvec.clear();
    wvec.resize(3,0.0);
    wvec[0] = wwx;
    wvec[1] = wwy;
    wvec[2] = wwz;

    // Now compute refined width

    surface_fine.clear();
    orderedpair.clear();
    orderedpair.resize(2,0.0);

    for (register int i = 1; i <= nang; i++) {
        theta1 = t1 + (((Real)(i - 10) * Pi) / ((Real)nang) / 10.0);
        orderedpair[0] = theta1;
        for (register int j = 1; j <= nang; j++) {
            phi1 = pp1 + (((Real)(j - 10)) * 2.0 * Pi / ((Real)nang) / 10.0);
            orderedpair[1] = phi1;
            surface_fine.push_back(orderedpair);
        }
    }
    computeCoords(surface_fine,x_fine,y_fine,z_fine,r_fine);

    surface_comp.clear();
    for (register int i = 1; i <= nang; i++) {
        theta1 = t2 + (((Real)(i - 10) * Pi) / ((Real)nang) / 10.0);
        orderedpair[0] = theta1;
        for (register int j = 1; j <= nang; j++) {
            phi1 = pp2 + (((Real)(j - 10)) * 2.0 * Pi / ((Real)nang) / 10.0);
            orderedpair[1] = phi1;
            surface_comp.push_back(orderedpair);
        }
    }
    computeCoords(surface_comp,x_comp,y_comp,z_comp,r_comp);

    tt1 = t1;
    tt2 = t2;
    pp1 = p1;
    pp2 = p2;

    rm = rmax;
    try {
        rmax = computeRefinedwidth(surface_fine, surface_comp, x_fine, y_fine, z_fine,
                               x_comp, y_comp, z_comp, r_fine, r_comp, lvec, wvec,
                               rm, ddot, tt1, tt2, pp1, pp2);
    }
    catch (tsquare::EOBException ex) { throw ex; }

    // TODO:  What does the rest of this function do?  Seems useless.

    t1 = tt1;
    t2 = tt2;
    p1 = pp1;
    p2 = pp2;

    dot = (wvec[0] * lvec[0]) + (wvec[1] * lvec[1]) + (wvec[2] * lvec[2]);
    dot = acos(dot);

    return rmax;
}

Real tsquare::Star::computeThickness(std::vector<std::vector<Real> > sc, std::vector<Real> xc,
                          std::vector<Real> yc, std::vector<Real> zc, std::vector<Real> rc,
                          std::vector<Real> lvec, std::vector<Real> wvec, std::vector<Real> &tvec,
                          Real &t1, Real &t2, Real &pp1, Real &pp2)
{
    using std::sqrt;
    using std::pow;
    using std::fabs;
    using std::acos;
    // using std::cout;
    // using std::endl;

    int nang = (int)(sqrt((Real)sc.size()));
    Real theta1,theta2,phi1,phi2;
    Real x1,y1,z1,x2,y2,z2;
    Real lx,ly,lz,wx,wy,wz,tx,ty,tz,ttx,tty,ttz;
    Real ppp1,ppp2,tt1,tt2;
    Real dotl,dotw,ddotl,ddotw;
    Real rm,rmax = 0.0;
    Real max_thickness;

    std::vector<std::vector<Real> > surface_fine, surface_comp;
    std::vector<Real> x_fine,y_fine,z_fine,x_comp,y_comp,z_comp;
    std::vector<Real> r_fine,r_comp,orderedpair;

    lx = lvec[0];
    ly = lvec[1];
    lz = lvec[2];
    
    wx = wvec[0];
    wy = wvec[1];
    wz = wvec[2];

    if ((sc.size() != xc.size()) || (sc.size() != yc.size())
         || (sc.size() != zc.size())) {
        throw tsquare::EOBException("Star","computeThickness","coarse coordinates",
                           sc.size(),xc.size());
    }
    for (register int i = 0; i < sc.size(); i++) {
        theta1 = sc[i][0];
        phi1 = sc[i][1];
        x1 = xc[i];
        y1 = yc[i];
        z1 = zc[i];
        for (register int j = 0; j < sc.size(); j++) {
            theta2 = sc[j][0];
            phi2 = sc[j][1];
            x2 = xc[j];
            y2 = yc[j];
            z2 = zc[j];

            rm = sqrt(pow(x1-x2,2.0)+pow(y1-y2,2.0)+pow(z1-z2,2.0));

            // Compute unit vector along trial thickness

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

                if ((fabs(dotl) < 0.1) && (fabs(dotw) < 0.1) && (rmax < rm)) {
                    rmax = rm;
                    t1 = theta1;
                    pp1 = phi1;
                    t2 = theta2;
                    pp2 = phi2;
                    ttx = tx;
                    tty = ty;
                    ttz = tz;
                }
            }
        }
    }

    tvec.clear();
    tvec.resize(3,0.0);
    tvec[0] = ttx;
    tvec[1] = tty;
    tvec[2] = ttz;

    /*
    cout << "Unrefined thickness = " << rmax << endl;
    cout.flush();
    */

    dotl = (ttx * lx) + (tty * ly) + (ttz * lz);
    dotw = (wx * ttx) + (wy * tty) + (wz * ttz);
    ddotl = fabs(dotl);
    ddotw = fabs(dotw);
    dotl = acos(dotl);
    dotw = acos(dotw);

    // Set final unit vector and angles for thickness to
    // first iteration of thickness

    tt1 = t1;
    ppp1 = pp1;
    tt2 = t2;
    ppp2 = pp2;

    // First set up the refined points to do pairwise search

    surface_fine.clear();
    orderedpair.clear();
    orderedpair.resize(2,0.0);

    for (register int i = 1; i <= nang; i++) {
        theta1 = t1 + (((Real)(i - 10) * Pi) / ((Real)nang) / 10.0);
        orderedpair[0] = theta1;
        for (register int j = 1; j <= nang; j++) {
            phi1 = pp1 + (((Real)(j - 10)) * 2.0 * Pi / ((Real)nang) / 10.0);
            orderedpair[1] = phi1;
            surface_fine.push_back(orderedpair);
        }
    }
    computeCoords(surface_fine,x_fine,y_fine,z_fine,r_fine);

    surface_comp.clear();
    for (register int i = 1; i <= nang; i++) {
        theta1 = t2 + (((Real)(i - 10) * Pi) / ((Real)nang) / 10.0);
        orderedpair[0] = theta1;
        for (register int j = 1; j <= nang; j++) {
            phi1 = pp2 + (((Real)(j - 10)) * 2.0 * Pi / ((Real)nang) / 10.0);
            orderedpair[1] = phi1;
            surface_comp.push_back(orderedpair);
        }
    }
    computeCoords(surface_comp,x_comp,y_comp,z_comp,r_comp);

    try {
        max_thickness = computeRefinedthickness(surface_fine,surface_comp,x_fine,
                                                y_fine,z_fine,x_comp,y_comp,z_comp,
                                                r_fine,r_comp,lvec,wvec,tvec,rmax,
                                                ddotl,ddotw,tt1,tt2,ppp1,ppp2);
    }
    catch (tsquare::EOBException ex) { throw ex; }

    t1 = tt1;
    t2 = tt2;
    pp1 = ppp1;
    pp2 = ppp2;

    return max_thickness;
}

Real tsquare::Star::computeTriaxialwidth(const int numpoints, const Real rstep,
                                    Real &rotate_angle, Real &mpa, Real &ampa,
                                    Real &amaxwidth, Real &aamw)
{
    using std::cout;
    using std::endl;
    using std::sqrt;
    using std::pow;
    using std::cos;
    using std::sin;

    Real theta1,phi1,phi2;
    Real area,rm,rmax = 0.0;
    Real testd2,testd2max;
    Real r1,r2,x1,y1;
    std::vector<Real> tvals,pvals;
    Real alpha,beta,gamma;
    Real rotate_step_size; // in degrees
    Real max_width = 0.0;

    std::vector<Real> r_comp;

    mpa = 0.0; // maximum projected area value
    ampa = 0.0; // angle rotated from the beginning of rotations at maximum projected area
    amaxwidth = 0.0; // absolute maximum width found, corresponding to W
    aamw = 0.0;      // angle rotated from the beginning of rotations where W was found

    // set the quadrature points for the projected areas

    try {
        setSurface(numpoints,nmax_from_anmfile_,true);
    }
    catch (tsquare::EOBException ex) { throw ex; }

    // begin compute max projected area and width

    rotate_step_size = rstep * Deg2Rad;

    tvals.clear();
    pvals.clear();
    for (register int i = 1; i <= numpoints; i++) {
        tvals.push_back((0.5 * Pi * xg_[i]) + (0.5 * Pi));
        pvals.push_back((Pi * xg_[i]) + Pi);
    }
    alpha = Pi/2.0;
    gamma = -alpha;

    beta = rotate_step_size;
    mpa = amaxwidth = 0.0;

    for (rotate_angle = 0.0; rotate_angle <= Pi; rotate_angle += rotate_step_size) {

        // begin project points onto xy plane

        r_comp.clear();
        r_comp.resize(numpoints,0.0);
        for (register int i = 0; i < numpoints; i++) {
            phi1 = pvals[i];
            rmax = 0.0;
            for (register int j = 0; j < numpoints; j++) {
                theta1 = tvals[j];
                rm = getR(theta1,phi1);
                x1 = rm * sin(theta1) * cos(phi1);
                y1 = rm * sin(theta1) * sin(phi1);
                rm = sqrt(pow(x1,2.0) + pow(y1,2.0));
                if (rm > rmax) rmax = rm;
            }
            r_comp[i] = rmax;
        }
        // end project points onto xy plane

        // begin find projected area
        testd2max = 0.0;
        area = 0.0;
        for (register int i = 0; i < numpoints; i++) {
            r1 = r_comp[i];
            r2 = r_comp[numpoints-i-1];
            phi1 = pvals[i];
            phi2 = pvals[numpoints-i-1];
            area += wg_[i+1] * r1 * r1;
            if (i > 0 && i < (numpoints/2)) {
                testd2 = (r1 * sin(phi1)) - (r2 * sin(phi2));
                if (testd2 > testd2max) {
                    testd2max = testd2;
                }
            }
        }
        area *= 0.5 * Pi;
        // end find projected area

        // begin update max area and width and angle
        if (area > mpa) {
            mpa = area;
            max_width = testd2max;
            if (max_width > amaxwidth) {
                amaxwidth = max_width;
                aamw = rotate_angle;
            }
            ampa = rotate_angle;
        }

        if (verbose_) {
            cout << "Rotating particle about x-axis by " << beta * Rad2Deg
                 << " degrees... ";
            cout.flush();
        }
        try {
            doRotate(alpha,beta,0.0);
            doRotate(0.0,0.0,gamma);
            setSurface(numpoints,nmax_from_anmfile_,true);
        }
        catch (tsquare::EOBException ex) { throw ex; }
        // end rotate the particle
    }
    
   

    if (verbose_) {
        cout << "Found maximum projected area of " << mpa << endl;
        cout << "Triaxial width (D2) is " << max_width << endl;
        cout << "Found D2 at angle " << ampa * Rad2Deg << endl;
        cout << "Maximum width (W) is " << amaxwidth << endl;
        cout << "Found W at angle " << aamw * Rad2Deg << endl;
        cout << "Current angle is " << rotate_angle * Rad2Deg << " deg" << endl;
    }

    return max_width;
}

Real tsquare::Star::computeTriaxialthickness(const int numpoints, const Real init_angle,
                                        const Real d2_angle, const Real w_angle,
                                        Real &absmax_thickness)
{
    using std::cout;
    using std::endl;
    using std::string;

    Real alpha,beta,gamma,theta1,phi1;
    Real t1,t2,pp1,pp2;
    Real thickness,triax_thickness;
    std::vector<std::vector<Real> > surface_coarse;
    std::vector<Real> x_coarse, y_coarse, z_coarse, r_coarse;
    std::vector<Real> lvec, wvec, tvec;
    std::vector<Real> rot_angle, orderedpair;

    // We are going to do two rotations about the x-axis, stored in rot_angle vector;
    // We first rotate the particle to the angle, d2_angle, where the triaxial width was found
    // The current angle is init_angle
    
    rot_angle.clear();
    rot_angle.resize(2,0.0);

    alpha = -Pi/2.0;
    gamma = -alpha;
    beta = d2_angle - init_angle;
    rot_angle[0] = beta;
    
    beta = w_angle - d2_angle;
    rot_angle[1] = beta;

    // For both these operations, the length axis is aligned parallel to the x-axis,
    // and the width axis, whether D2 or W, will be aligned parallel to the
    // y-axis, so we can set their unit vectors once and for all

    lvec.clear();
    lvec.resize(3,0.0);
    lvec[0] = 1.0;
    lvec[1] = lvec[2] = 0.0;

    wvec.clear();
    wvec.resize(3,0.0);
    wvec[1] = 1.0;
    wvec[0] = wvec[2] = 0.0;

    orderedpair.clear();
    orderedpair.resize(2,0.0);

    for (register int anum = 0; anum < 2; anum++) {
        if (verbose_) {
            cout << "Rotating particle about x-axis by "
                 << rot_angle[anum] * Rad2Deg << " degrees... ";
            cout.flush();
        }
        beta = -rot_angle[anum];
        doRotate(alpha,beta,0.0);
        doRotate(0.0,0.0,gamma);

         // Write the rotated SH coefficients to a file

        if (anum == 0) {
            try {
                string outname = name_ + "_rotated";
                std::ofstream out(outname.c_str());
                if (!out) {
                    throw tsquare::FileException("Star","calculateSHcoeffs",outname,
                                    "Could not open for writing");
                }
                for (register int n = 0; n <= nmax_; n++) {
                    for (register int m = n; m >= -n; m--) {
                        out << std::setfill(' ');
                        out << std::setw(5) << n << std::setw(5) << m;
                        out << std::setw(20) << std::setprecision(10) << std::fixed << std::real(a_[RDIM][n][m]);
                        out << std::setw(20) << std::setprecision(10) << std::fixed << std::imag(a_[RDIM][n][m]) << endl;
                    }
                }
                out.close();
            }
            catch (tsquare::FileException fex) {
                fex.printException();
            }
        }

        // end rotate the particle

        // First confirm that this is the correct orientation
        // Only as a debugging step, elminate this code block later on from here...

        /*
        tvals.clear();
        pvals.clear();
        for (register int i = 1; i <= numpoints; i++) {
            tvals.push_back((0.5 * Pi * xg_[i]) + (0.5 * Pi));
            pvals.push_back((Pi * xg_[i]) + Pi);
        }

        r_coarse.clear();
        r_coarse.resize(numpoints,0.0);
        for (register int i = 0; i < numpoints; i++) {
            phi1 = pvals[i];
            rmax = 0.0;
            for (register int j = 0; j < numpoints; j++) {
                theta1 = tvals[j];
                rm = getR(theta1,phi1);
                x1 = rm * sin(theta1) * cos(phi1);
                y1 = rm * sin(theta1) * sin(phi1);
                rm = sqrt(pow(x1,2.0) + pow(y1,2.0));
                if (rm > rmax) rmax = rm;
            }
            r_coarse[i] = rmax;
        }

        testd2max = 0.0;
        for (register int i = 0; i < numpoints; i++) {
            r1 = r_coarse[i];
            r2 = r_coarse[numpoints-i-1];
            phi1 = pvals[i];
            phi2 = pvals[numpoints-i-1];
            if (i > 0 && i < (numpoints/2)) {
                testd2 = (r1 * sin(phi1)) - (r2 * sin(phi2));
                if (testd2 > testd2max) {
                    testd2max = testd2;
                }
            }
        }
        if (anum == 0) {
            cout << "I think I found a D2 = " << testd2max << endl;
        } else {
            cout << "I think I found a W = " << testd2max << endl;
        }
        */
        // ... to here

        surface_coarse.clear();
        x_coarse.clear();
        y_coarse.clear();
        z_coarse.clear();

        for (register int i = 1; i <= numpoints; i++) {
            theta1 = 0.001 + ((((Real)(i - 1)) * Pi) / ((Real)numpoints));
            orderedpair[0] = theta1;
            for (register int j = 1; j <= numpoints; j++) {
                phi1 = ((Real)(j - 1)) * 2.0 * Pi / ((Real)numpoints);
                orderedpair[1] = phi1;
                surface_coarse.push_back(orderedpair);
            }
        }

        computeCoords(surface_coarse,x_coarse,y_coarse,z_coarse,r_coarse);

        // We now compute the triaxial thickness just like we would do the
        // ordinary thickness, only in this special orientation where the
        // projected area is parallel to the xy-plane
 
        try {
            thickness = computeThickness(surface_coarse,x_coarse,
                                         y_coarse,z_coarse,r_coarse,
                                         lvec,wvec,tvec,t1,t2,pp1,pp2);
        }
        catch (tsquare::EOBException ex) { throw ex; }

        if (anum == 0) {
            triax_thickness = thickness;
            cout << "Found D3 = " << thickness << endl;
            cout.flush();
        } else {
            absmax_thickness = thickness;
            cout << "Found T = " << thickness << endl;
            cout.flush();
        }

    }

    return triax_thickness;
}

void tsquare::Star::computeDimensions(bool triaxialcalc=true)
{
    using std::cout;
    using std::endl;
    using std::atan;
    using std::sin;
    using std::cos;
    using std::acos;
    using std::atan;
    using std::sqrt;
    using std::fabs;

    int nang = 50;
    Real theta1,phi1;
    Real rmax = 0.0;
    Real r1,r2;
    Real t1,t2,pp1,pp2;
    Real xn,yn,zn,dot;
    Real t1deg,t2deg,pp1deg,pp2deg;
    std::vector<Real> L1coords,L2coords;
    std::vector<Real> tvals,pvals;
    std::vector<Real> length_vec,width_vec,thickness_vec;
    std::vector<Real> L1temp,L2temp;
    Real alpha,beta;
    Real rotate_step_size = 2.0; // in degrees
    Real rotate_angle;
    Real max_projected_area = 0.0;
    Real absmax_width = 0.0;
    Real absmax_thickness = 0.0;
    Real max_width = 0.0;
    Real angle_at_max_projected_area = 0.0;
    Real angle_at_absmax_width = 0.0;
    Real max_thickness;

    // First store the old surface points, if there are any.

    std::vector<std::vector<Real> > surface_coarse,surface_fine,surface_comp;
    std::vector<Real> orderedpair;
    std::vector<Real> r_coarse,r_fine,r_comp;
    std::vector<Real> x_coarse,y_coarse,z_coarse;
    std::vector<Real> x_fine,y_fine,z_fine;
    std::vector<Real> x_comp,y_comp,z_comp;

    // begin compute coarse coordinates

    orderedpair.clear();
    orderedpair.resize(2,0.0);
    surface_coarse.clear();
    if (nang > 1) {
        for (register int i = 1; i <= nang; i++) {
            theta1 = 0.001 + ((((Real)(i - 1)) * Pi) / ((Real)nang));
            orderedpair[0] = theta1;
            for (register int j = 1; j <= nang; j++) {
                phi1 = ((Real)(j - 1)) * 2.0 * Pi / ((Real)nang);
                orderedpair[1] = phi1;
                surface_coarse.push_back(orderedpair);
            }
        }
    } else {
        throw tsquare::DataException("Star","computeDimensions",
                          "Number of search points must be greater than 1");
    }

    computeCoords(surface_coarse,x_coarse,y_coarse,z_coarse,r_coarse);

    // compute the length axis
    cout << "Computing length axis..." << endl;
    cout.flush();
    try {
        rmax = computeLength(surface_coarse,x_coarse,y_coarse,z_coarse,
                                length_vec,t1,t2,pp1,pp2);
    }
    catch (tsquare::EOBException ex) { throw ex; }

    dim_[0] = rmax;
    triaxialdim_[0] = rmax;

    cout << "Length (L = D1) = " << rmax << endl;
    cout.flush();

    r1 = getR(t1,pp1);
    r2 = getR(t2,pp2);

    t1deg = t1 * Rad2Deg / Pi;
    t2deg = t2 * Rad2Deg / Pi;
    pp1deg = pp1 * Rad2Deg / Pi;
    pp2deg = pp2 * Rad2Deg / Pi;

    L1coords.clear();
    L2coords.clear();
    L1coords.push_back(r1*sin(t1)*cos(pp1));
    L1coords.push_back(r1*sin(t1)*sin(pp1));
    L1coords.push_back(r1*cos(t1));
    L2coords.push_back(r2*sin(t2)*cos(pp2));
    L2coords.push_back(r2*sin(t2)*sin(pp2));
    L2coords.push_back(r2*cos(t2));

    if (verbose_) {
        cout << "End points are 1. r(" << t1deg << "," << pp1deg
             << ") = (" << L1coords[0];
        cout << "," << L1coords[1] << "," << L1coords[2]
             << "), distance from center = " << r1 << endl;
        cout << "               2. r(" << t2deg << "," << pp2deg
             << ") = (" << L2coords[0];
        cout << "," << L2coords[1] << "," << L2coords[2]
             << "), distance from center = " << r2 << endl;
        cout.flush();
    }

    // end compute length axis
    
    if (triaxialcalc) {

        // rotate the length axis parallel to x-axis

        L1temp = L1coords;
        L2temp = L2coords;
        if (fabs(L2coords[0]-L1coords[0]) > 0.0) {
            alpha = atan((L2coords[1]-L1coords[1])/(L2coords[0]-L1coords[0]));
            if (alpha > (Pi/2.0)) alpha -= Pi;
        } else {
            alpha = (Pi/2.0);
        }

        L1temp[0] = L1coords[0]*cos(alpha) + L1coords[1]*sin(alpha);
        L1temp[1] = L1coords[1]*cos(alpha) - L1coords[0]*sin(alpha);
        L2temp[0] = L2coords[0]*cos(alpha) + L2coords[1]*sin(alpha);
        L2temp[1] = L2coords[1]*cos(alpha) - L2coords[0]*sin(alpha);

        if (fabs(L2temp[0]-L1temp[0]) > 0.0) {
            beta = -atan((L2temp[2]-L1temp[2])/(L2temp[0]-L1temp[0]));
            if (beta < (-Pi/2.0)) beta += Pi;
        } else {
            beta = (Pi/2.0);
        }

        if (verbose_) {
            cout << "Rotating particle by (" << alpha * Rad2Deg << ","
                 << beta * Rad2Deg << ", 0.0)... ";
            cout.flush();
        }
        doRotate(alpha,beta,0.0);

        // done rotating the length axis

        // check on the new length orientation for verification
        // begin recomputing the coarse coordinates

        orderedpair.clear();
        orderedpair.resize(2,0.0);
        surface_coarse.clear();
        for (register int i = 1; i <= nang; i++) {
            theta1 = 0.001 + ((((Real)(i - 1)) * Pi) / ((Real)nang));
            orderedpair[0] = theta1;
            for (register int j = 1; j <= nang; j++) {
                phi1 = ((Real)(j - 1)) * 2.0 * Pi / ((Real)nang);
                orderedpair[1] = phi1;
                surface_coarse.push_back(orderedpair);
            }
        }

        computeCoords(surface_coarse,x_coarse,y_coarse,z_coarse,r_coarse);

        // compute the length axis
        try {
            rmax = computeLength(surface_coarse,x_coarse,y_coarse,z_coarse,
                                length_vec,t1,t2,pp1,pp2);
        }
        catch (tsquare::EOBException ex) { throw ex; }

        dim_[0] = rmax;
        triaxialdim_[0] = rmax;

        if (verbose_) {
            cout << "Length = " << dim_[0] << " at (" << t1 << ","
                 << pp1 << ") and (" << t2 << "," << pp2 << ")" << endl;
            cout.flush();
        }

        r1 = getR(t1,pp1);
        r2 = getR(t2,pp2);

        t1deg = t1 * Rad2Deg / Pi;
        t2deg = t2 * Rad2Deg / Pi;
        pp1deg = pp1 * Rad2Deg / Pi;
        pp2deg = pp2 * Rad2Deg / Pi;

        L1coords.clear();
        L2coords.clear();
        L1coords.resize(3,0.0);
        L2coords.resize(3,0.0);
        L1coords[0] = r1*sin(t1)*cos(pp1);
        L1coords[1] = r1*sin(t1)*sin(pp1);
        L1coords[2] = r1*cos(t1);
        L2coords[0] = r2*sin(t2)*cos(pp2);
        L2coords[1] = r2*sin(t2)*sin(pp2);
        L2coords[2] = r2*cos(t2);

        // end compute length axis

        zn = L1coords[2] - L2coords[2];
        yn = L1coords[1] - L2coords[1];
        xn = L1coords[0] - L2coords[0];
        r2 = sqrt((zn*zn) + (yn*yn) + (xn*xn));
        if (r2 > 0.0) {
            zn /= r2;
            yn /= r2;
            xn /= r2;
            cout << "Unit length vector = (" << xn << "," << yn << ","
                 << zn << ")" << endl;
            cout.flush();
            r1 = acos(fabs(xn)) * Rad2Deg;
            cout << "Note: Length and x-axis diverge by " << r1
                 << " degrees" << endl;
            cout.flush();
        } else {
            throw tsquare::DataException("Star","computeDimensions",
                          "Unit vector has zero length");
        }
        // end check length orientation
    }

    // begin to compute the width
    // we do this whether or not we are also computing triaxial dimensions
    // If triaxial dimensions also will be calculated, we use the orientation
    // of the width axis here to refine our search of projected areas, and
    // then we find an even more accurate width at the same time

    try {
        rmax = computeWidth(surface_coarse,x_coarse,y_coarse,z_coarse,
                                length_vec,width_vec,t1,t2,pp1,pp2);
    }
    catch (tsquare::EOBException ex) { throw ex; }

    dim_[1] = rmax;
    cout << "Width (W) = " << dim_[1] << endl;
    cout.flush();

    r1 = getR(t1,pp1);
    r2 = getR(t2,pp2);
    t1deg = t1 * Rad2Deg;
    t2deg = t2 * Rad2Deg;
    pp1deg = pp1 * Rad2Deg;
    pp2deg = pp2 * Rad2Deg;

    dot = (width_vec[0] * length_vec[0]) + (width_vec[1] * length_vec[1])
          + (width_vec[2] * length_vec[2]);
    dot = acos(dot);

    // end of computing width

    if (triaxialcalc) {

        // begin rotate width parallel to xy plane

        if (fabs(width_vec[1]) > 0.0) {
            beta = -atan(width_vec[2]/width_vec[1]);
            if (beta < (-Pi/2.0)) beta += Pi;
        } else {
            beta = (Pi/2.0);
        }

        if (verbose_) {
            cout << "Width axis makes an angle of " << beta * Rad2Deg
                 << " with xy plane" << endl;
            cout.flush();
        }

        // commenting this portion out for now, no rotation

        /*
        cout << "Rotating width axis 45 degrees beyone where width is parallel to xy plane" << endl;
     
        if (beta >= 0.0) {
            beta += (Pi/3.0);
        } else {
            beta -= (Pi/3.0);
        }

        cout << "Rotating width axis " << beta * Rad2Deg << " degrees" << endl;
        alpha = -Pi/2.0;
        gamma = -alpha;
        doRotate(alpha,beta,0.0);
        doRotate(0.0,0.0,gamma);
        */

        // end rotate width parallel to xy plane

        // begin compute triaxial width

        // Now that the particle is rotated into roughly the correct
        // position, we can scan the projected areas to find the
        // maximum projected area, the maximum width D2 of that
        // maximum projected area, and also W, which should be a refined,
        // more accurate value than the one already obtained

        // begin set quad points for projected area
        int numpoints = nang; // This needs to be even for some of
                             // the tricks used below

        rotate_step_size = 3.0;  // in degrees

        try {
            max_width = computeTriaxialwidth(numpoints,rotate_step_size,
                                             rotate_angle,max_projected_area,
                                             angle_at_max_projected_area,absmax_width,
                                             angle_at_absmax_width);
        }
        catch (tsquare::EOBException ex) { throw ex; }

        triaxialdim_[1] = max_width;
        dim_[1] = absmax_width;

        // end compute maximum projected area and D2 and W
        // end compute triaxial width

        width_vec[0] = width_vec[2] = 0.0;
        width_vec[1] = 1.0;
        length_vec[0] = 1.0;
        length_vec[1] = length_vec[2] = 0.0;

        // The vector max_thickness is a 2-element vector.  The first element holds
        // the triaxial thickness (perpendicular to both L and D2), and the second
        // element holds the absolute thickness (perpendicular to both L and W)

        try {
            max_thickness = computeTriaxialthickness(numpoints,rotate_angle,
                                                     angle_at_max_projected_area,
                                                     angle_at_absmax_width,
                                                     absmax_thickness);
        }
        catch (tsquare::EOBException ex) { throw ex; }

        triaxialdim_[2] = max_thickness;
        dim_[2] = absmax_thickness;

    } else {
        try {
            absmax_thickness = computeThickness(surface_coarse,x_coarse,
                                                y_coarse,z_coarse,r_coarse,
                                                length_vec,width_vec,
                                                thickness_vec,t1,t2,pp1,pp2);
        }
        catch (tsquare::EOBException ex) { throw ex; }
        dim_[2] = absmax_thickness;
        triaxialdim_[1] = dim_[1];
        triaxialdim_[2] = dim_[2];
    }

    ndim_[0] = dim_[0]/dim_[2];
    ndim_[1] = dim_[1]/dim_[2];
    ndim_[2] = 1.0;

    cout << endl;
    cout << "******************************************" << endl;
    cout << "Results: " << endl;
    cout << "      Length (L) = " << dim_[0] << endl;
    cout << "       Width (W) = " << dim_[1] << endl;
    cout << "   Thickness (T) = " << dim_[2] << endl;
    cout << endl;
    if (triaxialcalc) {
        cout << "              D1 = " << triaxialdim_[0] << endl;
        cout << "              D2 = " << triaxialdim_[1] << endl;
        cout << "              D3 = " << triaxialdim_[2] << endl;
    } else {
        cout << "  (Triaxial dimensions not calculated)" << endl;
    }
    cout << "******************************************" << endl;
    cout << endl;

    return;
}

void tsquare::Star::computeBoundingspheres(void)
{
    using std::cout;
    using std::endl;
    using std::sin;
    using std::cos;
    using std::sqrt;

    int count = 0;
    Real theta,phi,opposite_theta,opposite_phi,rmin,rmax,dcenter,rminit;
    Real stepsize,biggest_inscribed_so_far,smallest_enclosing_so_far;
    std::vector<Real> center;
    Real xdist,ydist,zdist,rtest,rtest_opposite;
    Real radius_here,is_r_thiscenter,es_r_thiscenter;
    Real st,sp,ct,cp;
    bool isdone_surfscan,esdone_surfscan,alldone_surfscan;
    bool isdone_radialscan,esdone_radialscan,alldone_radialscan;

    center.resize(3,0.0);

    rmin = 1.0e9;
    rmax = 0.0;
    cout << "  Initializing at center... ";
    for (register int i = 0; i < r_.size(); i++) {
        if (r_[i] < rmin) rmin = r_[i];
        if (r_[i] > rmax) rmax = r_[i];
    }
    biggest_inscribed_so_far = rmin;
    smallest_enclosing_so_far = rmax;
    stepsize = rmin / 100.0;
    rminit = rmin;
    if (verbose_) {
        cout << "rmin = " << rmin << ", rmax = " << rmax
             << ", stepsize = " << stepsize << endl;
        cout.flush();
    }
     
    tsquare::Sphere is(biggest_inscribed_so_far);
    tsquare::Sphere es(smallest_enclosing_so_far);

    try {
        for (register int i = 1; i <= ntheta_; i++) {
            for (register int j = 1; j <= nphi_; j++) {
                theta = surface_.at(count)[0];
                opposite_theta = Pi - theta;
                phi = surface_.at(count)[1];
                opposite_phi = (phi < Pi) ? (phi + Pi) : (phi - Pi);
                
                rtest = getR(theta,phi);
                rtest_opposite = getR(opposite_theta,opposite_phi);
                st = sin(theta);
                sp = sin(phi);
                ct = cos(theta);
                cp = cos(phi);
                isdone_radialscan = esdone_radialscan = alldone_radialscan = false;
                for (dcenter = stepsize; (dcenter < rtest) && !alldone_radialscan; dcenter += stepsize) {
                    center[0] = dcenter * st * cp;
                    center[1] = dcenter * st * sp;
                    center[2] = dcenter * ct;
                    is_r_thiscenter = 1.0e9;
                    es_r_thiscenter = 0.0;
                    if (rtest - dcenter < biggest_inscribed_so_far) isdone_radialscan = true;
                    if (rtest_opposite + dcenter > smallest_enclosing_so_far) esdone_radialscan = true;
                    alldone_radialscan = (isdone_radialscan && esdone_radialscan);

                    // This loop checks for updates to both the max inscribed sphere and the
                    // minimum enclosing sphere simultaneously.  Description below is just for
                    // the max inscribed sphere, but the description for the minimum enclosing
                    // sphere is analogous.
                    //
                    // The current largest circle radius, rmin, is checked against every
                    // surface point.  If every point is more than rmin away from the new circle
                    // center, then that becomes the updated max circle.  Otherwise, if a single
                    // point is found that is less than rmin from the new circle center, we finish
                    // and move the candidate circle center to the next position along this direction.

                    isdone_surfscan = isdone_radialscan;
                    esdone_surfscan = esdone_radialscan;
                    alldone_surfscan = alldone_radialscan;
                    for (register int ii = 1; ii < r_.size() && !alldone_surfscan; ii++) {
                        xdist = x_[ii] - center[0];
                        ydist = y_[ii] - center[1];
                        zdist = z_[ii] - center[2];
                        radius_here = sqrt((xdist*xdist) + (ydist*ydist) + (zdist*zdist));
                        if (!isdone_surfscan) {
                            if (radius_here < biggest_inscribed_so_far) {
                                // Already know that biggest inscribed sphere at this
                                // point is smaller than one I've already found, so don't
                                // check for inscribed spheres further along this direction
                                isdone_surfscan = true;
                            } else if (radius_here < is_r_thiscenter) {
                                is_r_thiscenter = radius_here;
                            }
                        }
                        if (!esdone_surfscan) {
                            if (radius_here > smallest_enclosing_so_far) {
                                esdone_surfscan = true;
                            } else if (radius_here > es_r_thiscenter) {
                                es_r_thiscenter = radius_here;
                            }
                        }
                        alldone_surfscan = (isdone_surfscan && esdone_surfscan);
    
                    }
                    if (!isdone_surfscan && (is_r_thiscenter > biggest_inscribed_so_far)) {
                        biggest_inscribed_so_far = is_r_thiscenter;
                        is.setRadius(biggest_inscribed_so_far);
                        is.setPosition(center[0],center[1],center[2]);
                    }
                    if (!esdone_surfscan && (es_r_thiscenter < smallest_enclosing_so_far)) {
                        smallest_enclosing_so_far = es_r_thiscenter;
                        es.setRadius(smallest_enclosing_so_far);
                        es.setPosition(center[0],center[1],center[2]);
                    }
                }

                if (verbose_) {
                    cout << "(" << surface_[count][0] * 180.0/Pi << ","
                         << surface_[count][1] * 180.0/Pi
                         << "): Max radius = " << biggest_inscribed_so_far
                         << " compared to " << rminit << endl;
                    cout.flush();
                }
                count++;
            }
        }
    }

    catch (std::out_of_range &oor) {
         throw tsquare::EOBException("Star","computeBoundingspheres","surface_",
                            surface_.size(),count);
    }

    inscribedsphere_ = is;
    enclosingsphere_ = es;
    cout << "Maximum inscribed sphere at (" << inscribedsphere_.getXr() << ","
         << inscribedsphere_.getYr() << "," << inscribedsphere_.getZr()
         << ") with radius " << inscribedsphere_.getRadius() << endl;
    cout << endl;
    cout << "Minimum enclosing sphere at (" << enclosingsphere_.getXr() << ","
         << enclosingsphere_.getYr() << "," << enclosingsphere_.getZr()
         << ") with radius " << enclosingsphere_.getRadius() << endl;
    cout << endl;
    return;
}

Real tsquare::Star::getIntegratedH(void)
{
    using std::sin;
    using std::sqrt;

    Real theta,phi,st,rp,rt,r;
    Real h = 0.0;
    Real area = 0.0;
    Real ddd,ddd2,f1;
    
    int count = 0;
    try {
        for (register int i = 1; i <= ntheta_; i++) {
            for (register int j = 1; j <= nphi_; j++) {
                theta = surface_.at(count)[0];
                phi = surface_.at(count)[1];
                st = sin(theta);
                r = getR(theta,phi);
                rp = getRp(theta,phi);
                rt = getRt(theta,phi);
                ddd = (rp * rp) + (rt * rt * st * st) + (r * r * st * st);
                ddd2 = sqrt(ddd);
                f1 = r * ddd2;
                area += (wg_.at(i) * wg_.at(j) * f1);
                h += (wg_.at(i) * wg_.at(j) * f1 * getH(theta,phi));
                count++;
            }
        }

        if (area > 0.0) {   
            h /= area;
        } else {
            throw tsquare::FloatException("Star","getIntegratedH","Divide by area = 0");
        }
        return h;
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getIntegratedH","surface_",
                            surface_.size(),count);
        ex.printException();
        exit(1);
    }
    catch (tsquare::FloatException flex) {
        flex.printException();
        exit(1);
    }
}

Real tsquare::Star::getIntegratedK(void)
{
    using std::sin;
    using std::sqrt;

    Real theta,phi,st,rp,rt,r;
    Real k = 0.0;
    Real area = 0.0;
    Real ddd,ddd2,f1;
    
    int count = 0;
    try {
        for (register int i = 1; i <= ntheta_; i++) {
            for (register int j = 1; j <= nphi_; j++) {
                theta = surface_.at(count)[0];
                phi = surface_.at(count)[1];
                st = sin(theta);
                r = getR(theta,phi);
                rp = getRp(theta,phi);
                rt = getRt(theta,phi);
                ddd = (rp * rp) + (rt * rt * st * st) + (r * r * st * st);
                ddd2 = sqrt(ddd);
                f1 = r * ddd2;
                area += (wg_.at(i) * wg_.at(j) * f1);
                k += (wg_.at(i) * wg_.at(j) * f1 * getK(theta,phi));
                count++;
            }
        }

        if (area > 0.0) {
            k /= area;
        } else {
            throw tsquare::FloatException("Star","getIntegratedK","Divide by area = 0");
        }
        return k;
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getIntegratedH","surface_",
                            surface_.size(),count);
        ex.printException();
        exit(1);
    }
    catch (tsquare::FloatException flex) {
        flex.printException();
        exit(1);
    }
}

int tsquare::Star::getI(std::vector<Real> &itensor, std::vector<Real> &rg2tensor, Real &isphere, Real &rg2sphere)
{
    using std::cout;
    using std::endl;
    using std::sin;
    using std::cos;
    using std::pow;
    using std::sqrt;
    using std::fabs;

    itensor.clear();
    rg2tensor.clear();
    isphere = rg2sphere = 0.0;
    
    Real i11,i22,i33,i12,i23,i13;
    Real rg11,rg22,rg33,rg12,rg23,rg13;
    Real theta,phi,st,ct,sp,cp,r,vol,fact;
    Real factor = (0.5 * Pi * Pi);
    int count = 0;

    i11 = i22 = i33 = i12 = i23 = i13 = fact = vol = 0.0;
    rg11 = rg22 = rg33 = rg12 = rg23 = rg13 = 0.0;
    try {
        for (register int i = 1; i <= ntheta_; i++) {
            for (register int j = 1; j <= nphi_; j++) {
                theta = surface_.at(count)[0];
                phi = surface_.at(count)[1];
                st = sin(theta);
                ct = cos(theta);
                sp = sin(phi);
                cp = cos(phi);
                r = getR(theta,phi);
                fact = 0.2 * wg_.at(i) * wg_.at(j) * pow(r,5.0) * st;
                vol += (wg_.at(i) * wg_.at(j) * pow(r,3.0) * st / 3.0);
                i11 += (fact * (1.0 - pow(cp*st,2.0)));
                i22 += (fact * (1.0 - pow(sp*st,2.0)));
                i33 += (fact * (1.0 - pow(ct,2.0)));
                i12 -= (fact * cp * sp * st * st);
                i13 -= (fact * cp * ct * st);
                i23 -= (fact * sp * ct * st);
                rg11 += (fact * pow(cp*st,2.0));
                rg22 += (fact * pow(sp*st,2.0));
                rg33 += (fact * pow(ct,2.0));
                rg12 += (fact * cp * sp * st * st);
                rg13 += (fact * cp * ct * st);
                rg23 += (fact * sp * ct * st);
                count++;
            }
        }

        vol *= factor;
        i11 *= factor;
        i22 *= factor;
        i33 *= factor;
        i12 *= factor;
        i13 *= factor;
        i23 *= factor;
        rg11 *= factor;
        rg22 *= factor;
        rg33 *= factor;
        rg12 *= factor;
        rg13 *= factor;
        rg23 *= factor;
        if (vol > 0.0) {
            i11 /= vol;
            i22 /= vol;
            i33 /= vol;
            i12 /= vol;
            i13 /= vol;
            i23 /= vol;
            itensor.push_back(i11);
            itensor.push_back(i22);
            itensor.push_back(i33);
            itensor.push_back(i12);
            itensor.push_back(i13);
            itensor.push_back(i23);
            rg11 /= vol;
            rg22 /= vol;
            rg33 /= vol;
            rg12 /= vol;
            rg13 /= vol;
            rg23 /= vol;
            rg2tensor.push_back(rg11);
            rg2tensor.push_back(rg22);
            rg2tensor.push_back(rg33);
            rg2tensor.push_back(rg12);
            rg2tensor.push_back(rg13);
            rg2tensor.push_back(rg23);
        } else {
            throw tsquare::FloatException("Star","getI","Divide by vol = 0");
        }

        // Find the radius of the volume equivalent sphere
        Real requiv = pow((3.0 * vol / 4.0 / Pi),(1.0/3.0));

        // Moment of inertia of this sphere
        isphere = 1.2 * requiv * requiv; // sum over all three principle axes

        // Square of radius of gyration of this sphere
        rg2sphere = 0.6 * requiv * requiv;

        // Now diagonalize the tensor, finding eigenvalues and eigenvectors
        // A is the original matrix, Q is a storage array, w stores the eigenvalues
        Real Q[3][3],A[3][3],w[3];
        Real rg2Q[3][3],rg2A[3][3],rg2w[3];
        const int n = 3;  // Dimension of the matrix

        Real sd,rg2sd,so;       // Sums of diagonal and off-diagonal elements
        Real tt;          // Tangent of rotation angle
        Real g,h,z,theta; // Temporary storage variables
        Real thresh;

        // Initialize A and rgA
        A[0][0] = itensor[0];
        A[1][1] = itensor[1];
        A[2][2] = itensor[2];
        A[0][1] = itensor[3];
        A[1][0] = itensor[3];
        A[0][2] = itensor[4];
        A[2][0] = itensor[4];
        A[1][2] = itensor[5];
        A[2][1] = itensor[5];
        rg2A[0][0] = rg2tensor[0];
        rg2A[1][1] = rg2tensor[1];
        rg2A[2][2] = rg2tensor[2];
        rg2A[0][1] = rg2tensor[3];
        rg2A[1][0] = rg2tensor[3];
        rg2A[0][2] = rg2tensor[4];
        rg2A[2][0] = rg2tensor[4];
        rg2A[1][2] = rg2tensor[5];
        rg2A[2][1] = rg2tensor[5];

        // Initialize Q and rgQ to the identity matrix
        for (int i = 0; i < n; i++) {
            Q[i][i] = 1.0;
            rg2Q[i][i] = 1.0;
            for (int j = 0; j < i; j++) {
                Q[i][j] = Q[j][i] = 0.0;
                rg2Q[i][j] = rg2Q[j][i] = 0.0;
            }
        }

        // Initialize w to diag(A)
        for (int i = 0; i < n; i++) {
            w[i] = A[i][i];
            rg2w[i] = rg2A[i][i];
        }
    
        // Calculate square of Tr(A) and of Tr(rgA);
        sd = rg2sd = 0.0;
        for (int i = 0; i < n; i++) {
            sd += fabs(w[i]);
            rg2sd += fabs(rg2w[i]);
        }
        sd = (sd * sd);
        rg2sd = (rg2sd * rg2sd);

        // Main loop for diagonalizing A
        so = 1.0e6;
        for (int nIter = 0; nIter < 50 && so != 0.0; nIter++) {
            // Test for convergence
            so = 0.0;
            for (int p = 0; p < n; p++) {
                for (int q = p+1; q < n; q++) {
                    so += fabs(A[p][q]);
                }
            }
            if (so != 0.0) {
                if (nIter < 4) {
                    thresh = 0.2 * so / ((Real)(n * n));
                } else {
                    thresh = 0.0;
                }
    
                // Do sweep
                for (int p = 0; p < n; p++) {
                    for (int q = p + 1; q < n; q++) {
                        g = 100.0 * fabs(A[p][q]);
                        if (nIter > 4 && fabs(w[p]) + g == fabs(w[p])
                            && fabs(w[q]) + g == fabs(w[q])) {
                            A[p][q] = 0.0;
                        } else if (fabs(A[p][q]) > thresh) {
                            //Calculate Jacobi transformation
                            h = w[q] - w[p];
                            if (fabs(h) + g == fabs(h)) {
                                tt = A[p][q] / h;
                            } else {
                                theta = 0.5 * h / A[p][q];
                                if (theta < 0.0) {
                                    tt = -1.0 / (sqrt(1.0 + (theta*theta)) - theta);
                                } else {
                                    tt = 1.0 / (sqrt(1.0 + (theta*theta)) + theta);
                                }
                            }
                            ct = 1.0 / sqrt(1.0 + (tt*tt));
                            st = tt * ct;
                            z = tt * A[p][q];
    
                            // Apply Jacobi transformation
                            A[p][q] = 0.0;
                            w[p] -= z;
                            w[q] += z;
                            for (int r = 0; r < p; r++) {
                                tt = A[r][p];
                                A[r][p] = (ct * tt) - (st * A[r][q]);
                                A[r][q] = (st * tt) + (ct * A[r][q]);
                            }
                            for (int r = p + 1; r < q; r++) {
                                tt = A[p][r];
                                A[p][r] = (ct * tt) - (st * A[r][q]);
                                A[r][q] = (st * tt) + (ct * A[r][q]);
                            }
                            for (int r = q + 1; r < n; r++) {
                                tt = A[p][r];
                                A[p][r] = (ct * tt) - (st * A[q][r]);
                                A[q][r] = (st * tt) + (ct * A[q][r]);
                            }
    
                            // Update eigenvectors
                            for (int r = 0; r < n; r++) {
                                tt = Q[r][p];
                                Q[r][p] = (ct * tt) - (st * Q[r][q]);
                                Q[r][q] = (st * tt) + (ct * Q[r][q]);
                            }
                        }
                    }
                }
            }
        }
    
        // Main loop for diagonalizing rg2A
        so = 1.0e6;
        for (int nIter = 0; nIter < 50 && so != 0.0; nIter++) {
            // Test for convergence
            so = 0.0;
            for (int p = 0; p < n; p++) {
                for (int q = p+1; q < n; q++) {
                    so += fabs(rg2A[p][q]);
                }
            }
            if (so != 0.0) {
                if (nIter < 4) {
                    thresh = 0.2 * so / ((Real)(n * n));
                } else {
                    thresh = 0.0;
                }
    
                // Do sweep
                for (int p = 0; p < n; p++) {
                    for (int q = p + 1; q < n; q++) {
                        g = 100.0 * fabs(rg2A[p][q]);
                        if (nIter > 4 && fabs(rg2w[p]) + g == fabs(rg2w[p])
                            && fabs(rg2w[q]) + g == fabs(rg2w[q])) {
                            rg2A[p][q] = 0.0;
                        } else if (fabs(rg2A[p][q]) > thresh) {
                            //Calculate Jacobi transformation
                            h = rg2w[q] - rg2w[p];
                            if (fabs(h) + g == fabs(h)) {
                                tt = rg2A[p][q] / h;
                            } else {
                                theta = 0.5 * h / rg2A[p][q];
                                if (theta < 0.0) {
                                    tt = -1.0 / (sqrt(1.0 + (theta*theta)) - theta);
                                } else {
                                    tt = 1.0 / (sqrt(1.0 + (theta*theta)) + theta);
                                }
                            }
                            ct = 1.0 / sqrt(1.0 + (tt*tt));
                            st = tt * ct;
                            z = tt * rg2A[p][q];
    
                            // Apply Jacobi transformation
                            rg2A[p][q] = 0.0;
                            rg2w[p] -= z;
                            rg2w[q] += z;
                            for (int r = 0; r < p; r++) {
                                tt = rg2A[r][p];
                                rg2A[r][p] = (ct * tt) - (st * rg2A[r][q]);
                                rg2A[r][q] = (st * tt) + (ct * rg2A[r][q]);
                            }
                            for (int r = p + 1; r < q; r++) {
                                tt = rg2A[p][r];
                                rg2A[p][r] = (ct * tt) - (st * rg2A[r][q]);
                                rg2A[r][q] = (st * tt) + (ct * rg2A[r][q]);
                            }
                            for (int r = q + 1; r < n; r++) {
                                tt = rg2A[p][r];
                                rg2A[p][r] = (ct * tt) - (st * rg2A[q][r]);
                                rg2A[q][r] = (st * tt) + (ct * rg2A[q][r]);
                            }
    
                            // Update eigenvectors
                            for (int r = 0; r < n; r++) {
                                tt = rg2Q[r][p];
                                rg2Q[r][p] = (ct * tt) - (st * rg2Q[r][q]);
                                rg2Q[r][q] = (st * tt) + (ct * rg2Q[r][q]);
                            }
                        }
                    }
                }
            }
        }

        // Finished diagonalizing.  Eigenvalues are in w vector, eigenvectors are in Q
        if (so != 0.0) {
            //Failed to converge; print error message of some sort
            cout << "Failed to diagonalize matrix, report non-diagonalized form." << endl;
            cout.flush();
        }

        itensor[0] = w[0];
        itensor[1] = w[1];
        itensor[2] = w[2];
        itensor[3] = itensor[4] = itensor[5] = 0.0;
    
        rg2tensor[0] = rg2w[0];
        rg2tensor[1] = rg2w[1];
        rg2tensor[2] = rg2w[2];
        rg2tensor[3] = rg2tensor[4] = rg2tensor[5] = 0.0;

        // Output eigenvector information
        cout << endl;
        cout << "Eigenvectors for I:" << endl;
        cout << " | " << Q[0][0] << " " << Q[1][0] << " " << Q[2][0] << " |" << endl;
        cout << " | " << Q[0][1] << " " << Q[1][1] << " " << Q[2][1] << " |" << endl;
        cout << " | " << Q[0][2] << " " << Q[1][2] << " " << Q[2][2] << " |" << endl;
        cout << endl;  
        cout << "Eigenvectors for S:" << endl;
        cout << " | " << rg2Q[0][0] << " " << rg2Q[1][0] << " " << rg2Q[2][0] << " |" << endl;
        cout << " | " << rg2Q[0][1] << " " << rg2Q[1][1] << " " << rg2Q[2][1] << " |" << endl;
        cout << " | " << rg2Q[0][2] << " " << rg2Q[1][2] << " " << rg2Q[2][2] << " |" << endl;
        cout << endl;  

        return 0;
    }
    catch (std::out_of_range &oor) {
        tsquare::EOBException ex("Star","getI","surface_",
                            surface_.size(),count);
        ex.printException();
        exit(1);
    }
    catch (tsquare::FloatException flex) {
        flex.printException();
        exit(1);
    }
}

int tsquare::Star::createVRML(const std::string vrmlname, bool render_inscribed_sphere,
             bool render_enclosing_sphere, bool render_hull,
             int surface_color_scheme, Real rotation)
{
    using std::cout;
    using std::endl;
    using std::cos;
    using std::sin;

    try {
        
      Real red,green,blue;
      Real zmax,zmin;
      Real xx,yy,zz,theta,phi;
      std::vector<Real> xvec;
      register int i,j;
      int num, number;
      Facet_iterator fi;
      Halfedge_around_facet_circulator hc;
      Real HULLRED = 0.923106;
      Real HULLGREEN = 0.923106;
      Real HULLBLUE = 0.923106;
      Real HULLTRANSPARENCY = 0.1;
        
      std::string filename = vrmlname + ".wrl";
      std::ofstream out(filename.c_str());
      if (!out) {
          throw tsquare::FileException("Star","createVRML",filename,"writing");
      }
        
      out << "#VRML V2.0 utf8" << endl;
      out << "Background {" << endl;
      out << "    skyColor [ 1 1 1 ]" << endl;
      out << "}" << endl;
      out << "NavigationInfo {" << endl;
      out << "\ttype [\"EXAMINE\",\"WALK\",\"FLY\"]" << endl;
      out << "}" << endl;
      out.flush();
      out << "Transform {" << endl;
      out << "  rotation 0 1 0 " << rotation * Pi/180.0 << endl;
      out << "  translation 0.0 0.0 -20.0" << endl;
      out << "  children [" << endl;
      out << "    Shape {" << endl;
      out << "      geometry IndexedFaceSet {" << endl;
      out << "        solid TRUE" << endl;
      out << "        ccw FALSE" << endl;
      out << "        coord Coordinate{" << endl;
      out << "          point [" << endl;
                      
      zmax = 0.0;
      zmin = 1.0e9;
      for (i = 0; i < surface_.size(); i++) {
          theta = surface_[i][0];
          phi = surface_[i][1];
          zz = r_[i] * cos(theta);
          if (zmax < zz) zmax = zz;
          if (zmin > zz) zmin = zz;
          if (i == 0) {
              out << "                    0.0000000000         0.0000000000 "
                  << zz << endl;
              out.flush();
          }
          xx = r_[i]*sin(theta)*cos(phi);
          yy = r_[i]*sin(theta)*sin(phi);
          out << "            " << xx << "        " << yy << "        "
              << zz << endl;
          out.flush();
          if (i == surface_.size() - 1) {
              out << "                    0.0000000000         0.0000000000 "
                  << zz << endl;
              out.flush();
          }
      }

      out.flush();
      out << "          ]" << endl;
      out << "        }" << endl;
      out << "        coordIndex [" << endl;
                  
      // To fill in the polygons making up the surface, we begin with the top
      // end cap triangles (i.e., all those that have a vertex in common with the
      // north pole).

      for (j = 1; j <= nphi_ - 1; j++) {
          out << "          0 " << j+1 << " " << j << " -1" << endl;
          out.flush();
      }
      out << "          0 1 " << ntheta_ << " -1" << endl;
      out.flush();

      for (i = 1; i <= ntheta_ - 1; i++) {
          for (j = 1; j <= nphi_; j++) {
              number = (ntheta_ * (i - 1)) + j;
              if (j == nphi_) {
                  num = number + 1 - nphi_;
                  out << "          " << number << " " << num << " ";
                  out << num+ntheta_ << " " << number+ntheta_ << " -1" << endl;
                  out.flush();
              } else {
                  out << "          " << number << " " << number+1 << " ";
                  out << number+1+ntheta_ << " " << number+ntheta_ << " -1" << endl;
                  out.flush();
              }
          }
      }

      number = ntheta_ * nphi_ + 1;
      for (j = 1; j <= nphi_ -1; j++) {
          num = ((ntheta_ - 1) * nphi_) + j;
          out << "          " << number << " " << num << " " << num+1 << " -1" << endl;
          out.flush();
      }
      out << "          " << number << " " << number-1 << " " << (ntheta_ - 1)*nphi_+1
          << " -1" << endl;
      out.flush();

      out << "        ]" << endl;

      // Code for color coding based on local roundness

      out << "        color Color {" << endl;
      out << "          color [" << endl;

      red = green = blue = 0.8;
      for (i = 0; i < surface_.size(); i++) {
          theta = surface_[i][0];
          phi = surface_[i][1];
          xx = r_[i]*sin(theta)*cos(phi);
          yy = r_[i]*sin(theta)*sin(phi);
          zz = r_[i] * cos(theta);
          if (surface_color_scheme != DOTPRODUCT && surface_color_scheme != WADELL) {
              red = blue = green = 0.8;  // uniform grey color
          } else {
              blue = 0.8 * roundness_[surface_color_scheme][i];
              red = 1.0 - roundness_[surface_color_scheme][i];
              if (red > 1.0) red = 1.0;
              green = red;
          }
          if (i == 0) out << " " << red << " " << green << " " << blue << endl;
          out << " " << red << " " << green << " " << blue << endl;
          if (i == (surface_.size() - 1)) out << " " << red << " " << green << " " << blue << endl;
      }
      out.flush();

      out << "          ]";
      out << "        }" << endl;
      out << "        colorPerVertex TRUE" << endl;
      out << "      }" << endl;
      out << "      appearance Appearance {" << endl;
      out << "        material Material {" << endl;
                  
      if (render_inscribed_sphere) {
        out << "          diffuseColor 0.923106 0.923106 0.923106" << endl;
        out << "          transparency 0.60000" << endl;
      } else {
        out << "          diffuseColor 0.8000 0.8000 0.8000" << endl;
      }
      out << "        }" << endl;
      out << "      }" << endl;
      out << "    }" << endl;
      out << "  ]" << endl;
      out << "}" << endl;

// Finished rendering the particle.  Work on the convex hull next.

      if (render_hull) {
    
        if (convexhull_.size_of_facets() == 0) {
            tsquare::DataException dex = tsquare::DataException("Star","createVRML",
                                 "Convex hull not yet created.  Doing it now.");
          dex.printException();
          getConvexity();
        }
        out << "Transform {" << endl;
        out << "  rotation 0 1 0 " <<  rotation * Pi/180.0 << endl;
        out << "  translation 0.0 0.0 -20.0" << endl;
        out << "  children [" << endl;
        out << "    Shape {" << endl;
        out << "      geometry IndexedFaceSet {" << endl;
        out << "        solid TRUE" << endl;
        out << "        ccw FALSE" << endl;
        out << "        coord Coordinate{" << endl;
        out << "          point [" << endl;
                      
        for (fi = convexhull_.facets_begin();
                  fi != convexhull_.facets_end(); ++fi) {
          hc = fi->facet_begin();
          do {
            xx = hc->vertex()->point().x();
            yy = hc->vertex()->point().y();
            zz = hc->vertex()->point().z();
            out << "                " << xx << " " << yy << " " << zz << endl;
            out.flush();
          } while ( ++hc != fi->facet_begin() );
        }
        out << "          ]" << endl;
        out << "        }" << endl;
        out << "        coordIndex [" << endl;
                  
        for (i = 0; i < convexhull_.size_of_facets(); i++) {
          num = 3 * i;
          out << "            " << num << " " << num + 1 << " "
              << num + 2 << " -1" << endl;
          out.flush();
        }

        out << "        ]" << endl;
        out << "      }" << endl;
        out << "      appearance Appearance {" << endl;
        out << "        material Material {" << endl;
                  
        out << "          diffuseColor " << HULLRED << " "
                                         << HULLGREEN << " "
                                         << HULLBLUE << endl;
        out << "          transparency " << HULLTRANSPARENCY << endl;
        out << "        }" << endl;
        out << "      }" << endl;
        out << "    }" << endl;
        out << "  ]" << endl;
        out << "}" << endl;
        out.flush();
      }

// Finished rendering the convex hull.  Work on the inscribed sphere next.

      if (render_inscribed_sphere) {
    
        if (inscribedsphere_.getRadius() <= 0.0) {
          cout << "Inscribed sphere not created yet.  Doing it now." << endl;
          computeBoundingspheres();
        }

        out << "Transform {" << endl;
        out << "  translation " << inscribedsphere_.getXr() << " " <<
        inscribedsphere_.getYr() << " " << inscribedsphere_.getZr() << endl;
        out << "    children [" << endl;
        out << "    Shape {" << endl;
        out << "      geometry Sphere { radius "
            << inscribedsphere_.getRadius() << " }" << endl;
        out << "      appearance Appearance {material Material "
            << "{diffuseColor 0 0 1}}" << endl;
        out << "    }" << endl;
        out << "    ]" << endl;
        out << "}" << endl;
        out.flush();
        out << "Transform {" << endl;
        out << "  translation " << inscribedsphere_.getXr() << " " <<
             inscribedsphere_.getYr() << " " << inscribedsphere_.getZr() << endl;
        out << "    children [" << endl;
        out << "    Shape {" << endl;
        out << "      geometry Sphere { radius "
            << inscribedsphere_.getRadius() << " }" << endl;
        out << "      appearance Appearance {material Material "
            << "{diffuseColor 0 0 1}}" << endl;
        out << "    }" << endl;
        out << "    ]" << endl;
        out << "}" << endl;
        out.flush();
      }

// Finished rendering the inscribed sphere.  Work on the enclosing sphere next.

      if (render_enclosing_sphere) {

        if (enclosingsphere_.getRadius() <= 0.0) {
          cout << "Enclosing sphere not created yet.  Doing it now." << endl;
          computeBoundingspheres();
        }
        out << "Transform {" << endl;
        out << "  translation " << enclosingsphere_.getXr() << " " <<
             enclosingsphere_.getYr() << " " << enclosingsphere_.getZr() << endl;
        out << "    children [" << endl;
        out << "    Shape {" << endl;
        out << "      geometry Sphere { radius "
            << enclosingsphere_.getRadius() << " }" << endl;
        out << "      appearance Appearance {material Material "
            << "{diffuseColor 1 0 0 transparency 0.6}}" << endl;
        out << "    }" << endl;
        out << "    ]" << endl;
        out << "}" << endl;
        out.flush();
        out << "Transform {" << endl;
        out << "  translation " << enclosingsphere_.getXr() << " "
            << enclosingsphere_.getYr() << " "
            << enclosingsphere_.getZr() << endl;
        out << "    children [" << endl;
        out << "    Shape {" << endl;
        out << "      geometry Sphere { radius "
            << enclosingsphere_.getRadius() << " }" << endl;
        out << "      appearance Appearance {material Material "
            << "{diffuseColor 1 0 0 transparency 0.6}}" << endl;
        out << "    }" << endl;
        out << "    ]" << endl;
        out << "}" << endl;
        out.flush();
      }
        
      out.close();
      return(0);
    }

    catch (tsquare::FileException fex) { throw fex; }
}

int tsquare::Star::createSHFile(const std::string shname)
{
    using std::endl;

    try {
        std::ofstream out(shname.c_str());
        if (!out) {
            throw tsquare::FileException("Star","createSHFile",shname,
                            "Could not open for writing");
        }
        for (register int n = 0; n <= nmax_; n++) {
            for (register int m = n; m >= -n; m--) {
                out << std::setfill(' ');
                out << std::setw(5) << n << std::setw(5) << m;
                out << std::setw(20) << std::setprecision(10) << std::fixed << std::real(a_[RDIM][n][m]);
                out << std::setw(20) << std::setprecision(10) << std::fixed << std::imag(a_[RDIM][n][m]) << endl;
            }
        }
        out.close();
        return(0);
    }
    catch (tsquare::FileException fex) {
        fex.printException();
    }
}

int tsquare::Star::createDigitizedImage(const std::string digitizedname, Real resolution)
{
    using std::cout;
    using std::endl;
    using std::cos;
    using std::acos;
    using std::sin;
    using std::atan2;
    using std::sqrt;

    try {
        std::vector<std::vector<std::vector<int> > > boundingBox;
        std::vector<int> oneDvector;
        std::vector<std::vector<int> > twoDvector;

        int partc = 0;
        Real rsurf,rtest,theta,phi;
        Real xt,yt,zt,xt2,yt2,zt2;

        std::string filename = digitizedname + ".txt";
        std::ofstream out(filename.c_str());
        if (!out) {
            throw tsquare::FileException("Star","createDigitizedImage",filename,"writing");
        }
        
        // First find box necessary box dimensions in all three Cartesian directions
        // Could probably find an optimum by making the length lie along the body diagonal
        // of the box.  Maybe will do that later.

        Real xmin,xmax,ymin,ymax,zmin,zmax;

        xmin = ymin = zmin = 0.0;
        xmax = ymax = zmax = 0.0;

        cout << "ntheta = " << ntheta_ << ", nphi = " << nphi_ << endl;
        cout.flush();
        for (unsigned int nt = 0; nt < ntheta_; nt++) {
            theta = 0.5 * Pi * (xg_[nt] + 1.0);
            for (unsigned int np = 0; np < nphi_; np++) {
                phi = Pi * (xg_[np] + 1.0);
                rsurf = getR(theta,phi);
                rsurf *= resolution;
                xt = rsurf * sin(theta) * cos(phi);
                yt = rsurf * sin(theta) * sin(phi);
                zt = rsurf * cos(theta);
                if (xt > xmax) xmax = xt;
                if (yt > ymax) ymax = yt;
                if (zt > zmax) zmax = zt;
                if (xt < xmin) xmin = xt;
                if (yt < ymin) ymin = yt;
                if (zt < zmin) zmin = zt;
            }
        }
    
        // Now size the 3D vector to fit the particle as is
        // It is better for each dimension to be an odd number so that there is
        // a well-defined body center voxel
    
        int xdim = (int)(xmax - xmin + 1);
        int ydim = (int)(ymax - ymin + 1);
        int zdim = (int)(zmax - zmin + 1);
        if (xdim%2 == 0) xdim += 1;
        if (ydim%2 == 0) ydim += 1;
        if (zdim%2 == 0) zdim += 1;
    
        // Find body center in the box
        int xc = xdim / 2;
        int yc = ydim / 2;
        int zc = zdim / 2;

        cout << endl;
        cout << "xdim = " << xdim << ", ydim = " << ydim << ", zdim = " << zdim << endl;
        cout << "xc = " << xc << ", yc = " << yc << ", zc = " << zc << endl;
        cout.flush();

        oneDvector.clear();
        oneDvector.resize(zdim,0);
        twoDvector.clear();
        twoDvector.resize(ydim,oneDvector);
        boundingBox.clear();
        boundingBox.resize(xdim,twoDvector);

        cout << endl << "boundingBox initialized:" << endl;
        cout << "boundingBox[" << boundingBox.size() << "]";
        cout.flush();
        cout << "[" << boundingBox[1].size() << "]";
        cout.flush();
        cout << "[" << boundingBox[1][1].size() << "]" << endl;
        cout.flush();

        int zupperbound = 0;
        int zlowerbound = boundingBox[1][1].size();
        for (int nz = 0; nz < zdim; nz++) {
            // Actual z-coordinate of body center of voxel in real space
            // Using box center as origin
            zt = (Real)(nz + 0.5) - (Real)(zc + 0.5);
            zt2 = zt * zt;
            for (int ny = 0; ny < ydim; ny++) {
                yt = (Real)(ny + 0.5) - (Real)(yc + 0.5);
                yt2 = yt * yt;
                for (int nx = 0; nx < xdim; nx++) {
                    xt = (Real)(nx + 0.5) - (Real)(xc + 0.5);
                    xt2 = xt * xt;
                    rtest = sqrt(xt2 + yt2 + zt2);
                    if (rtest > 0.0) {
                        theta = acos(zt/rtest);
                    } else {
                        theta = 0.0;
                    }
                    phi = atan2(yt,xt);
                    rsurf = getR(theta,phi);
                    rsurf *= resolution;
                    if (rtest <= rsurf) {
                        boundingBox[nx][ny][nz] = 1;
                        if (nz > zupperbound) zupperbound = nz;
                        if (nz < zlowerbound) zlowerbound = nz;
                    }
                }
            }
        }

        // Now the bounding box is populated, so we can output the image
        out << "Xsize: " << xdim << endl;
        out << "Ysize: " << ydim << endl;
        out << "Zsize: " << zupperbound - zlowerbound + 1 << endl;
        out << "Resolution: " << resolution;
        for (int nz = zlowerbound; nz <= zupperbound; nz++) {
            for (int ny = 0; ny < ydim; ny++) {
                for (int nx = 0; nx < xdim; nx++) {
                    out << endl << boundingBox[nx][ny][nz];
                }
            }
        }
        out.flush();
        out.close();

        return (0);
    }

    catch (tsquare::FileException fex) { throw fex; }
}

tsquare::Star& tsquare::Star::operator=(const tsquare::Star& that)
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
    x_ = that.getX();
    y_ = that.getY();
    z_ = that.getZ();

    return *this;
}
