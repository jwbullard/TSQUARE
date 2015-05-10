/**
@file
@brief Method definitions for Box class.
*/
#include "Box.h"

tsquare::Box::Box()
{
    halfdim_.clear();
    halfdim_.resize(3,1.0);
    volume_ = 8.0;
    area_ = 24.0;
    narea_ = 1.2407;
    ngc_ = diam_ = itrace_ = 0.0;
    dim_.clear();
    dim_.resize(3,2.0);
    ndim_ = dim_;
    triaxialdim_ = dim_;
    r_.clear();
    surface_.clear();
    xg_.clear();
    wg_.clear();
    hullvolume_ = hullarea_ = 0.0;
    verbose_ = false;

    inscribedsphere_.setXr(0.0);
    inscribedsphere_.setYr(0.0);
    inscribedsphere_.setZr(0.0);
    inscribedsphere_.setRadius(1.0);

    enclosingsphere_.setXr(0.0);
    enclosingsphere_.setYr(0.0);
    enclosingsphere_.setZr(0.0);
    enclosingsphere_.setRadius(std::sqrt(3.0));
}

tsquare::Box::Box(Real thickness, Real width, Real length, bool verbose) {
    verbose_ = verbose;

    using std::sqrt;

    if (thickness <= 0.0 || width <= 0.0 || length <= 0.0) {
        throw DataException("Box","Box",
                            "All dimensions must be positive numbers");
    }

    halfdim_.clear();
    halfdim_.resize(3);
    halfdim_[0] = 0.5 * length;
    halfdim_[1] = 0.5 * width;
    halfdim_[2] = 0.5 * thickness;
    std::sort(halfdim_.begin(),halfdim_.end());
    Real save = halfdim_[0];
    halfdim_[0] = halfdim_[2];
    halfdim_[2] = save;

    inscribedsphere_.setXr(0.0);
    inscribedsphere_.setYr(0.0);
    inscribedsphere_.setZr(0.0);
    inscribedsphere_.setRadius(halfdim_[0]);

    enclosingsphere_.setXr(0.0);
    enclosingsphere_.setYr(0.0);
    enclosingsphere_.setZr(0.0);

    Real rad2 = 0.0;
    for (int i = 0; i < 3; i++) rad2 += (halfdim_[i] * halfdim_[i]);
    enclosingsphere_.setRadius(sqrt(rad2));

    surface_.clear();
    r_.clear();
    xg_.clear();
    wg_.clear();
    dim_.clear();
    ndim_.clear();
    triaxialdim_.clear();
    dim_.resize(3,0.0);
    dim_[0] = 2.0 * halfdim_[0];
    dim_[1] = 2.0 * halfdim_[1];
    dim_[2] = 2.0 * halfdim_[2];
    ndim_.resize(3,0.0);
    ndim_[0] = dim_[0] / dim_[2];
    ndim_[1] = dim_[1] / dim_[2];
    ndim_[2] = 1.0;
    triaxialdim_ = dim_;
}

void tsquare::Box::setDimensions(Real thickness, Real width, Real length)
{
    using std::sqrt;

    if (thickness <= 0.0 || width <= 0.0 || length <= 0.0) {
        throw DataException("Box","setDimensions",
                            "All dimensions must be positive numbers");
    }

    halfdim_.clear();
    halfdim_.resize(3);
    halfdim_[0] = 0.5 * length;
    halfdim_[1] = 0.5 * width;
    halfdim_[2] = 0.5 * thickness;
    std::sort(halfdim_.begin(),halfdim_.end());
    Real save = halfdim_[0];
    halfdim_[0] = halfdim_[2];
    halfdim_[2] = save;

    halfdim_[0] *= 0.5;
    halfdim_[1] *= 0.5;
    halfdim_[2] *= 0.5;

    inscribedsphere_.setXr(0.0);
    inscribedsphere_.setYr(0.0);
    inscribedsphere_.setZr(0.0);
    inscribedsphere_.setRadius(halfdim_[0]);

    enclosingsphere_.setXr(0.0);
    enclosingsphere_.setYr(0.0);
    enclosingsphere_.setZr(0.0);

    Real rad2 = 0.0;
    for (int i = 0; i < 3; i++) rad2 += (halfdim_[i] * halfdim_[i]);
    enclosingsphere_.setRadius(sqrt(rad2));

    surface_.clear();
    r_.clear();
    xg_.clear();
    wg_.clear();
    dim_.clear();
    ndim_.clear();
    triaxialdim_.clear();
    dim_.resize(3,0.0);
    dim_[0] = 2.0 * halfdim_[0];
    dim_[1] = 2.0 * halfdim_[1];
    dim_[2] = 2.0 * halfdim_[2];
    ndim_.resize(3,0.0);
    ndim_[0] = dim_[0] / dim_[2];
    ndim_[1] = dim_[1] / dim_[2];
    ndim_[2] = 1.0;
    triaxialdim_ = dim_;
}

void tsquare::Box::setDimensions(std::vector<Real> dims)
{
    if (dims.size() < 3) {
        throw DataException("Box","setDimensions",
                            "All dimensions must be passed to this function");
    }
    
    setDimensions(dims[0],dims[1],dims[2]);
    return;
}

Real tsquare::Box::getR(const Real theta, const Real phi)
{
    using std::atan;
    using std::sin;
    using std::cos;

    Real pcrit = atan(halfdim_[1]/halfdim_[0]);
    Real tcrit,rad;
    Real t = theta;
    Real p = phi;
  
    if (p >= 1.5 * Pi) p = (2.0 * Pi) - p;
    if (p >= Pi) p -= Pi;
    if (p >= 0.5 * Pi) p = Pi - p;
    if (t > 0.5 * Pi) t = Pi - t;
  
    if (p <= pcrit) {
        tcrit = atan(halfdim_[0]/halfdim_[2]/cos(p));
        if (t <= tcrit) {
            rad = halfdim_[2] / cos(t);
        } else {
            rad = halfdim_[0] / sin(t) / cos(p);
        }
    } else {
        tcrit = atan(halfdim_[1]/halfdim_[2]/sin(p));
        if (t <= tcrit) {
            rad = halfdim_[2] / cos(t);
        } else {
            rad = halfdim_[1] / sin(t) / sin(p);
        }
    }
    return rad;

} 

Real tsquare::Box::getRp(Real theta, const Real phi)
{
    using std::atan;
    using std::sin;
    using std::cos;
    using std::fabs;

    Real pcrit = std::atan(halfdim_[1]/halfdim_[0]);
    Real tcrit;
    Real t = theta;
    Real p = phi;
    Real rp = 0;
  
    if (p >= 1.5 * Pi) p = (2.0 * Pi) - p;
    if (p >= Pi) p -= Pi;
    if (p >= 0.5 * Pi) p = Pi - p;
    if (t > 0.5 * Pi) t = Pi - t;
  
    if (p <= pcrit) {
        tcrit = atan(halfdim_[0]/halfdim_[2]/cos(p));
        if (t <= tcrit) {
            rp = 0.0;
        } else if (fabs(sin(t)) < 1.0e-9 || fabs(cos(p)) < 1.0e-9) {
            rp = 1.0e9;
            if (sin(t) < 0.0) rp = -1.0e9;
        } else {
            rp = halfdim_[0]  * sin(p) / sin(t) / cos(p) / cos(p);
        }
    } else {
        tcrit = atan(halfdim_[1]/halfdim_[2]/sin(p));
        if (t <= tcrit) {
            rp = 0.0;
        } else if (fabs(sin(t)) < 1.0e-9 || fabs(sin(p)) < 1.0e-9) {
            rp = 1.0e9;
            if (sin(t) < 0.0) rp = -1.0e9;
        } else {
            rp = -halfdim_[1]  * cos(p) / sin(t) / sin(p) / sin(p);
        }
    }
    return rp;

}

Real tsquare::Box::getRpp(Real theta, const Real phi)
{
    using std::atan;
    using std::sin;
    using std::cos;
    using std::fabs;

    Real pcrit = atan(halfdim_[1]/halfdim_[0]);
    Real tcrit,rad;
    Real t = theta;
    Real p = phi;
    Real rpp = 0.0;
  
    if (p >= 1.5 * Pi) p = (2.0 * Pi) - p;
    if (p >= Pi) p -= Pi;
    if (p >= 0.5 * Pi) p = Pi - p;
    if (t > 0.5 * Pi) t = Pi - t;
  
    if (p <= pcrit) {
        tcrit = atan(halfdim_[0]/halfdim_[2]/cos(p));
        if (t <= tcrit) {
            rpp = 0.0;
        } else if (fabs(sin(t)) < 1.0e-9 || fabs(cos(p)) < 1.0e-9) {
            rpp = 1.0e9;
            if (sin(t) * cos(p) < 0.0) rpp = -1.0e9;
        } else {
            rpp = (halfdim_[0] / sin(t) / cos(p)) * ( 1.0 + (2.0 * tan(p) * tan(p)) );
        }
    } else {
        tcrit = atan(halfdim_[1]/halfdim_[2]/sin(p));
        if (t <= tcrit) {
            rpp = 0.0;
        } else if (fabs(sin(t)) < 1.0e-9 || fabs(sin(p)) < 1.0e-9) {
            rpp = 1.0e9;
            if (sin(t) * sin(p) < 0.0) rpp = -1.0e9;
        } else {
            rpp = (halfdim_[1] / sin(t) / sin(p)) * ( 1.0 + (2.0 * cos(p) * cos(p) / sin(p) / sin(p)) );
        }
    }
    return rpp;

}

Real tsquare::Box::getRt(Real theta, const Real phi)
{
    using std::atan;
    using std::sin;
    using std::cos;
    using std::fabs;

    Real pcrit = atan(halfdim_[1]/halfdim_[0]);
    Real tcrit,rad;
    Real t = theta;
    Real p = phi;
    Real rt = 0.0;
  
    if (p >= 1.5 * Pi) p = (2.0 * Pi) - p;
    if (p >= Pi) p -= Pi;
    if (p >= 0.5 * Pi) p = Pi - p;
    if (t > 0.5 * Pi) t = Pi - t;
  
    if (p <= pcrit) {
        tcrit = atan(halfdim_[0]/halfdim_[2]/cos(p));
        if (t <= tcrit) {
            if (fabs(cos(t)) < 1.0e-9) {
                rt = 1.0e9 * sin(t);
            } else {
                rt = halfdim_[2]  * sin(t) / cos(t) / cos(t);
            }
        } else {
            if (fabs(cos(p)) < 1.0e-9 || fabs(sin(t)) < 1.0e-9) {
                rt = 1.0e9 * cos(t);
                if (cos(p) < 0.0) rt = -1.0e9 * cos(t);
            } else {
                rt = -halfdim_[0] * cos(t) / cos(p) / sin(t) / sin(t);
            }
        }
    } else {
        tcrit = atan(halfdim_[1]/halfdim_[2]/sin(p));
        if (t <= tcrit) {
            if (fabs(cos(t)) < 1.0e-9) {
                rt = 1.0e9 * sin(t);
            } else {
                rt = halfdim_[2]  * sin(t) / cos(t) / cos(t);
            }
        } else {
            if (fabs(sin(p)) < 1.0e-9 || fabs(sin(t)) < 1.0e-9) {
                rt = 1.0e9 * cos(t);
                if (sin(p) < 0.0) rt = -1.0e9 * cos(t);
            } else {
                rt = -halfdim_[1] * cos(t) / sin(p) / sin(t) / sin(t);
            }
        }
    }
    return rt;

}

Real tsquare::Box::getRtt(Real theta, const Real phi)
{
    using std::tan;
    using std::atan;
    using std::sin;
    using std::cos;
    using std::fabs;

    Real pcrit = atan(halfdim_[1]/halfdim_[0]);
    Real tcrit,rad;
    Real t = theta;
    Real p = phi;
    Real rtt = 0.0;
  
    if (p >= 1.5 * Pi) p = (2.0 * Pi) - p;
    if (p >= Pi) p -= Pi;
    if (p >= 0.5 * Pi) p = Pi - p;
    if (t > 0.5 * Pi) t = Pi - t;
  
    if (p <= pcrit) {
        tcrit = atan(halfdim_[0]/halfdim_[2]/cos(p));
        if (t <= tcrit) {
            if (fabs(cos(t)) < 1.0e-9) {
                rtt = 1.0e9;
                if (cos(t) < 0.0) rtt = -1.0e9;
            } else {
                rtt = (halfdim_[2] / cos(t)) * (1.0 + (2.0 * tan(t) * tan(t)));
            }
        } else {
            if (fabs(sin(t)) < 1.0e-9 || fabs(cos(p)) < 1.0e-9) {
                rtt = 1.0e9;
                if (sin(t) * cos(p) < 0.0) rtt = -1.0e9;
            } else {
                rtt = (halfdim_[0] / sin(t) / cos(p)) * (1.0 + (2.0 * cos(t) * cos(t) / sin(t) / sin(t)));
            }
        }
    } else {
        tcrit = atan(halfdim_[1]/halfdim_[2]/sin(p));
        if (t <= tcrit) {
            if (fabs(cos(t)) < 1.0e-9) {
                rtt = 1.0e9;
                if (cos(t) < 0.0) rtt = -1.0e9;
            } else {
                rtt = (halfdim_[2] / cos(t)) * (1.0 + (2.0 * tan(t) * tan(t)));
            }
        } else {
            if (fabs(sin(t)) < 1.0e-9 || fabs(sin(p)) < 1.0e-9) {
                rtt = 1.0e9;
                if (sin(t) * sin(p) < 0.0) rtt = -1.0e9;
            } else {
                rtt = (halfdim_[1] / sin(t) / sin(p)) * (1.0 + (2.0 * cos(t) * cos(t) / sin(t) / sin(t)));
            }
        }
    }
    return rtt;

}

Real tsquare::Box::getRtp(Real theta, const Real phi)
{
    using std::atan;
    using std::sin;
    using std::cos;
    using std::fabs;

    Real pcrit = atan(halfdim_[1]/halfdim_[0]);
    Real tcrit,rad;
    Real t = theta;
    Real p = phi;
    Real rtp = 0.0;
  
    if (p >= 1.5 * Pi) p = (2.0 * Pi) - p;
    if (p >= Pi) p -= Pi;
    if (p >= 0.5 * Pi) p = Pi - p;
    if (t > 0.5 * Pi) t = Pi - t;
  
    if (p <= pcrit) {
        tcrit = atan(halfdim_[0]/halfdim_[2]/cos(p));
        if (t <= tcrit) {
            rtp = 0.0;
        } else {
            if (fabs(sin(t)) < 1.0e-9 || fabs(cos(p)) < 1.0e-9) {
                rtp = 1.0e9;
                if (cos(t) * sin(p) < 0.0) rtp = -1.0e9;
            } else {
                rtp = -halfdim_[0] * cos(t) * sin(p) / sin(t) / sin(t) / cos(p) / cos(p);
            }
        }
    } else {
        tcrit = atan(halfdim_[1]/halfdim_[2]/sin(p));
        if (t <= tcrit) {
            rtp = 0.0;
        } else {
            if (fabs(sin(t)) < 1.0e-9 || fabs(sin(p)) < 1.0e-9) {
                rtp = 1.0e9;
                if (cos(t) * cos(p) < 0.0) rtp = -1.0e9;
            } else {
                rtp = halfdim_[1] * cos(t) * cos(p) / sin(t) / sin(t) / sin(p) / sin(p);
            }
        }
    }
    return rtp;

}

void tsquare::Box::doDilate(const Real dr)
{
    if (dr <= 0.0) {
        throw DataException("Box","doDilate",
                            "All dimensions must be positive numbers");
    }

    std::vector<Real> v(3);
    for (int i = 0; i < 3; i++) v[i] = halfdim_[i] * dr;
    setDimensions(v);

    std::complex<Real> cdr(dr,0.0);
    for (int i = 0; i <= nmax_; i++) {
        for (int j = -i; j <= i; j++) {
           a_[RDIM][i][j] *= cdr;
        }
    }

    return;
}

void tsquare::Box::computeVolume(void)
{
    volume_ = dim_[0] * dim_[1] * dim_[2];
    return;
}

void tsquare::Box::computeArea(void)
{
    area_ = 2.0 * ((dim_[0] * dim_[1]) + (dim_[0] * dim_[2]) + (dim_[1] * dim_[2]));
}

void tsquare::Box::computeBoundingspheres(void)
{
    using std::sqrt;

    inscribedsphere_.setRadius(halfdim_[0]);
    inscribedsphere_.setXr(0.0);
    inscribedsphere_.setYr(0.0);
    inscribedsphere_.setZr(0.0);

    enclosingsphere_.setXr(0.0);
    enclosingsphere_.setYr(0.0);
    enclosingsphere_.setZr(0.0);

    Real rad2 = 0.0;
    for (int i = 0; i < 3; i++) rad2 += (halfdim_[i] * halfdim_[i]);
    enclosingsphere_.setRadius(sqrt(rad2));

    return;
}

void tsquare::Box::computeDimensions(bool triaxialcalc)
{
    dim_[0] = 2.0 * halfdim_[0];
    dim_[1] = 2.0 * halfdim_[1];
    dim_[2] = 2.0 * halfdim_[2];
}

std::vector<Real> tsquare::Box::getI(void)
{
    std::vector<Real> itensor;
    itensor.clear();
    
    Real i11,i22,i33,i12,i23,i13;
    Real a = halfdim_[0];
    Real b = halfdim_[0];
    Real c = halfdim_[2];
    Real a2 = a * a;
    Real b2 = b * b;
    Real c2 = c * c;
    i11 = i22 = i33 = i12 = i23 = i13 = 0;

    i11 = (b2 + c2) / 12.0;
    i22 = (a2 + c2) / 12.0;
    i33 = (a2 + b2) / 12.0;

    itensor.push_back(i11);
    itensor.push_back(i22);
    itensor.push_back(i33);
    itensor.push_back(i12);
    itensor.push_back(i13);
    itensor.push_back(i23);

    return itensor;
}

void tsquare::Box::printProperties(void)
{
    using std::cout;
    using std::endl;

    cout << getType() << ":" << endl;
    cout << "a,b,c = " << halfdim_[0] << "," << halfdim_[1] << "," << halfdim_[2] << endl;
    cout << "volume = " << getVolume() << endl;
    cout << "surface area = " << getArea() << endl;
    cout << "normalized surface area = " << getNarea() << endl;
    cout << "diameter = " << getDiam() << endl;
    cout << "length, width, thickness = " << getLength() << ", "
         << getWidth() << ", " << getThickness() << endl;
    cout << "normalized length = " << getNlength() << endl;
    cout << "normalized width = " << getNwidth() << endl;
    cout << "Itrace = " << getItrace() << endl;
    cout << "Normalized Gaussian curvature = " << getNgc() << endl;
    cout.flush();
    cout << endl << endl;
    return;
}


tsquare::Box& tsquare::Box::operator=(const tsquare::Box& that)
{
    volume_ = that.getVolume();
    area_ = that.getArea();
    narea_ = that.getNarea();
    diam_ = that.getDiam();
    dim_ = that.getDim();
    ndim_ = that.getNdim();
    triaxialdim_ = that.getTriaxialdim();
    itrace_ = that.getItrace();
    ngc_ = that.getNgc();
    halfdim_[0] = 0.5 * that.getThickness();
    halfdim_[1] = 0.5 * that.getWidth();
    halfdim_[2] = 0.5 * that.getLength();
    ntheta_ = that.getNtheta();
    nphi_ = that.getNphi();
    surface_ = that.getSurface();
    r_ = that.getR();
    x_ = that.getX();
    y_ = that.getY();
    z_ = that.getZ();

    return *this;
}
