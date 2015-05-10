/**
@file
@brief Method definition for the Ellipsoid class.
*/
#include "Ellipsoid.h"

tsquare::Ellipsoid::Ellipsoid()
{
    adim_ = cdim_ = aspectratio_ = 1.0;
    volume_ = area_ = 0.0;
    narea_ = 0.0;
    ngc_ = diam_ = itrace_ = 0.0;
    dim_.clear();
    dim_.resize(3,1.0);
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
    enclosingsphere_.setRadius(1.0);
}

tsquare::Ellipsoid::Ellipsoid(Real a, Real c, bool verbose)
{
    verbose_ = verbose;
    adim_ = a;
    cdim_ = c;

    if (a <= 0.0 || c <= 0.0) {
        throw DataException("Ellipsoid","Ellipsoid",
                            "All dimensions must be positive numbers");
    }

    aspectratio_ = cdim_ / adim_;

    surface_.clear();
    r_.clear();
    xg_.clear();
    wg_.clear();
    dim_.clear();
    ndim_.clear();
    triaxialdim_.clear();
    dim_.resize(3,0.0);
    dim_[0] = dim_[1] = adim_;
    dim_[2] = cdim_;
    ndim_.resize(3,0.0);
    ndim_[0] = ndim_[1] = 1.0;
    ndim_[2] = aspectratio_;
    triaxialdim_ = dim_;

    inscribedsphere_.setXr(0.0);
    inscribedsphere_.setYr(0.0);
    inscribedsphere_.setZr(0.0);
    inscribedsphere_.setRadius(adim_);
    if (aspectratio_ < 1.0) inscribedsphere_.setRadius(cdim_);

    enclosingsphere_.setXr(0.0);
    enclosingsphere_.setYr(0.0);
    enclosingsphere_.setZr(0.0);
    enclosingsphere_.setRadius(cdim_);
    if (aspectratio_ < 1.0) enclosingsphere_.setRadius(adim_);
}

void tsquare::Ellipsoid::setA(const Real a)
{ 
    if (a <= 0) {
        throw DataException("Ellipsoid","setA",
                            "All dimensions must be positive numbers");
    }
    aspectratio_ *= (adim_ / a);
    adim_ = a;
}

void tsquare::Ellipsoid::setC(const Real c)
{
    if (c <= 0) {
        throw DataException("Ellipsoid","setC",
                            "All dimensions must be positive numbers");
    }
    aspectratio_ *= (c / cdim_);
    cdim_ = c;
}

void tsquare::Ellipsoid::setR(const Real theta, const Real phi)
{
    using std::pow;
    using std::sqrt;
    using std::sin;
    using std::cos;

    Real r1 = 1.0 / (sqrt(pow(cos(theta),2.0)+(pow(aspectratio_ * sin(theta),2.0))));

    // Check to see if this surface site is already listed.  If so,
    // then set the distance for that site. If not, add it now
 
    bool found = false;
    register int i;
    try {
        for (i = 0; (i < surface_.size() && !found); i++) {
            if (theta == surface_[i][0] && phi == surface_[i][1]) {
                found = true;
                r_.at(i) = r1;
            }
        }
    }
    catch (std::out_of_range &oor) {
        throw EOBException("Ellipsoid","setR","r_",r_.size(),i);
    }

    if (!found) {
        std::vector<Real> angles;
        angles.clear();
        angles.resize(2,0.0);
        angles[0] = theta;
        angles[1] = phi;
        surface_.push_back(angles);
        try {
            r_.at(surface_.size()-1) = r1;
        }
        catch (std::out_of_range &oor) {
            throw EOBException("Ellipsoid","setR","r_",
                                  r_.size(),surface_.size()-1);
        }
    }
    return;
}

void tsquare::Ellipsoid::setR(const unsigned int surfaceindex)
{
    using std::pow;
    using std::sqrt;
    using std::sin;
    using std::cos;

    Real t;
    try {
        t = surface_.at(surfaceindex).at(0);
    }
    catch (std::out_of_range &oor) {
        throw EOBException("Ellipsoid","setR","surface_",
                           surface_.size(),surfaceindex);
    }

    Real r1 = 1.0 / (sqrt(pow(cos(t),2.0)+(pow(aspectratio_ * sin(t),2.0))));
 
    if (surfaceindex < r_.size()) {
        r_[surfaceindex] = r1;
    } else {
        r_.push_back(r1);
    }
    return;
}

Real tsquare::Ellipsoid::getR(const Real t, const Real p)
{
    using std::pow;
    using std::sqrt;
    using std::sin;
    using std::cos;

    Real r1 = cdim_ / (sqrt(pow(cos(t),2.0)+(pow(aspectratio_ * sin(t),2.0))));
    return r1;
} 

void tsquare::Ellipsoid::doDilate(const Real dr)
{
    if (dr <= 0.0) {
        throw DataException("Ellipsoid","Ellipsoid",
                            "All dimensions must be positive numbers");
    }
    adim_ *= dr;
    cdim_ *= dr;

    std::complex<Real> cdr(dr,0.0);
    for (int i = 0; i <= nmax_; i++) {
        for (int j = -i; j <= i; j++) {
           a_[RDIM][i][j] *= cdr;
        }
    }

    return;
}
 
Real tsquare::Ellipsoid::getRt(Real t, const Real p)
{
    using std::pow;
    using std::sin;
    using std::cos;

    Real b = aspectratio_;
    if (cos(t) == 0.0 || sin(t) == 0.0 || b == 1.0) return 0.0;
    Real denom = pow(cos(t),2.0) + pow((b*sin(t)),2.0);
    Real rt = cdim_ * ((1.0 - (b*b)) * cos(t) * sin(t)) / (pow(denom,1.5));
    return rt;
}

Real tsquare::Ellipsoid::getRtt(Real t, const Real p)
{
    using std::pow;
    using std::sqrt;
    using std::sin;
    using std::cos;

    Real b = aspectratio_;
    Real rtt = cdim_ * (((1.0 - (b*b))*(4.0*(1.0 + (b*b))*cos(2.0*t) + 
         ((b*b) - 1.0)*(-5.0 + cos(4.0*t))))/
     (sqrt(2.0)*pow(1 + (b*b) - ((b*b) - 1.0)*cos(2.0*t),2.5)));
    return rtt;
}

void tsquare::Ellipsoid::computeVolume(void)
{
    volume_ = 4.0 * Pi * adim_ * adim_ * cdim_ / 3.0;
    return;
}

void tsquare::Ellipsoid::computeArea(void)
{
    using std::atanh;
    using std::sqrt;
    using std::asin;

    Real factor = 2.0 * Pi;
    Real e = 0.0;
    Real esquared = 0.0;

    if (cdim_ < adim_) {
        esquared = 1.0 - (cdim_ * cdim_ / adim_ / adim_);
        e = sqrt(esquared);
        area_ = factor * adim_ * adim_ * (1.0 + (((1.0 - esquared)/e) * atanh(e)));
    } else if (cdim_ > adim_) {
        esquared = 1.0 - (adim_ * adim_ / cdim_ / cdim_);
        e = sqrt(esquared);
        area_ = factor * adim_ * adim_ * (1.0 + ((cdim_ / adim_ / e) * asin(e)));
    } else {
        area_ = 4.0 * Pi * adim_ * adim_ * adim_ / 3.0;
    }
}

void tsquare::Ellipsoid::computeDimensions(bool triaxialcalc=true)
{
    dim_[0] = dim_[1] = adim_;
    dim_[2] = cdim_;
}

void tsquare::Ellipsoid::computeBoundingspheres(void)
{
    Real irad = adim_;
    if (aspectratio_ < 1.0) irad = cdim_;
    Real erad = cdim_;
    if (aspectratio_ < 1.0) erad = adim_;

    inscribedsphere_.setRadius(irad);
    inscribedsphere_.setXr(0.0);
    inscribedsphere_.setYr(0.0);
    inscribedsphere_.setZr(0.0);

    enclosingsphere_.setRadius(erad);
    enclosingsphere_.setXr(0.0);
    enclosingsphere_.setYr(0.0);
    enclosingsphere_.setZr(0.0);

    return;
}
std::vector<Real> tsquare::Ellipsoid::getI(void)
{
    std::vector<Real> itensor;
    itensor.clear();
    
    Real i11,i22,i33,i12,i23,i13;
    Real a = adim_;
    Real b = adim_;
    Real c = cdim_;
    Real a2 = a * a;
    Real b2 = b * b;
    Real c2 = c * c;
    i11 = i22 = i33 = i12 = i23 = i13 = 0;

    i11 = (a2 + c2);
    i22 = (b2 + c2);
    i33 = (a2 + b2);

    itensor.push_back(i11);
    itensor.push_back(i22);
    itensor.push_back(i33);
    itensor.push_back(i12);
    itensor.push_back(i13);
    itensor.push_back(i23);

    return itensor;
}

void tsquare::Ellipsoid::printProperties()
{
    using std::cout;
    using std:: endl;

    cout << getType() << ":" << endl;
    cout << "a,c = " << getA() << "," << getC() << endl;
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

tsquare::Ellipsoid& tsquare::Ellipsoid::operator=(const tsquare::Ellipsoid& that)
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
    adim_ = that.getA();
    cdim_ = that.getC();
    aspectratio_ = that.getAspectratio();
    ntheta_ = that.getNtheta();
    nphi_ = that.getNphi();
    surface_ = that.getSurface();
    r_ = that.getR();
    x_ = that.getX();
    y_ = that.getY();
    z_ = that.getZ();

    return *this;
}
