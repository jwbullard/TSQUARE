/**
@file
@brief Main program for differential geometry toolbox
*/
#include "tsquare.h"

using namespace tsquare;

int main (int argc, char *argv[])
{
    using std::cout;
    using std::cin;
    using std::endl;

    int choice = 1;

    while (choice != 0) {
        cout << "MAIN MENU:" << endl << endl;
        cout << "    0. Exit" << endl;
        cout << "    1. Differential geometry of an individual particle" << endl;
        cout << "    2. Shape statistics on a collectio nof particles" << endl << endl;
        cout << "Enter a menu number: ";
        cin >> choice;
        cout << choice << endl << endl;
        cout.flush();

        switch (choice) {
            case 0:         // Exit the program
                break;
            case 1:
                diffgeom();
                break;
            case 2:
                shapestats(argc,argv);
                break;
            default:
                break;
        }
    }

    cout << "Exiting now" << endl << endl;
}

int diffgeom(void)
{
    using std::cout;
    using std::cin;
    using std::endl;
    using std::sin;
    using std::cos;

    std::string buff;
    int npoints = 3;
    int choice = 1;
    int pointid = 0;
    int subchoice = 1;
    int dering_start_deg = 0;
    int numgausspoints = 106;
    int nmax = 1;
    std::vector<Real> pd1,pd2,unorm;
    Real alpha,beta,gamma,theta,phi,rval,isphere,rg2sphere;
    Real excludedvolume,convexity,volume,area,angularity,rotateangle;
    Real digitalResolution;
    bool gaussian = false;
    bool firsttime = true;
    bool render_hull,render_inscribed_sphere,render_enclosing_sphere;
    int surface_color_scheme;
    int digitalXdim,digitalYdim,digitalZdim;
    std::vector<Real> itensor,rg2tensor;
    std::vector<std::vector<std::vector<int> > > boundingBox;
    Star *s;
    char yes_no;
    std::string fname,outfilename,msg;
    std::string vrmlname,digitizedname,shname;

    // We are going to read in just one test particle and perform
    // unit tests on the Particle class

    s = NULL;

    // Menu for choices of properties to calculate

    cout << endl << endl;
    std::vector<Real> firstvector,secondvector,diffvector;
    firstvector.clear();
    firstvector.resize(3,0.0);
    secondvector.clear();
    secondvector.resize(3,0.0);
    firstvector[0] = 5.0;
    secondvector[0] = 1.0;
    firstvector[1] = 0.0;
    secondvector[1] = 3.0;
    firstvector[2] = 9.0;
    secondvector[2] = -3.0;
    diffvector = firstvector - secondvector;

    while (choice != 0) {
        cout << " 0. Back to main menu" << endl;
        cout << " 1. Read particle SH coefficients from file" << endl;
        cout << " 2. Calculate SH coefficients from binary image" << endl;
        cout << " 3. Print properties" << endl;
        cout << " 4. Set surface points" << endl;
        cout << " 5. Print radius of a surface point" << endl;
        cout << " 6. Print normal of a surface point" << endl;
        cout << " 7. Print curvature(s) of a surface point" << endl;
        cout << " 8. Print curvature(s) along a path" << endl;
        cout << " 9. Reduce ringing" << endl;
        cout << "10. Restore ringing" << endl;
        cout << "11. Rotate the particle" << endl;
        cout << "12. Compute convexity" << endl;
        cout << "13. Find stable rest position" << endl;
        cout << "14. Compute volume, area, and roundness" << endl;
        cout << "15. Find bounding spheres" << endl;
        cout << "16. Find length, width, and thickness" << endl;
        cout << "17. Compute shape factors" << endl;
        cout << "18. Compute moment of inertia and gyration tensor components" << endl;
        cout << "19. Output current SH coefficients to file" << endl;
        cout << "20. Output digitized image of particle" << endl;
        cout << "21. Output particle as VRML" << endl;

        cout << "Enter a menu number: ";
        cin >> choice;
        cout << choice << endl << endl;
        cout.flush();
        switch (choice) {
            case 0:  // Exit the program
                break;
            case 1:  // Read SH coefficients
                cout << "Enter name for file containing particle's SH coefficients: ";
                cin >> fname;
                cout << endl << "******" << fname << " ******" << endl;
                cout.flush();
                if (!firsttime) {
                    delete s;
                    firsttime = false;
                }
                try {
                    s = new Star(fname,true);
                }
                catch (std::bad_alloc &ba) {
                    msg = "Bad memory allocation in Particle constructor";
                    DataException("diffgeom","main",msg);
                    exit(1);
                }
                catch (FileException fex) {
                    fex.printException();
                    exit(1);
                }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                break;
            case 2:  // Calculate SH coefficients from binary image
                cout << "Enter name of file containing binary image:  ";
                cin >> fname;
                cout << "Enter name of output anm file to create:  ";
                cin >> outfilename;
                cout << "Enter number of Gaussian points to use:  ";
                cin >> npoints;
                cout << "Enter maximum degree to use for SH expansion:  ";
                cin >> nmax;
                cout << "Use surface mesh to calculate?  ";
                cin >> yes_no;
                if (!firsttime) {
                    delete s;
                    firsttime = false;
                }
                try {
                    s = new Star();
                }
                catch (std::bad_alloc &ba) {
                    msg = "Bad memory allocation in Particle constructor";
                    DataException("diffgeom","main",msg);
                    exit(1);
                }
                if (yes_no == 'N' || yes_no == 'n') {
                    s->calculateSHcoeffs(fname,outfilename,npoints,nmax,false);
                } else {
                    s->calculateSHcoeffs(fname,outfilename,npoints,nmax,true);
                }
                break;
            case 3:  // Print properties
                s->printProperties();
                break;
            case 4:  // Set points
                cout << "  1. Set Gaussian quadrature points" << endl;
                cout << "  2. Set regular array of points" << endl;
                cout << "  Enter a menu number: ";
                cin >>  subchoice;
                cout << subchoice << endl;
                cout << "  How many points? ";
                cin >> npoints;
                cout << npoints << endl;
                cout << "  Maximum degree of SH expansion ";
                cout << "(not greater than " << s->getNmax_from_anmfile() << "):";
                cin >> nmax;
                nmax = (nmax < s->getNmax_from_anmfile()) ? nmax : s->getNmax_from_anmfile();
                cout << nmax << endl;
                cout.flush();
                try {
                    switch (subchoice) {
                        case 1:
                            s->setSurface(npoints,nmax,true);
                            break;
                        default:
                            s->setSurface(npoints,nmax,false);
                            break;
                    }
                }
                catch (EOBException ex) {
                    ex.printException();
                    exit(1);
                }
                cout << endl << "nmax = " << nmax << ", mesh points = "
                     << npoints << ":" << endl;
                cout << "    Volume = " << std::setprecision(12) << s->getVolume() << endl;
                cout << "    Area = " << std::setprecision(12) << s->getArea() << endl << endl;
                break;
            case 5:  // Print radius
                cout << "  Pick theta: (0 to 180 degrees): ";
                cin >> theta;
                cout << "    " << theta << endl;
                cout << "  Pick phi: (0 to 360 degrees): ";
                cin >> phi;
                cout << "    " << phi << endl;
                theta *= (1.0/Rad2Deg);
                phi *= (1.0/Rad2Deg);
                rval = s->getR(theta,phi);
                cout << endl << "  r( " << Rad2Deg * theta << ","
                     << Rad2Deg * phi << ") = " << rval << endl;
                cout << endl << "  Cartesian coordinates: ("
                     << rval * sin(theta) * cos(phi) << ", ";
                cout << rval * sin(theta) * sin(phi) << ", "
                     << rval * cos(theta) << ")" << endl;
                cout << "Rt = " << s->getRt(theta,phi) << endl;
                cout << "Rp = " << s->getRt(theta,phi) << endl;
                cout << endl;
                break;
            case 6:  // Print normal vector
                cout << "  Theta = ";
                cin >> theta;
                cout << endl << "  Phi = ";
                cin >> phi;
                cout << endl << endl;
                try { unorm = s->getNormal(theta*Deg2Rad,phi*Deg2Rad); }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                cout << "  Point (" << theta << ","
                     << phi << "): n = ("
                     << std::setprecision(12) << unorm[0] << ","
                     << std::setprecision(12) << unorm[1] << ","
                     << std::setprecision(12) << unorm[2] << ")" << endl << endl;
                break;
            case 7:  // Print curvature measures at a single point
                cout << "  Pick a point (0 to "
                     << (npoints - 1) * (npoints - 1) - 1 << ")"
                     << endl;
                cin >> pointid;
                try { printCurvature(pointid,npoints,s,false); }
                catch (EOBException ex) {
                    ex.printException();
                    exit(1);
                }
                catch (FloatException flex) {
                    flex.printException();
                    exit(1);
                }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                break;
            case 8:  // Print curvature measures along a path
                try { printCurvatureOnAPath(s); }
                catch (EOBException ex) {
                    ex.printException();
                    exit(1);
                }
                catch (FloatException flex) {
                    flex.printException();
                    exit(1);
                }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                catch (FileException fex) {
                    fex.printException();
                    exit(1);
                }
                break;
            case 9:  // Reduce the amount of ringing on the surface (Gibbs phenomenon)
                pointid = -1;
                try {
                    cout << "  Choose Lanczos filter (0) or Hanning filter (1): ";
                    cin >> pointid;
                    cout << pointid << endl;
                    cout << "    Enter degree to begin de-ringing: ";
                    cin >> dering_start_deg;
                    cout << dering_start_deg << endl;
                    if (dering_start_deg < 0) {
                        msg = "De-ring start cannot be less than 0.";
                        throw DataException("diffgeom","main",msg);
                    } else if (dering_start_deg >= s->getNmax()) {
                        msg = "De-ring start degree exceeds maximum";
                        throw DataException("diffgeom","main",msg);
                    }
                    if (pointid == 0) {
                        cout << "    Enter exponent (typically 1 or 2): ";
                        cin >> subchoice;
                        cout << subchoice << endl;
                        cout.flush();
                        s->doLanczos(subchoice,dering_start_deg);
                    } else if (pointid == 1) {
                        cout.flush();
                        s->doHanning(dering_start_deg);
                    } else {
                        throw DataException("diffgeom","main",
                                            "Unrecognized de-ringing method");
                    }
                }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                break;
            case 10:  // Cancel the effects of option 7 (reduce ringing)
                pointid = -1;
                try {
                    cout << "  Undo Lanczos filter (0) or Hanning filter (1)? ";
                    cin >> pointid;
                    cout << endl << pointid << endl;
                    cout << "    Enter degree to begin undoing: ";
                    cin >> dering_start_deg;
                    cout << dering_start_deg << endl;
                    if (dering_start_deg < 0) {
                        msg = "De-ring start cannot be less than 0.";
                        throw DataException("diffgeom","main",msg);
                    } else if (dering_start_deg >= s->getNmax()) {
                        msg = "De-ring start degree exceeds maximum";
                        throw DataException("diffgeom","main",msg);
                    }
                    if (pointid == 0) {
                        cout << "    Enter exponent (typically 1 or 2): ";
                        cin >> subchoice;
                        s->undoLanczos(subchoice,dering_start_deg);
                    } else if (pointid == 1) {
                        s->undoHanning(dering_start_deg);
                    } else {
                        throw DataException("diffgeom","main",
                                            "Unrecognized de-ringing method");
                    }
                }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                break;
            case 11:  // Rotate the particle
                cout << "NOTE:  For each angle, a positive value "
                     << "means a COUNTERCLOCKWISE" << endl;
                cout << "       rotation OF THE PARTICLE..." << endl;
                cout << endl;
                cout << "       First angle: alpha (deg) about z-axis: ";
                cin >> alpha;
                cout << "    " << alpha << endl;
                cout << "       Second angle: beta (deg) about NEW y-axis: ";
                cin >> beta;
                cout << "    " << beta << endl;
                cout << "       Third angle: gamma (deg) about NEW z-axis: ";
                cin >> gamma;
                cout << "    " << gamma << endl;
                cout << "       Rotating particle by (" << alpha << ","
                     << beta << "," << gamma << ")... ";
                cout.flush();
                alpha *= Deg2Rad;
                beta *= (-1.0 * Deg2Rad);
                gamma *= Deg2Rad;
                s->doRotate(alpha,beta,gamma);
                beta *= -1.0;
                cout << "Done!" << endl;
                cout << endl;
                break;
            case 12:  // Compute convexity of the particle
                try { convexity = s->getConvexity(); }
                catch (FileException fex) {
                    fex.printException();
                    exit(1);
                }
                catch (FloatException flex) {
                    flex.printException();
                    exit(1);
                }
                catch (EOBException ex) {
                    ex.printException();
                    exit(1);
                }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                cout << "Convexity = " << convexity << endl;
                break;
            case 13:  // Find stable resting position
                try { s->getRestingConfiguration(); }
                catch (FileException fex) {
                    fex.printException();
                    exit(1);
                }
                catch (FloatException flex) {
                    flex.printException();
                    exit(1);
                }
                catch (EOBException ex) {
                    ex.printException();
                    exit(1);
                }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                cout << "Convexity = " << convexity << endl;
                break;
            case 14:  // Compute particle volume, area, and angularity
                cout << "Computing roundness factors... " << endl;
                cout.flush();
                try {
                    s->computeRoundness(ALLROUNDNESSES);
                }
                catch (EOBException ex) {
                    ex.printException();
                    exit(1);
                }
                catch (FloatException flex) {
                    flex.printException();
                    exit(1);
                }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                cout << "              Volume = " << std::setprecision(12) << s->getVolume() << endl;
                cout << "        Surface area = " << std::setprecision(12) << s->getArea() << endl;
                cout << "        DP Roundness = " << std::setprecision(12) << s->getRoundness(DOTPRODUCT) << endl;
                cout << "    Wadell Roundness = " << std::setprecision(12) << s->getRoundness(WADELL) << endl;
                cout.flush();
                break;
            case 15:  // Compute bounding spheres
                cout << "Calculating bounding spheres now..." << endl;
                cout.flush();
                try { s->computeBoundingspheres(); }
                catch (EOBException ex) {
                    ex.printException();
                    exit(1);
                }
                cout.flush();
                break;
            case 16:  // Get the three orthogonal dimensions of the particle
                cout << "Length = " << std::setprecision(12) << s->getLength() << endl;
                cout << "Width = " << std::setprecision(12) << s->getWidth() << endl;
                cout << "Thickness = " << std::setprecision(12) << s->getThickness() << endl << endl;
                break;
            case 17:  // Compute shape factors
                cout << "Computing dimensions... " << endl;
                cout.flush();
                try {
                    s->computeDimensions(true);
                    cout << "Length = " << std::setprecision(12) << s->getLength() << endl;
                    cout << "Width = " << std::setprecision(12) << s->getWidth() << endl;
                    cout << "Thickness = " << std::setprecision(12) << s->getThickness() << endl << endl;
                    cout << "Wadell Sphericity = "
                         << std::setprecision(12) << s->getWadellsphericity() << endl;
                    cout << "Ratio of bounding spheres = "
                         << std::setprecision(12) << s->getBoundingsphereratio() << endl;
                    cout << "Corey Shape Factor = " << std::setprecision(12) << s->getCSF() << endl;
                    cout << "Aschenbrenner Shape Factor = " << std::setprecision(12) << s->getASF() << endl;
                    cout << "Maximum Projection Sphericity = " << std::setprecision(12) << s->getMPS() << endl;
                    cout << "E Shape Factor = " << std::setprecision(12) << s->getESF() << endl;
                    cout << "Hoffman's Shape Entropy = "
                         << std::setprecision(12) << s->getShapeentropy() << endl;
                }
                catch (DataException dex) {
                    dex.printException();
                    exit(1);
                }
                catch (EOBException ex) {
                    ex.printException();
                    exit(1);
                }
                catch (FloatException flex) {
                    flex.printException();
                    exit(1);
                }
                break;
            case 18:  // Compute moment of inertia and gyration tensors
                s->getI(itensor,rg2tensor,isphere,rg2sphere);
                cout << endl;
                cout << "      | " << itensor[0] << " " << itensor[3] << " " << itensor[4] << " |" << endl;
                cout << "I =   | " << itensor[3] << " " << itensor[1] << " " << itensor[5] << " |" << endl;
                cout << "      | " << itensor[4] << " " << itensor[5] << " " << itensor[2] << " |" << endl;
                cout << endl;
                cout << "Tr(I) = " << itensor[0] + itensor[1] + itensor[2] << endl;
                cout << "**********************************************************************" << endl;
                cout << endl;
                cout << "      | " << rg2tensor[0] << " " << rg2tensor[3] << " " << rg2tensor[4] << " |" << endl;
                cout << "Rg2 = | " << rg2tensor[3] << " " << rg2tensor[1] << " " << rg2tensor[5] << " |" << endl;
                cout << "      | " << rg2tensor[4] << " " << rg2tensor[5] << " " << rg2tensor[2] << " |" << endl;
                cout << endl;
                cout << "Tr(Rg2) = " << rg2tensor[0] + rg2tensor[1] + rg2tensor[2] << endl;
                cout << "**********************************************************************" << endl;
                cout << "I for equivalent sphere = " << isphere << endl;
                cout << "Rg2 for equivalent sphere = " << rg2sphere << endl;
                cout << "**********************************************************************" << endl;
                if (isphere > 0.0) {
                    itensor[0] /= isphere;
                    itensor[1] /= isphere;
                    itensor[2] /= isphere;
                    itensor[3] /= isphere;
                    itensor[4] /= isphere;
                    itensor[5] /= isphere;
                }
                if (rg2sphere > 0.0) {
                    rg2tensor[0] /= rg2sphere;
                    rg2tensor[1] /= rg2sphere;
                    rg2tensor[2] /= rg2sphere;
                    rg2tensor[3] /= rg2sphere;
                    rg2tensor[4] /= rg2sphere;
                    rg2tensor[5] /= rg2sphere;
                }
                cout << endl;
                cout << "              | " << itensor[0] << " " << itensor[3] << " " << itensor[4] << " |" << endl;
                cout << "I/Isphere =   | " << itensor[3] << " " << itensor[1] << " " << itensor[5] << " |" << endl;
                cout << "              | " << itensor[4] << " " << itensor[5] << " " << itensor[2] << " |" << endl;
                cout << endl;
                cout << "Tr(I/Isphere) = " << itensor[0] + itensor[1] + itensor[2] << endl;
                cout << "**********************************************************************" << endl;
                cout << endl;
                cout << "                | " << rg2tensor[0] << " " << rg2tensor[3] << " " << rg2tensor[4] << " |" << endl;
                cout << "Rg2/Rg2sphere = | " << rg2tensor[3] << " " << rg2tensor[1] << " " << rg2tensor[5] << " |" << endl;
                cout << "                | " << rg2tensor[4] << " " << rg2tensor[5] << " " << rg2tensor[2] << " |" << endl;
                cout << endl;
                cout << "Tr(Rg2/Rg2sphere) = " << rg2tensor[0] + rg2tensor[1] + rg2tensor[2] << endl;
                cout << endl;
                cout.flush();
                break;

            case 19:  // Output SH coefficients
                cout << "Output file name: ";
                cin >> shname;
                cout << endl << shname << endl;
                try {
                    s->createSHFile(shname);
                }
                catch (FileException fex) {
                    fex.printException();
                    exit(1);
                }
                break;

            case 20:  // Output digitized image
                cout << "Digitized image file root name: ";
                cin >> digitizedname;
                cout << endl << digitizedname << endl;
                cout << "Voxel resolution: ";
                cin >> digitalResolution;
                cout << endl << digitalResolution << endl;
                if (digitalResolution < 0.1) digitalResolution = 0.1;
                try {
                    s->createDigitizedImage(digitizedname,digitalResolution);
                }
                catch (FileException fex) {
                    fex.printException();
                    exit(1);
                }
                break;

            case 21:  // Output VRML
                render_inscribed_sphere = false;
                render_enclosing_sphere = false;
                render_hull = false;
                cout << "    VRML file root name: ";
                cin >> vrmlname;
                cout << vrmlname << endl;
                cout << "Surface color map:" << endl;
                cout << "    Maximum curvature roundness (" << WADELL << ")" << endl;
                cout << "    Dot product roundness (" << DOTPRODUCT << ")" << endl;
                cout << "    Uniform color (" << ALLROUNDNESSES << ")" << endl;
                cout << "Selection: ";
                cin >> surface_color_scheme;
                cout << surface_color_scheme << endl;
                cout << "Render maximimum inscribed sphere (1) or not (0): ";
                cin >> subchoice;
                cout << subchoice << endl;
                if (subchoice == 1) render_inscribed_sphere = true;
                cout << "Render minimum enclosing sphere (1) or not (0): ";
                cin >> subchoice;
                cout << subchoice << endl;
                if (subchoice == 1) render_enclosing_sphere = true;
                cout << "Render convex hull (1) or not (0): ";
                cin >> subchoice;
                cout << subchoice << endl;
                cout << "Angle to rotate about y axis (degrees): ";
                cin >> rotateangle;
                cout << rotateangle << endl;
                if (subchoice == 1) render_hull = true;
                try { s->createVRML(vrmlname,render_inscribed_sphere,
                                  render_enclosing_sphere, render_hull,
                                  surface_color_scheme,rotateangle); }
                catch (FileException fex) {
                    fex.printException();
                    exit(1);
                }
                catch (EOBException ex) {
                    ex.printException();
                    exit(1);
                }
                cout << "Done!" << endl;
                break;
            default:
                DataException dex("diffgeom","main","Unrecognized choice in main menu");
                dex.printException();
                break;
        }
    }

    if (!firsttime) delete s;
}

Real getExcludedvolume(Particle *p)
{
    using std::cout;
    using std::endl;
    using std::pow;
    using std::sqrt;
    using std::sin;
    using std::cos;
    using std::acos;
    using std::asin;

    Real length = p->getLength();
    Real boxhalflength = length;
    Real Rtheta1,Rphi1,Rtheta2,Rphi2,Rmag;
    std::vector<Real> R;
    R.resize(3,0.0);
    std::vector<Real> theta1,theta2,phi1,phi2;
    theta1.resize(2,0.0);
    theta2.resize(2,0.0);
    phi1.resize(2,0.0);
    phi2.resize(2,0.0);
    Real halfpi = 0.5 * Pi;
    Real alpha,beta,gamma;
    Real cosbeta,sinbeta;
    std::vector<int> si1,si2;
    std::vector<Real> r2,point;
    Real pmag,th,ph,prad;
    bool done;
    Real ev = 0.0;
    // RanGen rg(-142234);

    // First set the number of points on the surface to 50 in each direction
    // p->setSurface(50,false);

    // Next, create an exact copy of the particle
    Particle *q = p;

    // The original particle, p, will remain stationary, but the copy particle
    // will be rotated and translated

    Real xmin = -boxhalflength;
    Real ymin = xmin;
    Real zmin = xmin;
    Real xmax = boxhalflength;
    Real ymax = xmax;
    Real zmax = xmax;
    Real boxvolume = pow((xmax - xmin),3.0);
    
    int numsamples = 100;
    Real tally = 0.0;
    for (int i = 0; i < numsamples; i++) {
        cout << "Sample event: " << i << endl;
    
        // R[0] = -(boxhalflength) + ((rg.Ran3()) * (2.0 * boxhalflength));
        // R[1] = -(boxhalflength) + ((rg.Ran3()) * (2.0 * boxhalflength));
        // R[2] = -(boxhalflength) + ((rg.Ran3()) * (2.0 * boxhalflength));

        // This is dummy information to avoid library linking error with rg object
        // that I still don't understand

        R[0] = 1.0;
        R[1] = 0.0;
        R[2] = 0.0;
        Rmag = sqrt((R[0]*R[0]) + (R[1]*R[1]) + (R[2]*R[2]));
        Rtheta1 = 0.5 * Pi;
        Rphi1 = 0.0;
        if (Rmag > 0) Rtheta1 = acos(R[2]/Rmag);
        if (Rtheta1 > 0.0 && Rtheta1 < Pi) Rphi1 = acos(R[0]/(Rmag * sin(Rtheta1)));
        Rphi2 = Rphi1 + Pi;
        if (Rphi2 > 2.0 * Pi) Rphi2 -= (2.0 * Pi);
        Rtheta2 = Pi - Rtheta1;
        cout << "    Radius vector is R = " << Rmag << "; theta = "
             << Rtheta1 * 180.0/Pi << ", phi = " << Rphi1 * 180.0/Pi << endl;
        // beta = Pi * rg.Ran3();
        // alpha = 2.0 * Pi * rg.Ran3();
        // gamma = 2.0 * Pi * rg.Ran3();
    
        // These are dummy values to avoid libary linker error to rg object above
    
        beta = 0.5 * Pi;
        alpha = Pi;
        gamma = Pi;
        cout << "Rotating (" << beta*180.0/Pi << "," << gamma*180.0/Pi << ","
             << alpha*180.0/Pi << ")... ";
        cout.flush();
        q->doRotate(beta,gamma,alpha);
        cout << "Done!" << endl;
        cout.flush();
        cout << "    Finding angle ranges" << endl;
        
        // Angle range for particle 1 (stationary)
        theta1[0] = Rtheta1 - halfpi;
        theta1[1] = Rtheta1 + halfpi;
        phi1[0] = Rphi1 - halfpi;
        phi1[1] = Rphi1 + halfpi;
        cout << "        Angle range for particle 1 is (" << theta1[0] * 180.0/Pi << ","
             << theta1[1] * 180.0/Pi << ") ";
        cout << " and (" << phi1[0] * 180.0/Pi << "," << phi1[1] * 180.0/Pi << ")" << endl;
    
        // Angle range for particle 2 (moving)
        theta2[0] = Rtheta2 - halfpi;
        theta2[1] = Rtheta2 + halfpi;
        phi2[0] = Rphi2 - halfpi;
        phi2[1] = Rphi2 + halfpi;
        cout << "        Angle range for particle 2 is ("
             << theta2[0] * 180.0/Pi << "," << theta2[1] * 180.0/Pi << ") ";
        cout << " and (" << phi2[0] * 180.0/Pi << "," << phi2[1] * 180.0/Pi
             << ")" << endl;
        cout << "    Looping over all surface sites" << endl;
        
        q->getSurfacesubset(theta2[0],theta2[1],phi2[0],phi2[1],si2);
        done = false;
        try {
            for (int k = 0; (k < si2.size() && !done); k++) {
                r2 = q->getX(si2[k]);
                getVectorsum(r2,R,point);
                pmag = sqrt((point[XDIM]*point[XDIM]) +
                            (point[YDIM]*point[YDIM]) +
                            (point[ZDIM]*point[ZDIM]));
                th = acos(point[ZDIM] / pmag);
                if (th > 0.0 && th < Pi) {
                    if (point[XDIM] * point[XDIM] > 0.0) {
                        ph = acos(point[XDIM] / pmag / sin(th));
                    } else {
                        ph = asin(point[YDIM] / pmag / sin(th));
                    }
                } else {
                    ph = 0.0;
                }
                prad = p->getR(th,ph);
                if (pmag < prad) { cout << "Hit!" << endl; tally += 1.0; done = true; }
            }
        }
        catch (EOBException ex) { throw ex; }
    }
    
    return (tally * boxvolume / ((Real)numsamples));
}

void getVectorsum(std::vector<Real> v1, std::vector<Real> v2, std::vector<Real> &result)
{
    result.clear();
    result.resize(v1.size(),0.0);
    if (v1.size() != v2.size()) return;

    for (int i = 0; i < v1.size(); i++) {
        result[i] = (v1[i] + v2[i]);
    }
    return;
}

void printCurvature(int pointid, int npoints, Particle *p, bool all=false)
{
    using std::cout;
    using std::endl;

    std::vector<Real> pd1,pd2;
    try {
        if (all) {
            for (int pointid = 0; pointid < 0.5 * ((npoints - 1) * (npoints - 1)); pointid++) {
                cout << "(" << Rad2Deg * p->getTheta(pointid) << ","
                     << Rad2Deg * p->getPhi(pointid) << "):" << endl;
                if (!(p->isUmbilic(pointid))) {
                    cout << "    k1 = " << std::setprecision(12)
                         << p->getPrincipalcurvature(0,pointid) << endl;
                    cout << "    k2 = " << std::setprecision(12)
                         << p->getPrincipalcurvature(1,pointid) << endl;
                    pd1 = p->getPrincipaldirection(0,pointid);
                    pd2 = p->getPrincipaldirection(1,pointid);
                    cout << "    Principal directions = (" << std::setprecision(12) << pd1[0]
                         << "," << std::setprecision(12) << pd1[1]
                         << "), (" << std::setprecision(12) << pd2[0] << "," << std::setprecision(12)
                         << pd2[1] << ")" << endl;
                    cout << "     H = " << std::setprecision(12) << p->getH(pointid) << endl;
                    cout << "     K = " << std::setprecision(12) << p->getK(pointid) << endl;
                } else if (!(p->isParabolicumbilic(pointid))) {
                    int kval = 0;
                    cout << "  Elliptic umbilic point: ";
                    cout << "    k1 = " << std::setprecision(12) << p->getEllipticumbiliccurvature(kval,pointid)
                         << endl;
                    cout << "    k2 = " << std::setprecision(12) << p->getEllipticumbiliccurvature(kval,pointid)
                         << endl;
                    cout << "     H = " << std::setprecision(12) << p->getEllipticumbiliccurvature(kval,pointid)
                         << endl;
                    cout << "     K = "
                         << std::setprecision(12) << pow((p->getEllipticumbiliccurvature(kval,pointid)),2.0) << endl;
                } else {
                    cout << " Parabolic elliptic point!";
                }
                cout << endl;
            }
        } else {
            cout << "(" << Rad2Deg * p->getTheta(pointid) << ","
                 << Rad2Deg * p->getPhi(pointid) << "):" << endl;
            if (!(p->isUmbilic(pointid))) {
                cout << "    k1 = " << std::setprecision(12) << p->getPrincipalcurvature(0,pointid) << endl;
                cout << "    k2 = " << std::setprecision(12) << p->getPrincipalcurvature(1,pointid) << endl;
                pd1 = p->getPrincipaldirection(0,pointid);
                pd2 = p->getPrincipaldirection(1,pointid);
                cout << "    Principal directions = (" << std::setprecision(12) << pd1[0]
                     << "," << std::setprecision(12) << pd1[1]
                     << "), (" << std::setprecision(12) << pd2[0] << "," << std::setprecision(12) << pd2[1] << ")" << endl;
                cout << "     H = " << std::setprecision(12) << p->getH(pointid) << endl;
                cout << "     K = " << std::setprecision(12) << p->getK(pointid) << endl;
            } else if (!(p->isParabolicumbilic(pointid))) {
                int kval = 0;
                cout << "  Elliptic umbilic point: ";
                cout << "    k1 = " << std::setprecision(12) << p->getEllipticumbiliccurvature(kval,pointid)
                     << endl;
                cout << "    k2 = " << std::setprecision(12) << p->getEllipticumbiliccurvature(kval,pointid)
                     << endl;
                cout << "     H = " << std::setprecision(12) << p->getEllipticumbiliccurvature(kval,pointid)
                     << endl;
                cout << "     K = "
                     << std::setprecision(12) << pow((p->getEllipticumbiliccurvature(kval,pointid)),2.0) << endl;
            } else {
                cout << " Parabolic elliptic point!";
            }
            cout << endl;
        }
        return;
    }
    catch (EOBException ex) { throw ex; }
    catch (FloatException flex) { throw flex; }
    catch (DataException dex) { throw dex; }
}

void printCurvatureOnAPath(Particle *p)
{
    using std::cout;
    using std::cin;
    using std::endl;

    int pathchoice = -1;
    Real ctheta,cphi,astart,aend,ainc,amin,amax,angle;
    Real cthetamin,cthetamax,cphimin,cphimax;
    std::string frname = "mypath";
    std::vector<Real> pd1,pd2;

    cout << "  Enter root name of output files to create: ";
    cin >> frname;
    cout << endl;

    cout << "  Constant theta (0) or phi (1) path? ";
    cin >> pathchoice;
    cout << endl;

    if (pathchoice != 0 && pathchoice != 1) {
        throw DataException("diffgeom","printCurvatureOnAPath",
                            "Unrecognized path choice");
    }

    cthetamin = 1.0e-5;
    cthetamax = Pi - cthetamin;
    cphimin = 0.0;
    cphimax = 2.0 * Pi;
    astart = aend = ainc = -1.0;
    
    switch (pathchoice) {
        case 0:
            amin = 0.0;
            amax = 2.0 * Pi;
            ctheta = -1.0;
            cout << "  What is the constant value of theta for path (degrees)? ";
            cin >> ctheta;
            ctheta *= Deg2Rad;
            cout << endl;
            if (ctheta < cthetamin || ctheta > cthetamax) {
                throw DataException("diffgeom","printCurvatureOnAPath",
                                    "Bad value for path theta");
            }
            cout << "  Enter beginning angle in degrees: ";
            cin >> astart;
            astart *= Deg2Rad;
            cout << endl;
            if (astart < amin || astart > amax) {
                throw DataException("diffgeom","printCurvatureOnAPath",
                                    "Bad value for path theta start");
            }
            cout << "  Enter ending angle in degrees: ";
            cin >> aend;
            aend *= Deg2Rad;
            cout << endl;
            if (aend < astart || aend > amax) {
                throw DataException("diffgeom","printCurvatureOnAPath",
                                    "Bad value for path theta end");
            }
            cout << "  Enter angle increment along path in degrees: ";
            cin >> ainc;
            ainc *= Deg2Rad;
            cout << endl;
            if (ainc < 0.0) {
                throw DataException("diffgeom","printCurvatureOnAPath",
                                    "Negative value for path theta increment");
            }
            try { printCurvatureOnALatitude(p,ctheta,astart,aend,ainc,frname); }
            catch (EOBException ex) { throw ex; }
            catch (FloatException flex) { throw flex; }
            catch (DataException dex) { throw dex; }
            catch (FileException fex) { throw fex; }
            break;
        default:
            amin = 1.0e-5;
            amax = Pi - 1.0e-5;
            cphi = -1.0;
            cout << "  What is the constant value of phi for path (degrees)? ";
            cin >> cphi;
            cphi *= Deg2Rad;
            cout << endl;
            if (cphi < cphimin || cphi > cphimax) {
                throw DataException("diffgeom","printCurvatureOnAPath",
                                    "Bad value for path phi");
            }
            cout << "  Enter beginning angle in degrees: ";
            cin >> astart;
            astart *= Deg2Rad;
            cout << endl;
            if (astart < amin || astart > amax) {
                throw DataException("diffgeom","printCurvatureOnAPath",
                                    "Bad value for path phi start");
            }
            cout << "  Enter ending angle in degrees: ";
            cin >> aend;
            aend *= Deg2Rad;
            cout << endl;
            if (aend < astart || aend > amax) {
                throw DataException("diffgeom","printCurvatureOnAPath",
                                    "Bad value for path phi end");
            }
            while (ainc < 0.0) {
                cout << "  Enter angle increment along path in degrees: ";
                cin >> ainc;
                ainc *= Deg2Rad;
                cout << endl;
            }
            if (ainc < 0.0) {
                throw DataException("diffgeom","printCurvatureOnAPath",
                                    "Negative value for path phi increment");
            }
            try { printCurvatureOnALongitude(p,cphi,astart,aend,ainc,frname); }
            catch (EOBException ex) { throw ex; }
            catch (FloatException flex) { throw flex; }
            catch (DataException dex) { throw dex; }
            catch (FileException fex) { throw fex; }
            break;
    }
    return;
}

void printCurvatureOnALatitude(Particle *p, Real ctheta, Real astart, Real aend,
                               Real ainc, std::string &frname)
{
    using std::cout;
    using std::endl;

    std::stringstream ssk1,ssk2,ssH,ssK;
    int itheta = int((ctheta * Rad2Deg) + 0.5);
    Real angle,k1,k2,H,K;
    
    // Open the input file
    ssk1 << frname << "-theta-" << itheta << "-k1.dat";
    ssk2 << frname << "-theta-" << itheta << "-k2.dat";
    ssH << frname << "-theta-" << itheta << "-H.dat";
    ssK << frname << "-theta-" << itheta << "-K.dat";

    std::string sk1 = ssk1.str();
    std::string sk2 = ssk2.str();
    std::string sH = ssH.str();
    std::string sK = ssK.str();

    std::ofstream osk1(sk1.c_str());
    if (!osk1) {
        throw FileException("diffgeom","printCurvatureOnALatitude",
                            sk1,"writing");
    }

    std::ofstream osk2(sk2.c_str());
    if (!osk2) {
        osk1.close();
        throw FileException("diffgeom","printCurvatureOnALatitude",
                            sk2,"writing");
    }

    std::ofstream osH(sH.c_str());
    if (!osH) {
        osk1.close();
        osk2.close();
        throw FileException("diffgeom","printCurvatureOnALatitude",
                            sH,"writing");
    }

    std::ofstream osK(sK.c_str());
    if (!osK) {
        osk1.close();
        osk2.close();
        osH.close();
        throw FileException("diffgeom","printCurvatureOnALatitude",
                            sK,"writing");
    }

    std::vector<Real> pd1,pd2;
    Real ncomp1,ncomp2;
    try {
        for (angle = astart; angle < aend; angle += ainc) {
            cout << "Point (" << Rad2Deg * ctheta << "," << Rad2Deg * angle << "): ";
            cout << "r = " << p->getR(ctheta,angle);
            cout.flush();
 
            k1 = p->computePrincipalcurvature(0,ctheta,angle);
            cout << "  k1 = " << std::setprecision(12) << k1;
            cout.flush();
            k2 = p->computePrincipalcurvature(1,ctheta,angle);
            cout << "  k2 = " << std::setprecision(12) << k2;
            cout.flush();
            H = p->getH(ctheta,angle);
            cout << "  H = " << std::setprecision(12) << H;
            cout.flush();
            K = p->getK(ctheta,angle);
            cout << "  K = " << std::setprecision(12) << K << endl;
            cout.flush();

            osk1 << Rad2Deg * angle << " " << std::setprecision(12) << k1 << endl;
            osk2 << Rad2Deg * angle << " " << std::setprecision(12) << k2 << endl;
            osH << Rad2Deg * angle << " " << std::setprecision(12) << H << endl;
            osK << Rad2Deg * angle << " " << std::setprecision(12) << K << endl;
        }
        osk1.close();
        osk2.close();
        osH.close();
        osK.close();
    }
    catch (EOBException ex) {
        osk1.close();
        osk2.close();
        osH.close();
        osK.close();
        throw ex;
    }
    catch (FloatException flex) {
        osk1.close();
        osk2.close();
        osH.close();
        osK.close();
        throw flex;
    }
    return;
}

void printCurvatureOnALongitude(Particle *p, Real cphi, Real astart, Real aend, Real ainc, std::string &frname)
{
    using std::cout;
    using std::endl;

    std::stringstream ssk1,ssk2,ssH,ssK;
    int iphi = int((cphi * Rad2Deg) + 0.5);
    Real angle,k1,k2,H,K;
    

    // Open the input file
    ssk1 << frname << "-phi-" << iphi << "-k1.dat";
    ssk2 << frname << "-phi-" << iphi << "-k2.dat";
    ssH << frname << "-phi-" << iphi << "-H.dat";
    ssK << frname << "-phi-" << iphi << "-K.dat";

    std::string sk1 = ssk1.str();
    std::string sk2 = ssk2.str();
    std::string sH = ssH.str();
    std::string sK = ssK.str();

    std::ofstream osk1(sk1.c_str());
    if (!osk1) {
        throw FileException("diffgeom","printCurvatureOnALongitude",
                            sk1,"writing");
    }

    std::ofstream osk2(sk2.c_str());
    if (!osk2) {
        osk1.close();
        throw FileException("diffgeom","printCurvatureOnALongitude",
                            sk2,"writing");
    }

    std::ofstream osH(sH.c_str());
    if (!osH) {
        osk1.close();
        osk2.close();
        throw FileException("diffgeom","printCurvatureOnALongitude",
                            sH,"writing");
    }

    std::ofstream osK(sK.c_str());
    if (!osK) {
        osk1.close();
        osk2.close();
        osH.close();
        throw FileException("diffgeom","printCurvatureOnALongitude",
                            sK,"writing");
    }

    try {
        for (angle = astart; angle < aend; angle += ainc) {
            cout << "Point (" << Rad2Deg * angle << "," << Rad2Deg * cphi << "): ";
            cout << "r = " << p->getR(angle,cphi);
            cout.flush();
            k1 = p->computePrincipalcurvature(0,angle,cphi);
            cout << "  k1 = " << std::setprecision(12) << k1;
            cout.flush();
            k2 = p->computePrincipalcurvature(1,angle,cphi);
            cout << "  k2 = " << std::setprecision(12) << k2;
            cout.flush();
            H = p->getH(angle,cphi);
            cout << "  H = " << std::setprecision(12) << H;
            cout.flush();
            K = p->getK(angle,cphi);
            cout << "  K = " << std::setprecision(12) << K << endl;
            cout.flush();

            osk1 << angle << " " << std::setprecision(12) << k1 << endl;
            osk2 << angle << " " << std::setprecision(12) << k2 << endl;
            osH << angle << " " << std::setprecision(12) << H << endl;
            osK << angle << " " << std::setprecision(12) << K << endl;
        }
        osk1.close();
        osk2.close();
        osH.close();
        osK.close();
    }
    catch (EOBException ex) {
        osk1.close();
        osk2.close();
        osH.close();
        osK.close();
        throw ex;
    }
    catch (FloatException flex) {
        osk1.close();
        osk2.close();
        osH.close();
        osK.close();
        throw flex;
    }
    return;
}

int shapestats(int argc, char *argv[])
{
    using std::cout;
    using std::cin;
    using std::endl;

    bool parallel_mode;
    std::string shapesetfolder;
    std::vector<std::string> shapeset;
    std::string buff,fname;
    std::ifstream in;
    std::string geomfilename,pgfname;
    std::string shapestatsfilename;
    int maxnamelength;
    int Lanczos_exponent,dstart;
    bool dogauss;
    std::vector<Sphere> maxinscribedsphere,minenclosingsphere;
    std::vector<Real> wadell_sphericity;
    std::vector<std::vector<Real> > Itensor,Rgtensor;
    std::vector<Real> pItensor,pRgtensor;
    std::vector<Real> Isphere,Rgsphere;
    Real pIsphere,pRgsphere;
    std::vector<Real> ratio_of_spheres;
    std::vector<std::vector<Real> > roundness;
    std::vector<Real> rvalues;
    std::vector<Real> length;
    std::vector<Real> width;
    std::vector<Real> thickness;
    std::vector<Real> D2;
    std::vector<Real> D3;
    std::vector<Real> corey;
    std::vector<Real> aschenbrenner;
    std::vector<Real> mps;
    std::vector<Real> esf;
    std::vector<Real> hofmann;
    std::vector<Real> convexity;
    std::ofstream out;
    
    #ifdef PARALLEL
        MPI::Init(argc,argv);
        int nbProcs = MPI::COMM_WORLD.Get_size(); // number of processes
        int rank = MPI::COMM_WORLD.Get_rank();    // rank of local process
        parallel_mode = true;

        // change some process options
        // silently fix unaligned user access
        // prctl(PR_SET_UNALIGN, PR_UNALIGN_NOPRINT, 0, 0, 0);
        // emulate floating point operations accesses
        // prctl(PR_SET_FPEMU, PR_FPEMU_NOPRINT, 0, 0, 0);
    #else
        int nbProcs = 1;
        int rank = 0;       // Default value when no parallel processing
        parallel_mode = false;
    #endif
    
    int npts;
    cout << "Enter number of Gaussian points to use:  ";
    cin >> npts;

    #ifdef PARALLEL
    if (rank == 0) { 
        cout << "Provide name of file containing shape sets to process: ";
        cin >> fname;
        cout << fname << endl;
        in.open(fname.c_str());
        if (!in) {
            cout << "ERROR:  Could not open " << fname << " for reading" << endl;
            throw FILE_OPEN_ERROR;
        }
        in >> shapesetfolder;
        shapeset.clear();
        while (!in.eof()) {
            in >> buff;
            cout << "Read shape set " << buff << endl;
            if (!in.eof()) {
                shapeset.push_back(buff);
            }
        }
        in.close();
    }
    #else
        cout << "Provide rank: ";
        cin >> rank;
        cout << rank << endl;
        cout << "Provide shape set folder name: ";
        cin >> shapesetfolder;
        cout << shapesetfolder << endl;
        cout << "Provide shape set name: ";
        cin >> buff;
        shapeset.clear();
        shapeset.push_back(buff);
    #endif
    
    std::ofstream gout;
    std::ostringstream partial_gfile_name;
    int linenum,numperfile;

    #ifdef PARALLEL
        if (rank == 0) { 

            // We are in parallel mode, so we split the geom file into
            // smaller geom files, one for each processor.

            // Loop over the number of particles (rows in the geom file)

            for (register int ss = 0; ss < shapeset.size(); ss++) {
                geomfilename = shapesetfolder + "/" + shapeset[ss] + "/" + shapeset[ss] + "-geom.dat";
                cout << "Geom file " << ss << ": " << geomfilename << endl;
                in.open(geomfilename.c_str(),std::ifstream::in);
                Lineitem oneline;
                line.clear();
                maxnamelength = 0;

                // Read one line of the geom file, corresponding to one particle

                while (!in.eof()) {
                    in >> oneline.name;   
                    if (oneline.name.length() > maxnamelength) maxnamelength = oneline.name.length();
                    cout << "  Shape " << oneline.name << endl;
                    if (!in.eof()) {
                        in >> oneline.xlow;   
                        in >> oneline.xhi;   
                        in >> oneline.ylow;   
                        in >> oneline.yhi;   
                        in >> oneline.zlow;   
                        in >> oneline.zhi;   
                        in >> oneline.volume;   
                        in >> oneline.surfarea;   
                        in >> oneline.nsurfarea;   
                        in >> oneline.diam;   
                        in >> oneline.Itrace;   
                        in >> oneline.Nnn;   
                        in >> oneline.NGC;   
                        in >> oneline.length;   
                        in >> oneline.width;   
                        in >> oneline.thickness;   
                        in >> oneline.nlength;   
                        in >> oneline.nwidth;   
                        line.push_back(oneline);
                    }
                }           // Done reading this line
                in.close();

                // The number of rows from the master geom file to include in
                // each smaller geom file

                numperfile = line.size() / nbProcs;

                while (linenum < line.size()) {

                    // create a geom file for each process

                    for (register int pn = 0; pn < nbProcs; pn++) {
                        partial_gfile_name << shapesetfolder << "/" << shapeset[ss]
                                           << "/" << shapeset[ss] << "-geom-rank-"
                                           << pn << ".dat";
                        pgfname = partial_gfile_name.str();
                        gout.open(pgfname.c_str());
                        for (register int ii = 0; ii < numperfile; ii++) {
                            gout << line[linenum].name << " "    
                                 << line[linenum].xlow << " "    
                                 << line[linenum].xhi << " "    
                                 << line[linenum].ylow << " "    
                                 << line[linenum].yhi << " "    
                                 << line[linenum].zlow << " "    
                                 << line[linenum].zhi << " "    
                                 << line[linenum].volume << " "    
                                 << line[linenum].surfarea << " "    
                                 << line[linenum].nsurfarea << " "    
                                 << line[linenum].diam << " "    
                                 << line[linenum].Itrace << " "    
                                 << line[linenum].Nnn << " "    
                                 << line[linenum].NGC << " "    
                                 << line[linenum].length << " "    
                                 << line[linenum].width << " "    
                                 << line[linenum].thickness << " "    
                                 << line[linenum].nlength << " "    
                                 << line[linenum].nwidth << endl;   
                                 linenum++;
                        }
                    }
                }
            }
        }

    #endif
    
    cout << "OK, done.  Shape set size is " << shapeset.size() << endl;

    // Master loop over all shape sets (i.e., particle directories) specified in the input.
    for (register int ss = 0; ss < shapeset.size(); ss++) {
        partial_gfile_name.str("");
        // partial_gfile_name << shapesetfolder << "/" << shapeset[ss]
        //                 << "/" << shapeset[ss] << "-geom-rank-"
        //                 << rank << ".dat";
        partial_gfile_name << shapesetfolder << "/" << shapeset[ss]
                           << "/" << shapeset[ss] << "-geom-rank-" << rank << ".dat";
        geomfilename = partial_gfile_name.str();
        cout << "Rank: " << rank << ", Geom file " << ss << ": " << geomfilename << endl;
        cout << "Geom file " << ss << ": " << geomfilename << endl;
        in.open(geomfilename.c_str(),std::ifstream::in);

        Lineitem oneline;
        std::vector<Lineitem> line;
        line.clear();
        maxnamelength = 0;
        while (!in.eof()) {
            in >> oneline.name;   
            if (oneline.name.length() > maxnamelength) maxnamelength = oneline.name.length();
            cout << "  Shape " << oneline.name << endl;
            if (!in.eof()) {
                in >> oneline.xlow;   
                in >> oneline.xhi;   
                in >> oneline.ylow;   
                in >> oneline.yhi;   
                in >> oneline.zlow;   
                in >> oneline.zhi;   
                in >> oneline.volume;   
                in >> oneline.surfarea;   
                in >> oneline.nsurfarea;   
                in >> oneline.diam;   
                in >> oneline.Itrace;   
                in >> oneline.Nnn;   
                in >> oneline.NGC;   
                in >> oneline.length;   
                in >> oneline.width;   
                in >> oneline.thickness;   
                in >> oneline.nlength;   
                in >> oneline.nwidth;   
                line.push_back(oneline);
            }
        }
        in.close();
    
        // The following vectors will hold the respective property for each
        // particle analyzed.  Clear them out here.

        maxinscribedsphere.clear();
        minenclosingsphere.clear();
        wadell_sphericity.clear();
        Itensor.clear();
        Rgtensor.clear();
        ratio_of_spheres.clear();
        roundness.clear();
        rvalues.clear();
        length.clear();
        width.clear();
        thickness.clear();
        D2.clear();
        D3.clear();
        corey.clear();
        aschenbrenner.clear();
        mps.clear();
        esf.clear();
        hofmann.clear();
        convexity.clear();
    
        Particle *p;
        std::string anmfilename;
        rvalues.resize(ALLROUNDNESSES,0.0);
        roundness.resize(line.size(),rvalues);
        std::ostringstream ssfname;
        ssfname << shapesetfolder << "/" << shapeset[ss] << "/" << shapeset[ss] 
               << "-rank-" << rank << "-shapestats.dat";
        // ssfname << shapeset[ss] << "-rank-" << rank << "-shapestats.dat";
        shapestatsfilename = ssfname.str();
        cout << "Rank: " << rank << ", Shape stats file " << ss << ": " << shapestatsfilename << endl;

        // Prepare the master output file by writing its header

        #ifdef PARALLEL
        if (rank == 0) {
            out.open(shapestatsfilename.c_str(),std::ios_base::app);
            buff = "Name";
            out << std::setw(maxnamelength+4) << std::left << buff;
            buff = "IS_Rad";
            out << std::setw(12) << std::right << buff;
            buff = "IS_x";
            out << std::setw(12) << std::right << buff;
            buff = "IS_y";
            out << std::setw(12) << std::right << buff;
            buff = "IS_z";
            out << std::setw(12) << std::right << buff;
            buff = "ES_Rad";
            out << std::setw(12) << std::right << buff;
            buff = "ES_x";
            out << std::setw(12) << std::right << buff;
            buff = "ES_y";
            out << std::setw(12) << std::right << buff;
            buff = "ES_z";
            out << std::setw(12) << std::right << buff;
            buff = "Acoeff";
            out << std::setw(12) << std::right << buff;
            buff = "Sphericity";
            out << std::setw(12) << std::right << buff;
            buff = "BSR";
            out << std::setw(12) << std::right << buff;
            buff = "Round DP";
            out << std::setw(12) << std::right << buff;
            buff = "Round MaxC";
            out << std::setw(12) << std::right << buff;
            buff = "Length";
            out << std::setw(12) << std::right << buff;
            buff = "Width";
            out << std::setw(12) << std::right << buff;
            buff = "Thick";
            out << std::setw(12) << std::right << buff;
            buff = "D2";
            out << std::setw(12) << std::right << buff;
            buff = "D3";
            out << std::setw(12) << std::right << buff;
            buff = "Corey SF";
            out << std::setw(12) << std::right << buff;
            buff = "ASF";
            out << std::setw(12) << std::right << buff;
            buff = "MPS";
            out << std::setw(12) << std::right << buff;
            buff = "ESF";
            out << std::setw(12) << rstd::ight << buff;
            buff = "Hofmann";
            out << std::setw(12) << std::right << buff;
            buff = "Volume";
            out << std::setw(12) << std::right << buff;
            buff = "Area";
            out << std::setw(12) << std::right << buff;
            buff = "I11";
            out << std::setw(12) << std::right << buff << endl;
            buff = "I22";
            out << std::setw(12) << std::right << buff << endl;
            buff = "I33";
            out << std::setw(12) << std::right << buff << endl;
            buff = "Isphere";
            out << std::setw(12) << std::right << buff << endl;
            buff = "Rg11";
            out << std::setw(12) << std::right << buff << endl;
            buff = "Rg22";
            out << std::setw(12) << std::right << buff << endl;
            buff = "Rg33";
            out << std::setw(12) << std::right << buff << endl;
            buff = "Rgsphere";
            out << std::setw(12) << std::right << buff << endl;
            out.close();
        }
        #else
            out.open(shapestatsfilename.c_str(),std::ios_base::app);
            buff = "Name";
            out << std::setw(maxnamelength+4) << std::left << buff;
            buff = "I11";
            out << std::setw(12) << std::right << buff;
            buff = "I22";
            out << std::setw(12) << std::right << buff;
            buff = "I33";
            out << std::setw(12) << std::right << buff;
            buff = "Isphere";
            out << std::setw(12) << std::right << buff;
            buff = "Rg11";
            out << std::setw(12) << std::right << buff;
            buff = "Rg22";
            out << std::setw(12) << std::right << buff;
            buff = "Rg33";
            out << std::setw(12) << std::right << buff;
            buff = "Rgsphere";
            out << std::setw(12) << std::right << buff << endl;
            out.close();
        #endif

        // Master loop over all particles within this shape set

        for (register int pnum = 0; pnum < line.size(); pnum++) {
    
            anmfilename = shapesetfolder + "/" + shapeset[ss] + "/" + line[pnum].name;
            cout << endl << "Now reading SH coeffs for shape " << pnum << ": " << anmfilename;
            cout.flush();
            bool verbose = false;
            p = new Star(anmfilename,verbose);
            cout << endl << " Done!" << endl;
            cout.flush();
    
            Lanczos_exponent = 1;
            dstart = 10;
            cout << "De-ringing the particle now... Exponent = " << Lanczos_exponent
                 << ", starting at degree " << dstart << endl;
            cout.flush();
            p->doLanczos(1,10);
            cout << "Done!" << endl;
            cout.flush();
    
//          p->computeDimensions(true);
    
//          npts = 50;
//          cout << "Setting " << npts << " points in each direction for bounding spheres";
//          p->setSurface(npts,p->getNmax(),false);
//          cout << "Computing bounding spheres ... ";
//          cout.flush();
//          p->computeBoundingspheres();
//          cout << "Done computing bounding spheres." << endl;
//          cout.flush();
    
            dogauss = true;
            npts = 120;
            cout << "Setting surface points (" << npts << " in each direction) ... " << endl;
            cout << "Nmax = " << line[pnum].Nnn << endl;
            cout.flush();
            p->setSurface(npts,line[pnum].Nnn,dogauss);
            cout << "Done!" << endl;
            cout.flush();
    
            // Calculate all shape factors for this particle

//          cout << "Compiling shape factors for shape " << pnum << endl;

//          maxinscribedsphere.push_back(p->getMaxinscribedsphere());
//          minenclosingsphere.push_back(p->getMinenclosingsphere());

//          cout << "    Max inscribed sphere radius " << (p->getMaxinscribedsphere()).getRadius() << " at ";
//          cout << "(" << (p->getMaxinscribedsphere()).getXr() << "," << (p->getMaxinscribedsphere()).getYr();
//          cout << "," << (p->getMaxinscribedsphere()).getZr() << ")" << endl;
//          cout << "    Min enclosing sphere radius " << (p->getMinenclosingsphere()).getRadius() << " at ";
//          cout << "(" << (p->getMinenclosingsphere()).getXr() << "," << (p->getMinenclosingsphere()).getYr();
//          cout << "," << (p->getMinenclosingsphere()).getZr() << ")" << endl;

//          wadell_sphericity.push_back(p->getWadellsphericity());
//          cout << "    Wadell sphericity = " << p->getWadellsphericity() << endl;
//
//          ratio_of_spheres.push_back(p->getBoundingsphereratio());
//          cout << "    Ratio of spheres = " << p->getBoundingsphereratio() << endl;

//          p->computeRoundness(ALLROUNDNESSES);
//          cout << "    Dot product roundness = " << p->getRoundness(DOTPRODUCT) << endl;
//          cout << "    Wadell roundness = " << p->getRoundness(WADELL) << endl;
//          roundness[pnum][DOTPRODUCT] = p->getRoundness(DOTPRODUCT);
//          roundness[pnum][WADELL] = p->getRoundness(WADELL);

//          cout << "L = " << p->getDim(0) << endl;
//          cout << "W = " << p->getDim(1) << endl;
//          cout << "T = " << p->getDim(2) << endl;
//          cout << "D2 = " << p->getTriaxialdim(1) << endl;
//          cout << "D3 = " << p->getTriaxialdim(2) << endl;

//          length.push_back(p->getDim(0));
//          width.push_back(p->getDim(1));
//          thickness.push_back(p->getDim(2));
//          D2.push_back(p->getTriaxialdim(1));
//          D3.push_back(p->getTriaxialdim(2));
//
//          corey.push_back(p->getCSF());
//          cout << "    Corey shape factor = " << p->getCSF() << endl;
//          aschenbrenner.push_back(p->getASF());
//          cout << "    Aschenbrenner shape factor = " << p->getASF() << endl;
//          mps.push_back(p->getMPS());
//          cout << "    Max Projected Shape factor = " << p->getMPS() << endl;
//          esf.push_back(p->getESF());
//          cout << "    ESF factor = " << p->getESF() << endl;
//          hofmann.push_back(p->getShapeentropy());
//          cout << "    Hofmann shape entropy = " << p->getShapeentropy() << endl;
//          convexity.push_back(p->getConvexity());
            /*
            p->computeVolume();
            p->computeArea();
            area.push_back(p->getArea());
            volume.push_back(p->getVolume());
            */
//          cout << "    Convexity = " << convexity[pnum] << endl;
    
            // Get diagonalized moment-of-inertia tensor and gyration tensor
            p->getI(pItensor,pRgtensor,pIsphere,pRgsphere);
            Itensor.push_back(pItensor);
            Rgtensor.push_back(pRgtensor);
            Isphere.push_back(pIsphere);
            Rgsphere.push_back(pRgsphere);

            // Output results for this particle to the master output file for this shape set

            out.open(shapestatsfilename.c_str(),std::ios_base::app);
            out << std::setw(maxnamelength+4) << std::left << line[pnum].name;
//          out << std::setprecision(3) << std::setw(12) << std::right << maxinscribedsphere[pnum].getRadius(); 
//          out << std::setprecision(3) << std::setw(12) << std::right << maxinscribedsphere[pnum].getXr(); 
//          out << std::setprecision(3) << std::setw(12) << std::right << maxinscribedsphere[pnum].getYr(); 
//          out << std::setprecision(3) << std::setw(12) << std::right << maxinscribedsphere[pnum].getZr(); 
//          out << std::setprecision(3) << std::setw(12) << std::right << minenclosingsphere[pnum].getRadius(); 
//          out << std::setprecision(3) << std::setw(12) << std::right << minenclosingsphere[pnum].getXr(); 
//          out << std::setprecision(3) << std::setw(12) << std::right << minenclosingsphere[pnum].getYr(); 
//          out << std::setprecision(3) << std::setw(12) << std::right << minenclosingsphere[pnum].getZr(); 
//          out << std::setprecision(3) << std::setw(12) << std::right << line[pnum].surfarea/(pow(line[pnum].volume,(2.0/3.0))); 
//          out << std::setprecision(3) << std::setw(12) << std::right << wadell_sphericity[pnum]; 
//          out << std::setprecision(3) << std::setw(12) << std::right << ratio_of_spheres[pnum]; 
//          out << std::setprecision(3) << std::setw(12) << std::right << roundness[pnum][DOTPRODUCT];
//          out << std::setprecision(3) << std::setw(12) << std::right << roundness[pnum][WADELL];
//          out << std::setprecision(3) << std::setw(12) << std::right << length[pnum];
//          out << std::setprecision(3) << std::setw(12) << std::right << width[pnum];
//          out << std::setprecision(3) << std::setw(12) << std::right << thickness[pnum];
//          out << std::setprecision(3) << std::setw(12) << std::right << D2[pnum];
//          out << std::setprecision(3) << std::setw(12) << std::right << D3[pnum];
//          out << std::setprecision(3) << std::setw(12) << std::right << corey[pnum];
//          out << std::setprecision(3) << std::setw(12) << std::right << aschenbrenner[pnum];
//          out << std::setprecision(3) << std::setw(12) << std::right << mps[pnum];
//          out << std::setprecision(3) << std::setw(12) << std::right << esf[pnum];
//          out << std::setprecision(3) << std::setw(12) << std::right << hofmann[pnum];
//          out << std::setprecision(7) << std::setw(12) << std::right << line[pnum].volume;
//          out << std::setprecision(7) << std::setw(12) << std::right << line[pnum].surfarea;
//          out << std::setprecision(7) << std::setw(12) << std::right << convexity[pnum];
            out << std::setprecision(7) << std::setw(12) << std::right << Itensor[pnum][0];
            out << std::setprecision(7) << std::setw(12) << std::right << Itensor[pnum][1];
            out << std::setprecision(7) << std::setw(12) << std::right << Itensor[pnum][2];
            out << std::setprecision(7) << std::setw(12) << std::right << Isphere[pnum];
            out << std::setprecision(7) << std::setw(12) << std::right << Rgtensor[pnum][0];
            out << std::setprecision(7) << std::setw(12) << std::right << Rgtensor[pnum][1];
            out << std::setprecision(7) << std::setw(12) << std::right << Rgtensor[pnum][2];
            out << std::setprecision(7) << std::setw(12) << std::right << Rgsphere[pnum] << endl;
            out.close();
    
            cout << "Deleting the particle memory" << endl;
            cout.flush();
            delete p;
            cout << "  Done!" << endl;
            cout.flush();

        } // End of loop over all particles within a shape set

    }     // End of master loop over all shape sets
    
    // Now the program is basically done.  Do final output and clean up

    #ifdef PARALLEL
        MPI::Finalize();

        std::ostringstream cmd,ssfname2,ssfname3;
        std::string ssname;
        for (register int ss = 0; ss < shapeset.size(); ss++) {
            cmd.str("");  // Clear the cmd string
            ssfname.str("");
            ssfname2.str("");
            ssfname3.str("");
            ssname = shapeset[ss];
            ssfname << shapesetfolder << "/" << ssname << "/" << ssname
                    << "-rank-0-shapestats.dat";
            ssfname2 << shapesetfolder << "/" << ssname << "/" << ssname
                        << "-rank-1-shapestats.dat";
            ssfname3 << shapesetfolder << "/" << ssname << "/" << ssname 
                        << "." << 100 << ".shapestats.dat";
            cmd << "cat " << ssfname.str() << " " << ssfname2.str() << " > " << ssfname3.str();
            std::system((cmd.str()).c_str());
            for (register int i = 1; i < nbProcs; ss++) {
                cmd.str("");  // Clear the cmd string
                ssfname.str("");
                ssfname2.str("");
                ssfname3.str("");
                ssfname << shapesetfolder << "/" << ssname << "/" << ssname
                        << "." << i + 99 << ".shapestats.dat";
                ssfname2 << shapesetfolder << "/" << ssname << "/" << ssname
                        << "-rank-" << i << "-shapestats.dat";
                ssfname3 << shapesetfolder << "/" << ssname << "/" << ssname
                        << "." << i + 100 << ".shapestats.dat";
                cmd << "cat " << ssfname.str() << " " << ssfname2.str() << " > " << ssfname3.str();
                std::system((cmd.str()).c_str());
            }
            cmd.str("");
            ssfname.str("");
            ssfname << shapesetfolder << "/" << ssname << "/" << ssname
                    << ".shapestats.dat";
            cmd << "mv " << ssfname3.str() << ssfname.str();
            std::system((cmd.str()).c_str());
        }
    #endif
}
