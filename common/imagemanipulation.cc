/**
@file
@brief Collection of method definitions for image manipulation using CGAL library.
*/
#include "imagemanipulation.h"

Gray_level_image binary2Gray(const std::string filename, Real iso_value)
{
    using std::cout;
    using std::endl;

    std::string buff;
    int id,xsize,ysize,zsize,voxelsize;
    Real voxelscale;
    std::vector<int> voxelValues;
    voxelValues.clear();

    char *data = NULL;

    try {
        std::ifstream in(filename.c_str());
        if (!in) {
            throw tsquare::FileException("imagemanipulation","binary2Gray",filename,
                                "Could not open for reading");
        }

        // read binary header one line at a time
        cout << "...Reading x size" << endl;
        in >> buff >> xsize;
        if (buff != "Xsize:") {
            throw tsquare::FileException("imagemanipulation","binary2Gray",filename,
                                "Header format is incorrect.");
        }
        cout << "...Reading y size" << endl;
        in >> buff >> ysize;
        if (buff != "Ysize:") {
            throw tsquare::FileException("imagemanipulation","binary2Gray",filename,
                                "Header format is incorrect.");
        }
        cout << "...Reading z size" << endl;
        in >> buff >> zsize;
        if (buff != "Zsize:") {
            throw tsquare::FileException("imagemanipulation","binary2Gray",filename,
                                "Header format is incorrect.");
        }
        cout << "...Reading voxel scale" << endl;
        in >> buff >> voxelscale;
        if (buff != "VoxelScale:") {
            throw tsquare::FileException("imagemanipulation","binary2Gray",filename,
                                "Header format is incorrect.");
        }

        voxelValues.resize(xsize*ysize*zsize,0);
     
        id = 0;
        cout << "...Reading all the voxels" << endl;
        for (unsigned int k = 0; k < zsize; k++) {
            for (unsigned int j = 0; j < ysize; j++) {
                for (unsigned int i = 0; i < xsize; i++) {
                    if (in.eof()) {
                        throw tsquare::FileException("imagemanipulation","binary2Gray",filename,
                                "Premature end of file.");
                    }
                    in >> voxelValues[id];
                    id++;
                }
            }
        }
      
        in.close();

        cout << "...Attempting to create INRIMAGE file" << endl;
        cout.flush();

        //output .inr header
        std::ostringstream header;
        header << "#INRIMAGE-4#{\nXDIM=" << xsize << "\nYDIM=" << ysize
               << "\nZDIM=" << zsize << "\nVDIM=1\nTYPE=unsigned fixed\nPIXSIZE=8 bits\nCPU=decm\nVX="
               << voxelscale << "\nVY=" << voxelscale << "\nVZ=" << voxelscale << "\n";
        std::string headerString = header.str();
        std::string headerEND = "##}\n";

        int hlen = 256 - headerString.length() - headerEND.length();
        for (unsigned int i = 0; i < hlen; i++) {
            headerString += '\n';
        }
        headerString += headerEND;
        
        cout << "...Header completed, storing raw data now." << endl;
        cout.flush();

        // Next, we store the raw data in an array of chars

        char *data = new char[xsize*ysize*zsize];

        for (id = 0; id < (xsize*ysize*zsize); id++) {
            data[id] = (char)(voxelValues[id]);
        }

        // Now we must create the inrimage file for the image

        cout << "...Raw data stored, opening output file now" << endl;
        cout.flush();

        std::ofstream fp("inrimage.inr",std::ios::out | std::ios::binary);
        if (!fp) {
            throw tsquare::FileException("imagemanipulation","binary2Gray","inrimage.inr",
                                "Could not open inrimage for writing");
        }

        cout << "...Output file open, writing data now." << endl;
        cout.flush();

        fp.write(headerString.c_str(),headerString.size());
        fp.write(data,xsize*ysize*zsize);
        fp.close();

        cout << "...Data written.  Freeing allocated memory now." << endl;
        cout.flush();

        // Clean up the allocated memory
        delete [] data;

        // Now create the CGAL gray level image from the inrimage

        cout << "...Memory freed.  Attempting to construct gray level image now." << endl;
        cout.flush();

        Gray_level_image image("inrimage.inr",iso_value);
        return (image);

    }
    catch (tsquare::FileException fex) { 
        fex.printException();
        if (data != NULL) delete [] data;
        exit(1);
    }
}

std::vector<std::vector<std::vector<Real> > > getVoxelImage(std::string &fname, int &xBB,
                                                       int &yBB, int &zBB, Real &vscale)
{
    using std::cout;
    using std::endl;

    std::string buff;
    std::vector<std::vector<std::vector<Real> > > voxel;
    Real inputVal;

    try {
        cout << "Opening " << fname << endl;
        cout.flush();

        std::ifstream in(fname.c_str());
        if (!in) {
            throw tsquare::FileException("Particle","getVoxelImage",fname,
                            "Could not open for reading");
        }

        // First read the X, Y, and Z dimensions of the bounding box

        in >> buff >> xBB;
        cout << "...X size = " << xBB << endl;
        in >> buff >> yBB;
        cout << "...Y size = " << yBB << endl;
        in >> buff >> zBB;
        cout << "...Z size = " << zBB << endl;
        in >> buff >> vscale;
        cout << "...Voxel Scale = " << vscale << endl;
        cout.flush();

        voxel.clear();
        std::vector<Real> zdim;
        zdim.clear();
        zdim.resize(zBB,0.0);
        std::vector<std::vector<Real> > ydim;
        ydim.clear();
        ydim.resize(yBB,zdim);
        voxel.clear();
        voxel.resize(xBB,ydim);

        // Next, read in and store the Cartesian coordinates of the reference

        cout << "...Reading rest of the data" << endl;
        cout.flush();

        for (unsigned int k = 0; k < zBB; k++) {
          for (unsigned int j = 0; j < yBB; j++) {
            for (unsigned int i = 0; i < xBB; i++) {
              in >> inputVal;
              voxel[i][j][k] = inputVal;
            }
          }
        }

        in.close();

        cout << "Closed file" << endl;
        cout.flush();

    }
    catch (tsquare::FileException fex) {
        fex.printException();
        exit(1);
    }

    return voxel;
}

std::vector<Real> getImageCentroid(std::vector<std::vector<std::vector<Real> > > &voxel)
{
    using std::cout;
    using std::endl;

    std::vector<Real> centroid;
    centroid.clear();
    centroid.resize(3,0.0);
    Real voxsum = 0.0;

    int xsize = voxel.size();
    int ysize = voxel[0].size();
    int zsize = voxel[0][0].size();

    cout << "Getting Image Centroid:  xsize = " << xsize
         << ", ysize = " << ysize << ", zsize = " << zsize << endl;
    cout.flush();

    for (unsigned int i = 0; i < xsize; i++) {
      for (unsigned int j = 0; j < ysize; j++) {
        for (unsigned int k = 0; k < zsize; k++) {
            centroid[XDIM] += voxel[i][j][k] * (Real)(i);          
            centroid[YDIM] += voxel[i][j][k] * (Real)(j);          
            centroid[ZDIM] += voxel[i][j][k] * (Real)(k);          
            voxsum += voxel[i][j][k];
        }
      }
    }

    centroid[XDIM] *= (1.0/voxsum);
    centroid[YDIM] *= (1.0/voxsum);
    centroid[ZDIM] *= (1.0/voxsum);

    return centroid;
}
