/**
@file
@brief Functions for manipulating images using CGAL library
*/
#ifndef IMAGEMANIPULATIONH
#define IMAGEMANIPULATIONH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "global.h"
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

/**
@typedef
@brief Shorthand for a CGAL default surface triangulation.
*/
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

/**
@typedef
@brief Collection of geometric traits for a surface triangulation.
*/
typedef Tr::Geom_traits GT;

/**
@typedef
@brief Shorthand for a CGAL gray level image object.
*/
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;

/**
@brief Converts a binary digital image to gray level image.

The following function takes a file filled with voxel information and
returns a [[CGAL]] gray level image for further processing.  The format
file filled with the voxel information must have the following heading:

\verbatim
Xsize: 
Ysize: 
Zsize:
VoxelScale:
\endverbatim

Following this heading, one voxel integer id per pixel is required.  Non-zero numbers are all
treated as 1.  The voxels are listed one <i>z</i> plane at a time, starting with <i>z</i> = 0.
Within each plane, the voxels are listed one <i>y</i> column at a time, starting with <i>y</i> = 0,
and within each <i>y</i> column, the voxels are listed one <i>x</i> row at a time, starting
with <i>x</i> = 0.

@param filename is the name of the binary image file
@param iso_value is the isosurface value to use for the conversion (WHY WHAT IS THIS?)
@return the gray level image object
*/
Gray_level_image binary2Gray(const std::string filename, Real iso_value);

/**
@brief Gets a voxel image from an input file.

@param &fname is the name of the input file holding the image
@param &xBB will hold the <i>x</i> dimension of the image bounding box
@param &yBB will hold the <i>y</i> dimension of the image bounding box
@param &zBB will hold the <i>z</i> dimension of the image bounding box
@param &vscale is (WHY WHAT IS THIS?)
@return a 3D vector array of the voxels in the image
*/
std::vector<std::vector<std::vector<Real> > > getVoxelImage(std::string &fname, int &xBB,
                                             int &yBB, int &zBB, Real &vscale);
/**
@brief Finds the centroid of the image.

@param &voxel is a 3D voxel array (nested STL vectors)
@return the three coordinates of the centroid stored in an STL vector
*/
std::vector<Real> getImageCentroid(std::vector<std::vector<std::vector<Real> > > &voxel);

#endif
