/**
@file
@brief Global parameters used by all classes.
*/

#ifndef GLOBALH
#define GLOBALH

// #define DEBUG
// #define VERBOSE
// #define PARALLEL

/**
@typedef
@brief Defines the `Real` data type as `double`.

This makes it easy to quickly switch back and forth between
`double` and `float` if necessary.
*/
typedef double Real;

const int RETURN_NORMAL = 0;              // Normal return condition
const int CONTAINER_EOB = 1;              // Container element out of bounds
const int FILE_OPEN_ERROR = 2;            // Could not open a requested file
const int BAD_SURFACE_ELEMENT = 3;        // Surface element index out of range
const int PREMATURE_EOF = 4;              // Premature end of file
const int INVALID_INPUT = 5;              // Bad data
const int BAD_DIMENSION = 6;              // Bad dimension

const int MAXSTRING = 128;

// Saturation value for colors
const Real COLORSATVAL = 255.0;

// Value of pi
const Real Pi = 3.141592653589793;
const Real Deg2Rad = Pi/180.0;
const Real Rad2Deg = 180.0/Pi;

const int DOTPRODUCT = 0;
const int WADELL = 1;
const int ALLROUNDNESSES = 2;

const int RDIM = 0;
const int XDIM = 0;
const int YDIM = 1;
const int ZDIM = 2;

#include <stdexcept>
#include "Exceptions.h"
#endif
