#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([tsquare], [1.0], [jeffrey.bullard@nist.gov])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/diffgeom.cc])
AC_CONFIG_HEADERS([config.h])

# Tell autoconf where to find macros
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Check for Doxygen
AC_CHECK_PROGS([DOXYGEN], [doxygen], [AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])])

AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([docs/Doxyfile])])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

AC_PROG_RANLIB

# Checks for libraries
AC_CHECK_LIB([m],[atan2])
AX_BOOST_BASE([1.58], [], [AC_MSG_ERROR([Sorry, boost not found])])
#AC_CHECK_LIB(gmp, __gmpz_init2, ,[AC_MSG_ERROR("Linking GMP failed.")])
#AC_CHECK_LIB(mpfr, mpfr__init2, ,[AC_MSG_ERROR("Linking MPFR failed.")])
AX_LIB_CGAL_CORE([], [AC_MSG_ERROR([Sorry, cgal not found])])

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt exp cos sin acos atan atan2])

AC_CONFIG_FILES([Makefile
                 common/Makefile
                 src/Makefile]
		docs/Makefile)
AC_OUTPUT
