# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                              Sundance
#                 Copyright (2005) Sandia Corporation
# 
# Copyright (year first published) Sandia Corporation.  Under the terms 
# of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
# retains certain rights in this software.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#                                                                                 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA                                                                                
# Questions? Contact Kevin Long (krlong@sandia.gov), 
# Sandia National Laboratories, Livermore, California, USA
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT(Sundance, 2.0, krlong@sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running Sundance Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(/src-core/UserExprs/SundanceExpr.hpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/Support/SundanceConfig.h 
# from src/Support/SundanceConfig.h.in

AM_CONFIG_HEADER(src-solvers/TSF_config.h:src-solvers/TSF_config.h.in)
AM_CONFIG_HEADER(src-utils/SundanceConfig.h:src-utils/SundanceConfig.h.in)


# 
#AC_DISABLE_SHARED
#AC_PROG_LIBTOOL


# Allow users to specify their own "install" command.  If none is specified,
# the default is install-sh found in the config subdirectory.
									
AC_ARG_WITH(install,
[AC_HELP_STRING([--with-install=INSTALL_PROGRAM],
[Use the installation program INSTALL_PROGRAM rather the default that is provided.  For example --with-install="/path/install -p"])],
[
INSTALL=$withval
INSTALL_PROGRAM=$withval
INSTALL_SCRIPT=$withval
INSTALL_DATA="$withval -m 644"
],)
     
# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.

AM_MAINTAINER_MODE

# Define $build, $host, $target, etc

AC_CANONICAL_TARGET

# Use automake

AM_INIT_AUTOMAKE(1.8.2)

# Specify required version of autoconf.

AC_PREREQ(2.52)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for compilers
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB


# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package 
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# Check for Chaco
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_CHACO


# ------------------------------------------------------------------------
# Check for Trilinos
# ------------------------------------------------------------------------


AC_ARG_WITH(trilinos,
AC_HELP_STRING([--with-trilinos],[root directory of Trilinos installation]),
[
if test "X$withval" = "Xyes"; then
  TRILINOS_DIR=`pwd`/../..
else
  TRILINOS_DIR=${withval}
fi
AC_MSG_CHECKING(Trilinos root directory)
AC_MSG_RESULT([${TRILINOS_DIR}])
])

# -- Set path to Trilinos data


AC_ARG_WITH(trilinosdatadir,
AC_HELP_STRING([--with-trilinosdatadir],[location of Trilinos data]),
[
  TRILINOS_DATA_DIR=${withval}
  AC_MSG_CHECKING(Trilinos data directory)
  AC_MSG_RESULT([${TRILINOS_DATA_DIR}])
])




# -- Set Trilinos version

AC_ARG_WITH(trilinos_version,
AC_HELP_STRING([--with-trilinos-version],[Which version of Trilinos to link against]),
ac_cv_use_trilinos_version=$withval, ac_cv_use_trilinos_version=dev)

AC_MSG_CHECKING(which version of Trilinos to link against)

if test "X$ac_cv_use_trilinos_version" = "Xdev"; then
AC_DEFINE(TRILINOS_DEV,,[define if we are using the development version of Trilinos])
trilinosVersion=dev
AC_MSG_RESULT($ac_cv_use_trilinos_version)
else if test "X$ac_cv_use_trilinos_version" = "X7"; then
AC_DEFINE(TRILINOS_7,,[define if we are using Trilinos 7.0.x])
AC_MSG_RESULT($ac_cv_use_trilinos_version)
else 
AC_DEFINE(TRILINOS_8,,[define if we are using Trilinos 8.0.x])
if test "X$ac_cv_use_trilinos_version" = "X8"; then
AC_MSG_RESULT($ac_cv_use_trilinos_version)
else
trilinosVersion=dev
AC_MSG_RESULT(unknown or no version specified.  Defaulting to Trilinos Devel)
fi
fi
fi

AM_CONDITIONAL(TRILINOS_IS_EXTERNAL, test "X$trilinosVersion" != "Xdev" )



# -- Set path to Trilinos source

AC_ARG_WITH(trilinossrcdir,
AC_HELP_STRING([--with-trilinossrcdir],[location of Trilinos source code]),
[
  TRILINOS_SRC_DIR=${withval}
  AC_MSG_CHECKING(Trilinos src directory)
  AC_MSG_RESULT([${TRILINOS_SRC_DIR}])
])

echo "trilinos src dir=(" $TRILINOS_SRC_DIR ")"
echo "trilinos version=(" $trilinosVersion ")"



# ------------------------------------------------------------------------
# Check for FIAT
# ------------------------------------------------------------------------

AC_ARG_WITH(fiat,
[AC_HELP_STRING([--with-fiat=FIAT_BUILD_DIR],[use FIAT build directory (enables FIAT)])],
[
FIAT_STR=""
FIAT_BUILD_DIR=${withval}
AC_DEFINE(HAVE_FIAT,,[define if we want to use FIAT])
AC_MSG_CHECKING(FIAT directory)
AC_MSG_RESULT([${FIAT_BUILD_DIR}])
],
[
FIAT_STR="#"
AC_MSG_CHECKING(FIAT directory)
AC_MSG_RESULT([no])
]
)


# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Check for Sundance tests
# ------------------------------------------------------------------------

AC_ARG_ENABLE(sundance-tests, 
AC_HELP_STRING([--enable-sundance-tests],[whether to build Sundance tests]),
[case "${enableval}" in
  yes) sundanceTests=yes ;;
  no) sundanceTests=no;;
  *) AC_MSG_ERROR([bad value for ${enableval} for --enable-sundance-tests]);;
esac],[sundanceTests=yes])

AC_MSG_CHECKING(whether to enable Sundance tests)
AC_MSG_RESULT([$sundanceTests])
AM_CONDITIONAL(BUILD_TESTS, test "X$sundanceTests" = "Xyes")




# ------------------------------------------------------------------------
# Check for Sundance examples
# ------------------------------------------------------------------------

AC_ARG_ENABLE(sundance-examples, 
AC_HELP_STRING([--enable-sundance-examples],[whether to build Sundance examples]),
[case "${enableval}" in
  yes) sundanceExamples=yes ;;
  no) sundanceExamples=no;;
  *) AC_MSG_ERROR([bad value for ${enableval} for --enable-sundance-examples]);;
esac],[sundanceExamples=yes])

AC_MSG_CHECKING(whether to enable Sundance examples)
AC_MSG_RESULT([$sundanceExamples])
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$sundanceExamples" = "Xyes")

# These will never be defined, but are needed for building the tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")


# ------------------------------------------------------------------------
# Check for whether to enable moocho in Sundance
# ------------------------------------------------------------------------

AC_ARG_ENABLE(moocho, 
AC_HELP_STRING([--enable-moocho],[whether to build moocho]),
[case "${enableval}" in
  yes) sundanceMoocho=yes ;;
  no) sundanceMoocho=no;;
  *) AC_MSG_ERROR([bad value for ${enableval} for --enable-moocho]);;
esac],[sundanceMoocho=no])

AC_MSG_CHECKING(whether to enable moocho)
AC_MSG_RESULT([$sundanceMoocho])


if test "X$sundanceMoocho" = "Xyes"  ; then
  MOOCHO_STR=""
  AC_DEFINE(HAVE_ENABLED_MOOCHO,,[define if we want to use MOOCHO])
else
  MOOCHO_STR="#"
fi

echo "sundanceMoocho=(" $sundanceMoocho ")"
echo "moochoStr=(" $MOOCHO_STR ")"




# ---------------0--------------------------------------------------------
# Check for SUNDIALS
# ------------------------------------------------------------------------

AC_ARG_WITH(sundials,
[AC_HELP_STRING([--with-sundials=SUNDIALS_BUILD_DIR],[use SUNDIALS build directory (enables SUNDIALS)])],
[
SUNDIALS_STR=""
SUNDIALS_BUILD_DIR=${withval}
AC_DEFINE(HAVE_SUNDIALS,,[define if we want to use SUNDIALS])
AC_MSG_CHECKING(SUNDIALS directory)
AC_MSG_RESULT([${SUNDIALS_BUILD_DIR}])
],
[
SUNDIALS_STR="#"
AC_MSG_CHECKING(SUNDIALS directory)
AC_MSG_RESULT([no])
]
)

# ------------------------------------------------------------------------
# Check for Exodus
# ------------------------------------------------------------------------

AC_ARG_WITH(sundance-exodus,
[AC_HELP_STRING([--with-sundance-exodus],[use Exodus (enables Exodus I/O in Sundance)])],
[
EXODUS_STR=""
AC_DEFINE(HAVE_EXODUS,,[define if we want to use EXODUS I/O in Sundance])
AC_MSG_CHECKING(have EXODUS)
AC_MSG_RESULT([yes])
],
[
EXODUS_STR="#"
AC_MSG_CHECKING(have EXODUS)
AC_MSG_RESULT([no])
]
)


# ------------------------------------------------------------------------
# Get the absolute path to the build directory.  
# ------------------------------------------------------------------------

AC_MSG_CHECKING([finding build directory])
if test "X" = "X"; then
SUNDANCE_BUILD_DIRECTORY=`pwd`
AC_MSG_RESULT($SUNDANCE_BUILD_DIRECTORY)
fi

# ------------------------------------------------------------------------
# Enable optional python wrappers
# ------------------------------------------------------------------------
TAC_ENABLE_PYTHON

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)
AC_SUBST(prefix)
AC_SUBST(SUNDANCE_BUILD_DIRECTORY)
AC_SUBST(TRILINOS_SRC_DIR)
AC_SUBST(MOOCHO_STR)
AC_SUBST(FIAT_STR)
AC_SUBST(FIAT_BUILD_DIR)
AC_SUBST(SUNDIALS_STR)
AC_SUBST(SUNDIALS_BUILD_DIR)
AC_SUBST(EXODUS_STR)
AC_SUBST(TRILINOS_DIR)


# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
Makefile  
Makefile.internal
Makefile.export
Makefile.export.TSFExtended
Makefile.site.sundance
examples-tutorial/Makefile
examples-tutorial/c++/Makefile
examples-tutorial/python/Makefile
src-utils/Makefile
src-interop/Makefile
src-core/Makefile
src-std-mesh/Makefile
src-std-framework/Makefile
src-solvers/Makefile
src-solvers/SundanceDefaultPath.hpp
python/Makefile
python/src/Makefile
python/src/SUNDANCE_HOME_DIR
python/src/__init__.py
python/src/.deps/dummy
python/example/Makefile
test/Makefile
tests-utils/Makefile
tests-utils/Combinatorics/Makefile
  tests-utils/Functors/Makefile
  tests-utils/Maps/Makefile
  tests-core/Makefile
  tests-core/SymbolicTests/Makefile
  tests-solvers/Makefile
  tests-solvers/VectorTests/Makefile
  tests-solvers/OperatorTests/Makefile
  tests-solvers/SolverTests/Makefile
  tests-std-mesh/Makefile
  tests-std-mesh/Simple/Makefile
  tests-std-mesh/Readers/Makefile
  tests-std-framework/Makefile
  tests-std-framework/Assembly/Makefile
  tests-std-framework/Problem/Makefile
])

AC_SUBST(CHACO_LIB)
AC_SUBST(TRILINOS_DIR)
AC_SUBST(TRILINOS_DATA_DIR)

AC_OUTPUT()


# Bye World!
echo "----------------------------------------"
echo "Finished Running Sundance Configure Script"
echo "----------------------------------------"
