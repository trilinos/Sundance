# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                               Sundance
#                 Copyright (2002) Sandia Corporation
# 
#
# Copyright (2001) Sandia Corportation. Under the terms of Contract 
# DE-AC04-94AL85000, there is a non-exclusive license for use of this 
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
# 
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in ths data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
# 
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT(Sundance, 2.0, krlong@sandia.gov)

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(/src-core/UserExprs/SundanceExpr.hpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/Support/SundanceConfig.h 
# from src/Support/SundanceConfig.h.in

AM_CONFIG_HEADER(src-utils/SundanceConfig.h:src-utils/SundanceConfig.h.in)



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
TRILINOS_DIR=${withval}
AC_MSG_CHECKING(Trilinos root directory)
AC_MSG_RESULT([${TRILINOS_DIR}])
]
)

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

#  #np# - These options can disable the tests and examples of a package.
#  #np# - Packages that do not have tests or examples should #-out the 
#  #np# - option(s) that does (do) not apply.

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for Sundance], TESTS, yes)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for Sundance], EXAMPLES, yes)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_examples" != "Xno")

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

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
  Makefile  
  examples-tutorial/Makefile
  src-utils/Makefile
  src-core/Makefile
  src-std-mesh/Makefile
  src-std-framework/Makefile
  tests-core/Makefile
  tests-core/SymbolicTests/Makefile
  tests-std-mesh/Makefile
  tests-std-mesh/Simple/Makefile
  tests-std-mesh/Readers/Makefile
  tests-std-framework/Makefile
  tests-std-framework/HomogeneousDOFMap/Makefile
  tests-std-framework/Assembly/Makefile
  tests-std-framework/Problem/Makefile
  doc/Makefile
  doc/Doxygen/Makefile
  doc/Doxygen/Doxyfile.utils
  doc/Doxygen/Doxyfile.core
  doc/Doxygen/Doxyfile.trilinos

])

AC_SUBST(CHACO_LIB)
AC_SUBST(TRILINOS_DIR)

AC_OUTPUT()


# Bye World!
echo "----------------------------------------"
echo "Finished Running Sundance Configure Script"
echo "----------------------------------------"
