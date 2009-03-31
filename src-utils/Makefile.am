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

#np# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
#np# There is, at present, no distinction between most of the categories, 
#np# other than those which start with EXTRA_.  All categories which do not
#np# start with EXTRA_ are included in all builds.  Categories which start
#np# with EXTRA_ are included in some builds (typically controlled by 
#np# --with-xyz or --without-xyz)
#
#np# Make sure to list all source files in one of the following categories.
SUNDANCE_UTILITIES_H = \
$(top_srcdir)/src-utils/Utilities/SundanceADReal.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceArrayOfTuples.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceCollectiveExceptionCheck.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceCombinatorialUtils.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceDebug.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceDefs.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceExceptions.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceFunctorDomain.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceIntHashSet.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceMap.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceMultiSet.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceNoncopyable.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceObjectWithInstanceID.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceOrderedHandle.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceOrderedTuple.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceOut.hpp \
$(top_srcdir)/src-utils/Utilities/SundancePoint.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceSet.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceStdMathFunctors.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceTabs.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceUnaryFunctor.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceVersionString.hpp \
$(top_srcdir)/src-utils/Utilities/SundanceWatchFlag.hpp 

SUNDANCE_UTILITIES = \
$(top_srcdir)/src-utils/Utilities/SundanceADReal.cpp \
$(top_srcdir)/src-utils/Utilities/SundanceCollectiveExceptionCheck.cpp \
$(top_srcdir)/src-utils/Utilities/SundanceCombinatorialUtils.cpp \
$(top_srcdir)/src-utils/Utilities/SundanceExceptions.cpp \
$(top_srcdir)/src-utils/Utilities/SundanceFunctorDomain.cpp \
$(top_srcdir)/src-utils/Utilities/SundanceIntHashSet.cpp \
$(top_srcdir)/src-utils/Utilities/SundancePoint.cpp \
$(top_srcdir)/src-utils/Utilities/SundanceStdMathFunctors.cpp \
$(top_srcdir)/src-utils/Utilities/SundanceTabs.cpp \
$(top_srcdir)/src-utils/Utilities/SundanceUnaryFunctor.cpp \
$(top_srcdir)/src-utils/Utilities/SundanceWatchFlag.cpp  





UTIL_SOURCES = \
$(SUNDANCE_UTILITIES) 

UTIL_SOURCES_H = \
$(SUNDANCE_UTILITIES_H) 

include $(top_srcdir)/Makefile.common

#
#np#  Files which are included in select builds are typically have
#np#  names that start with EXTRA
#
#np#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#np#  set up by TAC_ARC_CONFIG_MPI 
#
EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif

AM_CPPFLAGS = $(UTIL_INCLUDES) $(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# Sundance library specifications
# ------------------------------------------------------------------------

UTIL_LIB = libsundanceUtils.a


#np# replace new_package with the name of the package being autotool'ed here
libsundanceUtils_a_SOURCES = $(UTIL_SOURCES)

#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libsundanceUtils_a_SOURCES = \
$(EXTRA_MPI) \
$(EXTRA_MPI_H) 

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libsundance_a_AR = $(ALTERNATE_AR)
else

libsundance_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libsundanceUtils_a_LIBADD = $(UTIL_TEMPLATE_OBJS) $(XTRALDADD) 

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(UTIL_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(UTIL_SOURCES_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: SundanceConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(includedir)
nodist_execinclude_HEADERS = SundanceConfig.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






