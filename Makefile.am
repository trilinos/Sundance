#np# Typically the only line that needs to be edited in this file is the last 
#np# line.  Set 'SUBDIRS= (list all subdirectories that are autotool'ed here)'
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.

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

## #######################################################################
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Aclocal command (rarely used - don't worry about it)
## #######################################################################
ACLOCAL	= aclocal -I config

#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = bootstrap \
	config/ac_cxx_have_std.m4 config/ac_cxx_have_stl.m4 \
	config/ac_cxx_namespaces.m4 config/ac_cxx_mutable.m4 \
	config/depcomp

#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs 
#
#  Again, I hope that AUX_DIR(../../config) eliminates these
#  config/install-sh config/missing config/mkinstalldirs 

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	src/common/config-h.in src/common/stamp-h.in \
	$(AUX_DIST) 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################
#  #np# - Most packages should have an example and a full suite 
#  of tests.  At present, new_package only offers src and example.
# SUBDIRS = src example test 

if BUILD_EXAMPLES
EXAMPLES_SUBDIR=examples-tutorial
endif

if PYTHON_USE
PYTHON_SUBDIRS=python
PYTESTS=(cd python/example; $(MAKE) tests)
PYSERIALTESTS=(cd python/example; $(MAKE) serialtests)
endif

if BUILD_CXXTESTS

SUBDIRS = src-solvers tests-solvers \
src-utils tests-utils \
src-core tests-core \
src-std-mesh tests-std-mesh \
src-std-framework tests-std-framework \
$(PYTHON_SUBDIRS) src-xml  $(EXAMPLES_SUBDIR)

CXX_SOLVER_TESTS = (cd tests-solvers; $(MAKE) tests) 
CXX_SOLVER_SERIALTESTS = (cd tests-solvers; $(MAKE) serialtests) 

CXX_SYMBOLIC_TESTS = (cd tests-core; $(MAKE) tests) 
CXX_SYMBOLIC_SERIALTESTS = (cd tests-core; $(MAKE) serialtests) 

CXX_FWK_TESTS = (cd tests-std-framework; $(MAKE) tests) 
CXX_FWK_SERIALTESTS = (cd tests-std-framework; $(MAKE) serialtests) 

else
SUBDIRS = src-solvers src-utils src-core src-std-mesh src-std-framework \
$(PYTHON_SUBDIRS) src-xml  $(EXAMPLES_SUBDIR)
endif



## #######################################################################
## Make libraries w/o making tests
## #######################################################################

libonly:
	(cd src-solvers; $(MAKE) install)
	(cd src-utils; $(MAKE) install)
	(cd src-core; $(MAKE) install)
	(cd src-std-mesh; $(MAKE) install)
	(cd src-std-framework; $(MAKE) install)

## #######################################################################
## emacs tags
## #######################################################################

tags:
	etags */*/*.hpp */*/*.cpp

tests:
	@echo "--------------------------------------------------------------------"
	@echo "running all tests"
	@echo "--------------------------------------------------------------------"
	$(CXX_SOLVER_TESTS)
	$(CXX_SYMBOLIC_TESTS)
	$(CXX_FWK_TESTS)
	$(PYTESTS)
	@echo "====================================================================="
	@echo "Summary of serial tests:"
	@echo "====================================================================="
	$(top_srcdir)/etc/checkTests.pl tests-*/*/*.tst python/example/*.tst
	@echo "====================================================================="
	@echo "Summary of parallel tests:"
	@echo "====================================================================="
	$(top_srcdir)/etc/checkTests.pl tests-*/*/*.ptst python/example/*.ptst

serialtests:
	$(CXX_SOLVER_SERIALTESTS)
	$(CXX_SYMBOLIC_SERIALTESTS)
	$(CXX_FWK_SERIALTESTS)
	$(PYSERIALTESTS)
	@echo "====================================================================="
	@echo "Summary of serial tests:"
	@echo "====================================================================="
	$(top_srcdir)/etc/checkTests.pl tests-*/*/*.tst python/example/*.tst
