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

#The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package Sundance built successfully."
	@echo ""

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################
#  #np# - Most packages should have an example and a full suite 
#  of tests.  At present, new_package only offers src and example.
# SUBDIRS = src example test 

if BUILD_EXAMPLES
EXAMPLES_SUBDIR=examples-tutorial
endif



# ----------------------------------------------------------------------
if BUILD_TESTS

SOLVERS_TESTS = tests-solvers
UTILS_TESTS = tests-utils
CORE_TESTS = tests-core
STD_MESH_TESTS = tests-std-mesh
STD_FRAMEWORK_TESTS = tests-std-framework

tests: 
	@echo ""
	@echo "Now building Sundance tests"
	@echo ""
	cd $(top_builddir)/tests-solvers && $(MAKE)
	cd $(top_builddir)/tests-utils && $(MAKE)
	cd $(top_builddir)/tests-core && $(MAKE)
	cd $(top_builddir)/tests-std-mesh && $(MAKE)
	cd $(top_builddir)/tests-std-framework && $(MAKE)
else
tests:
	@echo ""
	@echo "Sundance tests disabled at configure time"
	@echo ""
endif
## --------------- end of BUILD_TESTS section --------------------------


if PYTHON_USE
PYTHON_SUBDIRS=python
endif


SUBDIRS = src-solvers src-utils src-interop src-core src-std-mesh src-std-framework $(PYTHON_SUBDIRS) 



## #######################################################################
## Make libraries w/o making tests
## #######################################################################

libonly:
	(cd src-solvers; $(MAKE) install)
	(cd src-utils; $(MAKE) install)
	(cd src-interop; $(MAKE) install)
	(cd src-core; $(MAKE) install)
	(cd src-std-mesh; $(MAKE) install)
	(cd src-std-framework; $(MAKE) install)


## #######################################################################
## emacs tags
## #######################################################################

tags:
	etags */*/*.hpp */*/*.cpp



## #######################################################################
## hooks for putting data files (meshes, solver params, etc) in the
## installation directory
## #######################################################################

install-data-hook:
	mkdir -p $(prefix)/etc/SolverParameters
	mkdir -p $(prefix)/examples-tutorial/meshes
	mkdir -p $(prefix)/examples-tutorial/data
	cp $(top_srcdir)/etc/SolverParameters/*.xml $(prefix)/etc/SolverParameters
	cp $(top_srcdir)/etc/*.xml $(prefix)/etc
	cp $(top_srcdir)/examples-tutorial/data/*.* $(prefix)/examples-tutorial/data


## #######################################################################
## The below targets allow you to use the new
## testharness to run the test suite as make targets
## #######################################################################


if HAVE_MPI
TEUCHOS_CHECK_COMM=mpi
else
TEUCHOS_CHECK_COMM=serial
endif




old-tests:
	(cd tests-solvers; $(MAKE) old-tests)
	(cd tests-core; $(MAKE) old-tests)
	(cd tests-std-framework; $(MAKE) old-tests)

old-serialtests:
	(cd tests-solvers; $(MAKE) old-serialtests)
	(cd tests-core; $(MAKE) old-serialtests)
	(cd tests-std-framework; $(MAKE) old-serialtests)

if TRILINOS_IS_EXTERNAL

else
TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..

TRILINOS_TEST_CATEGORY=INSTALL


runtests-serial :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=serial \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --category=$(TRILINOS_TEST_CATEGORY) \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=Sundance
	$(RUN_PYSUNDANCE_SERIAL_TESTS)

runtests-mpi :
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=mpi \
  --mpi-go=$(TRILINOS_MPI_GO) \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --max-proc=$(TRILINOS_MPI_MAX_PROC) \
  --category=$(TRILINOS_TEST_CATEGORY) \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=Sundance
	$(RUN_PYSUNDANCE_MPI_TESTS)


runtests-python-serial:
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=serial \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --category=PYSUNDANCE \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=Sundance


runtests-python-mpi:
	$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
  --trilinos-dir=$(TRILINOS_HOME_DIR) \
  --comm=mpi \
  --mpi-go=$(TRILINOS_MPI_GO) \
  --build-dir=$(TRILINOS_BUILD_DIR) \
  --max-proc=$(TRILINOS_MPI_MAX_PROC) \
  --category=PYSUNDANCE \
  --output-dir=@abs_top_builddir@/test/runtests-results \
  --verbosity=1 \
  --packages=Sundance

endif



