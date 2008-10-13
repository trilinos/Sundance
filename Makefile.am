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
	config/strip_dup_incl_paths.pl config/strip_dup_libs.pl \
	config/ac_cxx_have_std.m4 config/ac_cxx_have_stl.m4 \
	config/ac_cxx_namespaces.m4 config/ac_cxx_mutable.m4 \
	config/depcomp \
	etc/SolverParameters/amesos.xml etc/SolverParameters/aztec-ifpack.xml \
	etc/SolverParameters/aztec-ml.xml \
	etc/SolverParameters/aztec-native.xml \
	etc/SolverParameters/aztec.xml etc/SolverParameters/bicgstab.xml \
	etc/SolverParameters/gmres.xml etc/SolverParameters/nox.xml \
	etc/SolverParameters/poissonParams.xml \
	etc/SolverParameters/userPrecParams.xml etc/SundanceConfig.xml \
	etc/checkTests.pl

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
#  We now build tests and examples through separate make targets, rather 
#  than during "make".  We still need to conditionally include the test
#  and example subdirectories in SUBDIRS, even though SUB_TEST and 
#  SUB_EXAMPLE will never be defined, so that the tests and examples
#  are included in the distribution tarball. 

if SUB_EXAMPLE
EXAMPLE_SUBDIR=examples-tutorial
endif



# ----------------------------------------------------------------------
if SUB_TEST
TEST_SUBDIR = tests-solvers tests-utils tests-core tests-std-mesh tests-std-framework test
endif

if BUILD_TESTS
tests: 
	@echo ""
	@echo "Now building Sundance tests"
	@echo ""
	cd $(top_builddir)/tests-solvers && $(MAKE)
	cd $(top_builddir)/tests-utils && $(MAKE)
	cd $(top_builddir)/tests-core && $(MAKE)
	cd $(top_builddir)/tests-std-mesh && $(MAKE)
	cd $(top_builddir)/tests-std-framework && $(MAKE)
	@echo ""
	@echo "Now copying data needed for Sundance tests to build directory"
	@echo ""
	$(MAKE) get-data-needed-for-tests
	@echo ""
	@echo "Finished building Sundance tests."
	@echo ""

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


SUBDIRS = src-solvers src-utils src-interop src-core src-std-mesh src-std-framework $(PYTHON_SUBDIRS) $(TEST_SUBDIR) $(EXAMPLE_SUBDIR)



## #######################################################################
## Make libraries w/o making tests
## #######################################################################

libonly:
	(cd src-solvers; $(MAKE))
	(cd src-utils; $(MAKE))
	(cd src-interop; $(MAKE))
	(cd src-core; $(MAKE))
	(cd src-std-mesh; $(MAKE))
	(cd src-std-framework; $(MAKE))
	@echo ""
	@echo "Now copying data needed for Sundance tests to build directory"
	@echo ""
	$(MAKE) get-data-needed-for-tests


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
	mkdir -p $(DESTDIR)$(sysconfdir)/SolverParameters
	mkdir -p $(DESTDIR)$(datadir)/Sundance/examples-tutorial/meshes
	mkdir -p $(DESTDIR)$(datadir)/Sundance/examples-tutorial/data
	cp $(top_srcdir)/etc/SolverParameters/*.xml $(DESTDIR)$(sysconfdir)/SolverParameters
	cp $(top_srcdir)/etc/*.xml $(DESTDIR)$(sysconfdir)
	cp $(top_srcdir)/examples-tutorial/data/*.* $(DESTDIR)$(datadir)/Sundance/examples-tutorial/data



## #######################################################################
## hooks for putting data files (meshes, solver params, etc) in the
## build directory
## #######################################################################

get-data-needed-for-tests:
	mkdir -p $(abs_top_builddir)/etc/SolverParameters
	mkdir -p $(abs_top_builddir)/examples-tutorial/meshes
	mkdir -p $(abs_top_builddir)/examples-tutorial/data
	cp $(top_srcdir)/etc/SolverParameters/*.xml $(abs_top_builddir)/etc/SolverParameters
	cp $(top_srcdir)/etc/*.xml $(abs_top_builddir)/etc
	cp $(top_srcdir)/examples-tutorial/data/*.* $(abs_top_builddir)/examples-tutorial/data


## #######################################################################
## The below targets allow you to use the new
## testharness to run the test suite as make targets
## #######################################################################


if HAVE_MPI
TEUCHOS_CHECK_COMM=mpi
else
TEUCHOS_CHECK_COMM=serial
endif


if BUILD_EXAMPLES
examples:
	@echo ""
	@echo "Now building Sundance examples."
	@echo ""
	cd $(top_builddir)/examples-tutorial && $(MAKE)
	@echo ""
	@echo "Finished building Sundance examples."
	@echo ""
install-examples:
	cd $(top_builddir)/examples-tutorial && $(MAKE) install
else
examples:
	@echo "Sundance examples were disabled at configure time"

install-examples:
	@echo "Sundance examples were disabled at configure time"
endif

clean-tests:
	cd $(top_builddir)/tests-solvers && $(MAKE) clean
	cd $(top_builddir)/tests-core && $(MAKE) clean
	cd $(top_builddir)/tests-std-framework && $(MAKE) clean

clean-examples:
	cd $(top_builddir)/examples-tutorial && $(MAKE) clean

everything:
	$(MAKE) && $(MAKE) examples && $(MAKE) tests

clean-everything:
	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

install-everything:
	$(MAKE) install && $(MAKE) install-examples




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





