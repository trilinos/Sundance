# @HEADER
# ************************************************************************
# 
#              PyTrilinos.Epetra: Python Interface to Epetra
#                   Copyright (2005) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER


all: $(ALL_SRC)
	if [ $(top_srcdir) != @abs_top_builddir@ ]; then /bin/cp $(srcdir)/$(SETPATH_PY) $(PARAM_FILES) . ; fi


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

PARAM_FILES = $(srcdir)/aztecSolver.py $(srcdir)/bicgstabSolver.py \
$(srcdir)/noxSolver.py $(srcdir)/amesosSolver.py $(srcdir)/blockTriangularSolver.py 

EXTRA_DIST = $(PARAM_FILES) setpath.py $(ALL_PROGS)

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in

# Define the suite of examples
PARALLEL_SRC = exPoisson1D exPoisson2D exHelmholtz1D exRadDiff1D exConvTest2D exBlockTest \
exVorticityStokes2D exComplexHelmholtz1D 

SERIAL_SRC = exPoisson3D exCylPoisson3D exTwoZoneLaplace1D exNonlinearPartialDomain \
exAToC

ALL_SRC = $(PARALLEL_SRC) $(SERIAL_SRC)



SETPATH_PY = setpath.py 


TESTS = $(ALL_PROGS:%.py=%.tst) 

PTESTS = $(PARALLEL_PROGS:%.py=%.ptst)

ALL_PROGS = $(PARALLEL_PROGS) $(SERIAL_PROGS)

PARALLEL_PROGS = $(PARALLEL_SRC:%=%.py)
SERIAL_PROGS = $(SERIAL_SRC:%=%.py)

@HAVE_MPI_TRUE@CHECKTEST_ARGS = $(TESTS) $(PTESTS)
@HAVE_MPI_FALSE@CHECKTEST_ARGS = $(TESTS)


NPROC = 2

@HAVE_MPI_TRUE@MPI_COMMAND = mpirun -np $(NPROC)
@HAVE_MPI_FALSE@MPI_COMMAND = 
   

%.tst: %.py
	./$< > $@ 2>&1
	echo "done" $@

%.ptst: %.py
	$(MPI_COMMAND) ./$< > $@ 2>&1
	echo "done" $@

tests:
	$(MAKE) all
	$(MAKE) $(TESTS)
	@HAVE_MPI_TRUE@$(MAKE) $(PTESTS)
	$(top_srcdir)/etc/checkTests.pl $(CHECKTEST_ARGS) 

serialtests:
	$(MAKE) all
	$(MAKE) $(TESTS)
	$(top_srcdir)/etc/checkTests.pl $(CHECKTEST_ARGS) 

ex%: $(srcdir)/ex%.py
	if [ $(top_srcdir) != @abs_top_builddir@ ]; then /bin/cp $(srcdir)/$@.py . ; fi
	chmod +x $@.py

clean:
	if [ $(top_srcdir) != @abs_top_builddir@ ]; then /bin/rm -f *.py ; fi
	/bin/rm -f *.pyc


