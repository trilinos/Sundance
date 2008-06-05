#@HEADER
# ************************************************************************
# 
#           TSFExtended: Trilinos Solver Framework Extended
#                 Copyright (2004) Sandia Corporation
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
#@HEADER

include $(top_srcdir)/Makefile.common

EXEEXT = .exe

PARALLEL_SRC = CompoundOpTest BlockOpTest LinearCombinationTest LTITest MultiVectorOpTest

SERIAL_SRC =

ALL_SRC = $(PARALLEL_SRC) $(SERIAL_SRC)

PARALLEL_PROGS = $(PARALLEL_SRC:%=%.exe)

SERIAL_PROGS = $(SERIAL_SRC:%=%.exe)

noinst_PROGRAMS = $(ALL_SRC:%=%.exe)




CompoundOpTest_SOURCES = $(srcdir)/CompoundOpTest.cpp 

CompoundOpTest_DEPENDENCIES = $(top_builddir)/src-solvers/libTSFExtended.a

CompoundOpTest_CXXFLAGS = $(CXXFLAGS) $(SOLVERS_INCLUDES) 

CompoundOpTest_LDADD = $(LDFLAGS) $(SOLVERS_LIBS) 



MultiVectorOpTest_SOURCES = $(srcdir)/MultiVectorOpTest.cpp 

MultiVectorOpTest_DEPENDENCIES = $(top_builddir)/src-solvers/libTSFExtended.a

MultiVectorOpTest_CXXFLAGS = $(CXXFLAGS) $(SOLVERS_INCLUDES) 

MultiVectorOpTest_LDADD = $(LDFLAGS) $(SOLVERS_LIBS) 




BlockOpTest_SOURCES = $(srcdir)/BlockOpTest.cpp 

BlockOpTest_DEPENDENCIES = $(top_builddir)/src-solvers/libTSFExtended.a

BlockOpTest_CXXFLAGS = $(CXXFLAGS) $(SOLVERS_INCLUDES) 

BlockOpTest_LDADD = $(LDFLAGS) $(SOLVERS_LIBS) 



LinearCombinationTest_SOURCES = $(srcdir)/LinearCombinationTest.cpp 

LinearCombinationTest_DEPENDENCIES = $(top_builddir)/src-solvers/libTSFExtended.a

LinearCombinationTest_CXXFLAGS = $(CXXFLAGS) $(SOLVERS_INCLUDES) 

LinearCombinationTest_LDADD = $(LDFLAGS) $(SOLVERS_LIBS) 




LTITest_SOURCES = $(srcdir)/LTITest.cpp 

LTITest_DEPENDENCIES = $(top_builddir)/src-solvers/libTSFExtended.a

LTITest_CXXFLAGS = $(CXXFLAGS) $(SOLVERS_INCLUDES) 

LTITest_LDADD = $(LDFLAGS) $(SOLVERS_LIBS) 





STESTS = $(ALL_SRC:%=%.tst)

PTESTS = $(PARALLEL_SRC:%=%.ptst)


@HAVE_MPI_TRUE@CHECKTEST_ARGS = $(STESTS) $(PTESTS)
@HAVE_MPI_FALSE@CHECKTEST_ARGS = $(STESTS)


NPROC = 2

@HAVE_MPI_TRUE@MPI_COMMAND = mpirun -np $(NPROC)
@HAVE_MPI_FALSE@MPI_COMMAND = 

%.tst : %.exe
	./$< > $@ 2>&1 
	echo "done" $@

%.ptst : %.exe
	$(MPI_COMMAND) ./$< > $@ 2>&1 
	echo "done" $@

old-tests:
	$(MAKE) $(STESTS)
	@HAVE_MPI_TRUE@$(MAKE) $(PTESTS)
	$(top_srcdir)/etc/checkTests.pl $(CHECKTEST_ARGS) 

old-serialtests:
	$(MAKE) $(STESTS)
	$(top_srcdir)/etc/checkTests.pl $(CHECKTEST_ARGS) 



