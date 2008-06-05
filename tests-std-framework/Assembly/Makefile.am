#@HEADER
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
#@HEADER

include $(top_srcdir)/Makefile.common

EXEEXT = .exe

SERIAL_SRC = BasisCheck TransformedIntegral2D  QuadratureTest 

ALL_SRC = $(PARALLEL_SRC) $(SERIAL_SRC) QuadratureIntegral2D

PARALLEL_PROGS = $(PARALLEL_SRC:%=%.exe)

SERIAL_PROGS = $(SERIAL_SRC:%=%.exe)

noinst_PROGRAMS = $(ALL_SRC:%=%.exe)




BasisCheck_SOURCES = $(srcdir)/BasisCheck.cpp

BasisCheck_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-interop/libsundanceInterop.a \
$(top_builddir)/src-utils/libsundanceUtils.a

BasisCheck_CXXFLAGS = $(STD_FWK_INCLUDES)

BasisCheck_LDADD = $(STD_FWK_LIBS) 



TransformedIntegral2D_SOURCES = $(srcdir)/TransformedIntegral2D.cpp

TransformedIntegral2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-interop/libsundanceInterop.a \
$(top_builddir)/src-utils/libsundanceUtils.a

TransformedIntegral2D_CXXFLAGS = $(STD_FWK_INCLUDES)

TransformedIntegral2D_LDADD = $(STD_FWK_LIBS) 




QuadratureIntegral2D_SOURCES = $(srcdir)/QuadratureIntegral2D.cpp

QuadratureIntegral2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-interop/libsundanceInterop.a \
$(top_builddir)/src-utils/libsundanceUtils.a

QuadratureIntegral2D_CXXFLAGS = $(STD_FWK_INCLUDES)

QuadratureIntegral2D_LDADD = $(STD_FWK_LIBS) 




QuadratureTest_SOURCES = $(srcdir)/QuadratureTest.cpp

QuadratureTest_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-interop/libsundanceInterop.a \
$(top_builddir)/src-utils/libsundanceUtils.a

QuadratureTest_CXXFLAGS = $(STD_FWK_INCLUDES)

QuadratureTest_LDADD = $(STD_FWK_LIBS) 


STESTS = $(SERIAL_SRC:%=%.tst)

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




