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


EXEEXT = .exe

SERIAL_SRC = TestEval VecOpTest UserDefTest 

ALL_SRC = $(SERIAL_SRC) SetTest ChainRule SimpleEval SensitivityTest \
PolynomialTest StdTransformation ScalabilityTest DepthScalabilityTest \
VecTimings VariationTest testExpr



PARALLEL_PROGS = $(PARALLEL_SRC:%=%.exe)

SERIAL_PROGS = $(SERIAL_SRC:%=%.exe)

noinst_PROGRAMS = $(ALL_SRC:%=%.exe)



LIBDEPS = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-interop/libsundanceInterop.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

SimpleEval_SOURCES = $(srcdir)/SimpleEval.cpp 

SimpleEval_DEPENDENCIES = $(LIBDEPS)

SimpleEval_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

SimpleEval_LDADD = $(CORE_LIBS) 



testExpr_SOURCES = $(srcdir)/testExpr.cpp 

testExpr_DEPENDENCIES = $(LIBDEPS)

testExpr_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

testExpr_LDADD = $(CORE_LIBS) 



ChainRule_SOURCES = $(srcdir)/ChainRule.cpp 

ChainRule_DEPENDENCIES = $(LIBDEPS)

ChainRule_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

ChainRule_LDADD = $(CORE_LIBS) 



UserDefTest_SOURCES = $(srcdir)/UserDefTest.cpp 

UserDefTest_DEPENDENCIES = $(LIBDEPS)

UserDefTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

UserDefTest_LDADD = $(CORE_LIBS) 




PolynomialTest_SOURCES = $(srcdir)/PolynomialTest.cpp 

PolynomialTest_DEPENDENCIES = $(LIBDEPS)

PolynomialTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

PolynomialTest_LDADD = $(CORE_LIBS) 



VecOpTest_SOURCES = $(srcdir)/VecOpTest.cpp 

VecOpTest_DEPENDENCIES = $(LIBDEPS)

VecOpTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

VecOpTest_LDADD = $(CORE_LIBS) 


VariationTest_SOURCES = $(srcdir)/VariationTest.cpp 

VariationTest_DEPENDENCIES = $(LIBDEPS)

VariationTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

VariationTest_LDADD = $(CORE_LIBS) 




VecTimings_SOURCES = $(srcdir)/VecTimings.cpp 

VecTimings_DEPENDENCIES = $(LIBDEPS)

VecTimings_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

VecTimings_LDADD = $(CORE_LIBS) 






ScalabilityTest_SOURCES = $(srcdir)/ScalabilityTest.cpp 

ScalabilityTest_DEPENDENCIES = $(LIBDEPS)

ScalabilityTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

ScalabilityTest_LDADD = $(CORE_LIBS) 



DepthScalabilityTest_SOURCES = $(srcdir)/DepthScalabilityTest.cpp 

DepthScalabilityTest_DEPENDENCIES = $(LIBDEPS)

DepthScalabilityTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

DepthScalabilityTest_LDADD = $(CORE_LIBS) 




TestEval_SOURCES = $(srcdir)/TestEval.cpp 

TestEval_DEPENDENCIES = $(LIBDEPS)

TestEval_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

TestEval_LDADD = $(CORE_LIBS) 


SetTest_SOURCES = $(srcdir)/SetTest.cpp 

SetTest_DEPENDENCIES = $(LIBDEPS)

SetTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

SetTest_LDADD = $(CORE_LIBS) 



SensitivityTest_SOURCES = $(srcdir)/SensitivityTest.cpp 

SensitivityTest_DEPENDENCIES = $(LIBDEPS)

SensitivityTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

SensitivityTest_LDADD = $(CORE_LIBS) 



StdTransformation_SOURCES = $(srcdir)/StdTransformation.cpp 

StdTransformation_DEPENDENCIES = $(LIBDEPS)

StdTransformation_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

StdTransformation_LDADD = $(CORE_LIBS) 


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


include $(top_srcdir)/Makefile.common
