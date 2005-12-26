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

noinst_PROGRAMS = $(PARALLEL_SRC) $(SERIAL_SRC)

SERIAL_SRC = Poisson3D PoissonOnDisk Eikonal2D BesselTest2D FunctionalTest \
Kepler Stokes2D

PARALLEL_SRC = Poisson1D Poisson2D HighOrderPoisson2D \
NoxPoisson1D PoissonBoltzmann2D \
RadDiff1D NoxAlgebraic \
Coupled1D Coupled2D Helmholtz1D VorticityNavierStokes2D \
VorticityStokes2D AdvectionDiffusion2D LinearVariationalTest OptPoisson1D

PARALLEL_PROGS = $(PARALLEL_SRC:%=%.exe)






Poisson1D_SOURCES = $(srcdir)/Poisson1D.cpp

Poisson1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Poisson1D_CXXFLAGS = $(STD_FWK_INCLUDES)

Poisson1D_LDADD = $(STD_FWK_LIBS) 



OptPoisson1D_SOURCES = $(srcdir)/OptPoisson1D.cpp

OptPoisson1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

OptPoisson1D_CXXFLAGS = $(STD_FWK_INCLUDES)

OptPoisson1D_LDADD = $(STD_FWK_LIBS) 





Poisson2D_SOURCES = $(srcdir)/Poisson2D.cpp

Poisson2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Poisson2D_CXXFLAGS = $(STD_FWK_INCLUDES)

Poisson2D_LDADD = $(STD_FWK_LIBS) 





HighOrderPoisson2D_SOURCES = $(srcdir)/HighOrderPoisson2D.cpp

HighOrderPoisson2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

HighOrderPoisson2D_CXXFLAGS = $(STD_FWK_INCLUDES)

HighOrderPoisson2D_LDADD = $(STD_FWK_LIBS) 






Poisson3D_SOURCES = $(srcdir)/Poisson3D.cpp

Poisson3D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Poisson3D_CXXFLAGS = $(STD_FWK_INCLUDES)

Poisson3D_LDADD = $(STD_FWK_LIBS) 







Stokes2D_SOURCES = $(srcdir)/Stokes2D.cpp

Stokes2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Stokes2D_CXXFLAGS = $(STD_FWK_INCLUDES)

Stokes2D_LDADD = $(STD_FWK_LIBS) 










Coupled1D_SOURCES = $(srcdir)/Coupled1D.cpp

Coupled1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Coupled1D_CXXFLAGS = $(STD_FWK_INCLUDES)

Coupled1D_LDADD = $(STD_FWK_LIBS) 




Coupled2D_SOURCES = $(srcdir)/Coupled2D.cpp

Coupled2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Coupled2D_CXXFLAGS = $(STD_FWK_INCLUDES)

Coupled2D_LDADD = $(STD_FWK_LIBS) 


RadDiff1D_SOURCES = $(srcdir)/RadDiff1D.cpp

RadDiff1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

RadDiff1D_CXXFLAGS = $(STD_FWK_INCLUDES)

RadDiff1D_LDADD = $(STD_FWK_LIBS) 





NoxAlgebraic_SOURCES = $(srcdir)/NoxAlgebraic.cpp

NoxAlgebraic_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

NoxAlgebraic_CXXFLAGS = $(STD_FWK_INCLUDES)

NoxAlgebraic_LDADD = $(STD_FWK_LIBS) 



NoxPoisson1D_SOURCES = $(srcdir)/NoxPoisson1D.cpp

NoxPoisson1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

NoxPoisson1D_CXXFLAGS = $(STD_FWK_INCLUDES)

NoxPoisson1D_LDADD = $(STD_FWK_LIBS) 



Kepler_SOURCES = $(srcdir)/Kepler.cpp

Kepler_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Kepler_CXXFLAGS = $(STD_FWK_INCLUDES)

Kepler_LDADD = $(STD_FWK_LIBS) 

PoissonOnDisk_SOURCES = $(srcdir)/PoissonOnDisk.cpp

PoissonOnDisk_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

PoissonOnDisk_CXXFLAGS = $(STD_FWK_INCLUDES)

PoissonOnDisk_LDADD = $(STD_FWK_LIBS) 










Eikonal2D_SOURCES = $(srcdir)/Eikonal2D.cpp $(srcdir)/bessi0.f

Eikonal2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Eikonal2D_CXXFLAGS = $(STD_FWK_INCLUDES)

Eikonal2D_LDADD = $(STD_FWK_LIBS) 









#TCATA1D_SOURCES = $(srcdir)/TCATA1D.cpp
#
#TCATA1D_DEPENDENCIES = \
#$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
#$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
#$(top_builddir)/src-core/libsundanceCore.a \
#$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a
#
#TCATA1D_CXXFLAGS = $(STD_FWK_INCLUDES)
#
#TCATA1D_LDADD = $(STD_FWK_LIBS) 










VorticityStokes2D_SOURCES = $(srcdir)/VorticityStokes2D.cpp

VorticityStokes2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

VorticityStokes2D_CXXFLAGS = $(STD_FWK_INCLUDES)

VorticityStokes2D_LDADD = $(STD_FWK_LIBS) 











Helmholtz1D_SOURCES = $(srcdir)/Helmholtz1D.cpp

Helmholtz1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Helmholtz1D_CXXFLAGS = $(STD_FWK_INCLUDES)

Helmholtz1D_LDADD = $(STD_FWK_LIBS) 










VorticityNavierStokes2D_SOURCES = $(srcdir)/VorticityNavierStokes2D.cpp

VorticityNavierStokes2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

VorticityNavierStokes2D_CXXFLAGS = $(STD_FWK_INCLUDES)

VorticityNavierStokes2D_LDADD = $(STD_FWK_LIBS) 






BesselTest2D_SOURCES = $(srcdir)/BesselTest2D.cpp 

BesselTest2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

BesselTest2D_CXXFLAGS = $(STD_FWK_INCLUDES)

BesselTest2D_LDADD = $(STD_FWK_LIBS) 



AdvectionDiffusion2D_SOURCES = $(srcdir)/AdvectionDiffusion2D.cpp 

AdvectionDiffusion2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

AdvectionDiffusion2D_CXXFLAGS = $(STD_FWK_INCLUDES)

AdvectionDiffusion2D_LDADD = $(STD_FWK_LIBS) 


FunctionalTest_SOURCES = $(srcdir)/FunctionalTest.cpp 

FunctionalTest_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

FunctionalTest_CXXFLAGS = $(STD_FWK_INCLUDES)

FunctionalTest_LDADD = $(STD_FWK_LIBS) 





LinearVariationalTest_SOURCES = $(srcdir)/LinearVariationalTest.cpp 

LinearVariationalTest_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

LinearVariationalTest_CXXFLAGS = $(STD_FWK_INCLUDES)

LinearVariationalTest_LDADD = $(STD_FWK_LIBS) 




PoissonBoltzmann2D_SOURCES = $(srcdir)/PoissonBoltzmann2D.cpp

PoissonBoltzmann2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

PoissonBoltzmann2D_CXXFLAGS = $(STD_FWK_INCLUDES)

PoissonBoltzmann2D_LDADD = $(STD_FWK_LIBS) 


TESTS = $(noinst_PROGRAMS:%.exe=%.tst)

PTESTS = $(PARALLEL_PROGS:%.exe=%.ptst)

@HAVE_MPI_TRUE@CHECKTEST_ARGS = $(TESTS) $(PTESTS)
@HAVE_MPI_FALSE@CHECKTEST_ARGS = $(TESTS)


NPROC = 2

@HAVE_MPI_TRUE@MPI_COMMAND = mpirun -np $(NPROC)
@HAVE_MPI_FALSE@MPI_COMMAND = 

%.tst : %.exe
	./$< > $@ 2>&1 
	echo "done" $@

%.ptst : %.exe
	$(MPI_COMMAND) ./$< > $@ 2>&1 
	echo "done" $@

tests:
	$(MAKE) $(TESTS)
	@HAVE_MPI_TRUE@$(MAKE) $(PTESTS)
	$(top_srcdir)/etc/checkTests.pl $(CHECKTEST_ARGS) 

serialtests:
	$(MAKE) $(TESTS)
	$(top_srcdir)/etc/checkTests.pl $(CHECKTEST_ARGS) 


include $(top_srcdir)/Makefile.common
