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



noinst_PROGRAMS = $(TESTS_SRC:%=%.exe)

TESTS_SRC = PostPotentialFlow RadiationDiffusion1D PressureVessel2D \
VSNS2D Exo2Triangle PressureStabilizedNS2D PressureStabilizedNS3D \
ReadField CubicCDR ControlCDR MixedOrderPrec ThermalLTI

#  get_matrices_for_ross

TESTS = $(noinst_PROGRAMS:%.exe=%.tst)

REQUIREMENTS = $(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a \
$(top_builddir)/src-solvers/libTSFExtended.a



PostPotentialFlow_SOURCES = $(srcdir)/PostPotentialFlow.cpp

PostPotentialFlow_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

PostPotentialFlow_CXXFLAGS = $(STD_FWK_INCLUDES)

PostPotentialFlow_LDADD = $(STD_FWK_LIBS) 


Exo2Triangle_SOURCES = $(srcdir)/Exo2Triangle.cpp

Exo2Triangle_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

Exo2Triangle_CXXFLAGS = $(STD_FWK_INCLUDES)

Exo2Triangle_LDADD = $(STD_FWK_LIBS) 




RadiationDiffusion1D_SOURCES = $(srcdir)/RadiationDiffusion1D.cpp

RadiationDiffusion1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

RadiationDiffusion1D_CXXFLAGS = $(STD_FWK_INCLUDES)

RadiationDiffusion1D_LDADD = $(STD_FWK_LIBS) 



PressureVessel2D_SOURCES = $(srcdir)/PressureVessel2D.cpp

PressureVessel2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

PressureVessel2D_CXXFLAGS = $(STD_FWK_INCLUDES)

PressureVessel2D_LDADD = $(STD_FWK_LIBS) 




VSNS2D_SOURCES = $(srcdir)/VSNS2D.cpp

VSNS2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

VSNS2D_CXXFLAGS = $(STD_FWK_INCLUDES)

VSNS2D_LDADD = $(STD_FWK_LIBS) 





PressureStabilizedNS2D_SOURCES = $(srcdir)/PressureStabilizedNS2D.cpp

PressureStabilizedNS2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

PressureStabilizedNS2D_CXXFLAGS = $(STD_FWK_INCLUDES)

PressureStabilizedNS2D_LDADD = $(STD_FWK_LIBS) 




PressureStabilizedNS3D_SOURCES = $(srcdir)/PressureStabilizedNS3D.cpp

PressureStabilizedNS3D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

PressureStabilizedNS3D_CXXFLAGS = $(STD_FWK_INCLUDES)

PressureStabilizedNS3D_LDADD = $(STD_FWK_LIBS) 



ReadField_SOURCES = $(srcdir)/ReadField.cpp

ReadField_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

ReadField_CXXFLAGS = $(STD_FWK_INCLUDES)

ReadField_LDADD = $(STD_FWK_LIBS) 




CubicCDR_SOURCES = $(srcdir)/CubicCDR.cpp

CubicCDR_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

CubicCDR_CXXFLAGS = $(STD_FWK_INCLUDES)

CubicCDR_LDADD = $(STD_FWK_LIBS) 




ControlCDR_SOURCES = $(srcdir)/ControlCDR.cpp

ControlCDR_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

ControlCDR_CXXFLAGS = $(STD_FWK_INCLUDES)

ControlCDR_LDADD = $(STD_FWK_LIBS) 




# get_matrices_for_ross_SOURCES = $(srcdir)/get_matrices_for_ross.cpp
# 
# get_matrices_for_ross_DEPENDENCIES = \
# $(top_builddir)/src-std-framework/libsundanceStdFwk.a \
# $(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
# $(top_builddir)/src-core/libsundanceCore.a \
# $(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a
# 
# get_matrices_for_ross_CXXFLAGS = $(STD_FWK_INCLUDES)
# 
# get_matrices_for_ross_LDADD = $(STD_FWK_LIBS) 




MixedOrderPrec_SOURCES = $(srcdir)/MixedOrderPrec.cpp

MixedOrderPrec_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

MixedOrderPrec_CXXFLAGS = $(STD_FWK_INCLUDES)

MixedOrderPrec_LDADD = $(STD_FWK_LIBS) 




ThermalLTI_SOURCES = $(srcdir)/ThermalLTI.cpp

ThermalLTI_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

ThermalLTI_CXXFLAGS = $(STD_FWK_INCLUDES)

ThermalLTI_LDADD = $(STD_FWK_LIBS) 





include $(top_srcdir)/Makefile.common
