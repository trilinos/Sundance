#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = PostPotentialFlow RadiationDiffusion1D PressureVessel2D \
VSNS2D





PostPotentialFlow_SOURCES = $(srcdir)/PostPotentialFlow.cpp

PostPotentialFlow_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

PostPotentialFlow_CXXFLAGS = $(STD_FWK_INCLUDES)

PostPotentialFlow_LDADD = $(STD_FWK_LIBS) 




RadiationDiffusion1D_SOURCES = $(srcdir)/RadiationDiffusion1D.cpp

RadiationDiffusion1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

RadiationDiffusion1D_CXXFLAGS = $(STD_FWK_INCLUDES)

RadiationDiffusion1D_LDADD = $(STD_FWK_LIBS) 



PressureVessel2D_SOURCES = $(srcdir)/PressureVessel2D.cpp

PressureVessel2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

PressureVessel2D_CXXFLAGS = $(STD_FWK_INCLUDES)

PressureVessel2D_LDADD = $(STD_FWK_LIBS) 




VSNS2D_SOURCES = $(srcdir)/VSNS2D.cpp

VSNS2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

VSNS2D_CXXFLAGS = $(STD_FWK_INCLUDES)

VSNS2D_LDADD = $(STD_FWK_LIBS) 





include $(top_srcdir)/Makefile.common
