#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = NonlinAlgebraic Poisson1D Poisson2D




NonlinAlgebraic_SOURCES = $(srcdir)/NonlinAlgebraic.cpp

NonlinAlgebraic_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

NonlinAlgebraic_CXXFLAGS = $(STD_FWK_INCLUDES)

NonlinAlgebraic_LDADD = $(STD_FWK_LIBS) 





Poisson1D_SOURCES = $(srcdir)/Poisson1D.cpp

Poisson1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

Poisson1D_CXXFLAGS = $(STD_FWK_INCLUDES)

Poisson1D_LDADD = $(STD_FWK_LIBS) 





Poisson2D_SOURCES = $(srcdir)/Poisson2D.cpp

Poisson2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

Poisson2D_CXXFLAGS = $(STD_FWK_INCLUDES)

Poisson2D_LDADD = $(STD_FWK_LIBS) 







include $(top_srcdir)/Makefile.common
