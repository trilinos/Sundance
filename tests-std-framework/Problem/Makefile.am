#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = NonlinAlgebraic




NonlinAlgebraic_SOURCES = $(srcdir)/NonlinAlgebraic.cpp

NonlinAlgebraic_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

NonlinAlgebraic_CXXFLAGS = $(STD_FWK_INCLUDES)

NonlinAlgebraic_LDADD = $(STD_FWK_LIBS) 







include $(top_srcdir)/Makefile.common
