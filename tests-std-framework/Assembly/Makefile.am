#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = AssemblerTest1D AssemblerTest2D


AssemblerTest1D_SOURCES = $(srcdir)/AssemblerTest1D.cpp

AssemblerTest1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

AssemblerTest1D_CXXFLAGS = $(STD_FWK_INCLUDES)

AssemblerTest1D_LDADD = $(STD_FWK_LIBS) 






AssemblerTest2D_SOURCES = $(srcdir)/AssemblerTest2D.cpp

AssemblerTest2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

AssemblerTest2D_CXXFLAGS = $(STD_FWK_INCLUDES)

AssemblerTest2D_LDADD = $(STD_FWK_LIBS) 











include $(top_srcdir)/Makefile.common
