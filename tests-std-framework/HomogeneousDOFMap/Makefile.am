#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = HomogDOFMapTest1D HomogDOFMapTest2D HomogDOFMapTest3D 


HomogDOFMapTest1D_SOURCES = $(srcdir)/HomogDOFMapTest1D.cpp

HomogDOFMapTest1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

HomogDOFMapTest1D_CXXFLAGS = $(STD_FWK_INCLUDES)

HomogDOFMapTest1D_LDADD = $(STD_FWK_LIBS) 



HomogDOFMapTest2D_SOURCES = $(srcdir)/HomogDOFMapTest2D.cpp

HomogDOFMapTest2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

HomogDOFMapTest2D_CXXFLAGS = $(STD_FWK_INCLUDES)

HomogDOFMapTest2D_LDADD = $(STD_FWK_LIBS) 




HomogDOFMapTest3D_SOURCES = $(srcdir)/HomogDOFMapTest3D.cpp

HomogDOFMapTest3D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

HomogDOFMapTest3D_CXXFLAGS = $(STD_FWK_INCLUDES)

HomogDOFMapTest3D_LDADD = $(STD_FWK_LIBS) 




include $(top_srcdir)/Makefile.common
