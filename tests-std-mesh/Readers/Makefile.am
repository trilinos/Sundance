#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = ExodusNetCDFTest2D ExodusNetCDFTest3D


ExodusNetCDFTest2D_SOURCES = $(srcdir)/ExodusNetCDFTest2D.cpp

ExodusNetCDFTest2D_DEPENDENCIES = $(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-utils/libsundanceUtils.a

ExodusNetCDFTest2D_CXXFLAGS = $(STD_MESH_INCLUDES)

ExodusNetCDFTest2D_LDADD = $(STD_MESH_LIBS) 


ExodusNetCDFTest3D_SOURCES = $(srcdir)/ExodusNetCDFTest3D.cpp

ExodusNetCDFTest3D_DEPENDENCIES = $(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-utils/libsundanceUtils.a

ExodusNetCDFTest3D_CXXFLAGS = $(STD_MESH_INCLUDES)

ExodusNetCDFTest3D_LDADD = $(STD_MESH_LIBS) 


include $(top_srcdir)/Makefile.common
