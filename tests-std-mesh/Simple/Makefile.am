#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = MeshTest1D MeshTest2D MeshTest3D

MeshTest1D_SOURCES = $(srcdir)/MeshTest1D.cpp

MeshTest1D_DEPENDENCIES = $(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-utils/libsundanceUtils.a

MeshTest1D_CXXFLAGS = $(STD_MESH_INCLUDES)

MeshTest1D_LDADD = $(STD_MESH_LIBS) 



MeshTest2D_SOURCES = $(srcdir)/MeshTest2D.cpp

MeshTest2D_DEPENDENCIES = $(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-utils/libsundanceUtils.a

MeshTest2D_CXXFLAGS = $(STD_MESH_INCLUDES)

MeshTest2D_LDADD = $(STD_MESH_LIBS) 




MeshTest3D_SOURCES = $(srcdir)/MeshTest3D.cpp

MeshTest3D_DEPENDENCIES = $(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-utils/libsundanceUtils.a

MeshTest3D_CXXFLAGS = $(STD_MESH_INCLUDES)

MeshTest3D_LDADD = $(STD_MESH_LIBS) 

include $(top_srcdir)/Makefile.common
