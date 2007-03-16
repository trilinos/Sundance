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

noinst_PROGRAMS = MeshTest1D MeshTest2D MeshTest3D

LIBDEPS =$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-interop/libsundanceInterop.a \
$(top_builddir)/src-utils/libsundanceUtils.a

MeshTest1D_SOURCES = $(srcdir)/MeshTest1D.cpp

MeshTest1D_DEPENDENCIES = $(LIBDEPS)

MeshTest1D_CXXFLAGS = $(STD_MESH_INCLUDES)

MeshTest1D_LDADD = $(STD_MESH_LIBS) 



MeshTest2D_SOURCES = $(srcdir)/MeshTest2D.cpp

MeshTest2D_DEPENDENCIES = $(LIBDEPS)

MeshTest2D_CXXFLAGS = $(STD_MESH_INCLUDES)

MeshTest2D_LDADD = $(STD_MESH_LIBS) 




MeshTest3D_SOURCES = $(srcdir)/MeshTest3D.cpp

MeshTest3D_DEPENDENCIES = $(LIBDEPS)

MeshTest3D_CXXFLAGS = $(STD_MESH_INCLUDES)

MeshTest3D_LDADD = $(STD_MESH_LIBS) 

include $(top_srcdir)/Makefile.common
