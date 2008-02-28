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

noinst_PROGRAMS = ExodusNetCDFTest2D ExodusNetCDFTest3D ExodusTest3D 

bin_PROGRAMS = partitionExo

LIBDEPS =$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-interop/libsundanceInterop.a \
$(top_builddir)/src-utils/libsundanceUtils.a



ExodusNetCDFTest2D_SOURCES = $(srcdir)/ExodusNetCDFTest2D.cpp

ExodusNetCDFTest2D_DEPENDENCIES = $(LIBDEPS)

ExodusNetCDFTest2D_CXXFLAGS = $(STD_MESH_INCLUDES)

ExodusNetCDFTest2D_LDADD = $(STD_MESH_LIBS) 


ExodusNetCDFTest3D_SOURCES = $(srcdir)/ExodusNetCDFTest3D.cpp

ExodusNetCDFTest3D_DEPENDENCIES = $(LIBDEPS)

ExodusNetCDFTest3D_CXXFLAGS = $(STD_MESH_INCLUDES)

ExodusNetCDFTest3D_LDADD = $(STD_MESH_LIBS) 



ExodusTest3D_SOURCES = $(srcdir)/ExodusTest3D.cpp

ExodusTest3D_DEPENDENCIES = $(LIBDEPS)

ExodusTest3D_CXXFLAGS = $(STD_MESH_INCLUDES)

ExodusTest3D_LDADD = $(STD_MESH_LIBS) 



partitionExo_SOURCES = $(srcdir)/partitionExo.cpp

partitionExo_DEPENDENCIES = $(LIBDEPS)

partitionExo_CXXFLAGS = $(STD_MESH_INCLUDES)

partitionExo_LDADD = $(STD_MESH_LIBS) 


include $(top_srcdir)/Makefile.common
