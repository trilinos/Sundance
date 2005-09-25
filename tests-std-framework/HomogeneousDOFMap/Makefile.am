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

noinst_PROGRAMS = HomogDOFMapTest1D HomogDOFMapTest2D HomogDOFMapTest3D \
EqnSetTest1D


HomogDOFMapTest1D_SOURCES = $(srcdir)/HomogDOFMapTest1D.cpp

HomogDOFMapTest1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.la \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.la \
$(top_builddir)/src-core/libsundanceCore.la \
$(top_builddir)/src-utils/libsundanceUtils.la

HomogDOFMapTest1D_CXXFLAGS = $(STD_FWK_INCLUDES)

HomogDOFMapTest1D_LDADD = $(STD_FWK_LIBS) 



HomogDOFMapTest2D_SOURCES = $(srcdir)/HomogDOFMapTest2D.cpp

HomogDOFMapTest2D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.la \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.la \
$(top_builddir)/src-core/libsundanceCore.la \
$(top_builddir)/src-utils/libsundanceUtils.la

HomogDOFMapTest2D_CXXFLAGS = $(STD_FWK_INCLUDES)

HomogDOFMapTest2D_LDADD = $(STD_FWK_LIBS) 




HomogDOFMapTest3D_SOURCES = $(srcdir)/HomogDOFMapTest3D.cpp

HomogDOFMapTest3D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.la \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.la \
$(top_builddir)/src-core/libsundanceCore.la \
$(top_builddir)/src-utils/libsundanceUtils.la

HomogDOFMapTest3D_CXXFLAGS = $(STD_FWK_INCLUDES)

HomogDOFMapTest3D_LDADD = $(STD_FWK_LIBS) 





EqnSetTest1D_SOURCES = $(srcdir)/EqnSetTest1D.cpp

EqnSetTest1D_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.la \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.la \
$(top_builddir)/src-core/libsundanceCore.la \
$(top_builddir)/src-utils/libsundanceUtils.la

EqnSetTest1D_CXXFLAGS = $(STD_FWK_INCLUDES)

EqnSetTest1D_LDADD = $(STD_FWK_LIBS) 





include $(top_srcdir)/Makefile.common
