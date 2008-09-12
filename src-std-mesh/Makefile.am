# @HEADER
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
# @HEADER

#np# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
#np# There is, at present, no distinction between most of the categories, 
#np# other than those which start with EXTRA_.  All categories which do not
#np# start with EXTRA_ are included in all builds.  Categories which start
#np# with EXTRA_ are included in some builds (typically controlled by 
#np# --with-xyz or --without-xyz)
#
#np# Make sure to list all source files in one of the following categories.

SUNDANCE_CELLFILTERS_H = \
$(top_srcdir)/src-std-mesh/CellFilter/SundanceBinaryCellFilter.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceBoundaryCellFilter.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceBoundaryCellPredicate.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellFilterBase.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellFilter.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellIterator.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellPredicateBase.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellPredicate.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellSetBase.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellSet.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCFMeshPair.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceDimensionalCellFilter.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceExplicitCellIterator.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceExplicitCellSet.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceImplicitCellSet.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceLabelCellPredicate.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceMaximalCellFilter.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceNullCellFilter.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundancePositionalCellPredicate.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceSideCellFilter.hpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceSubsetCellFilter.hpp

SUNDANCE_CELLFILTERS = \
$(top_srcdir)/src-std-mesh/CellFilter/SundanceBinaryCellFilter.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceBoundaryCellPredicate.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellFilterBase.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellFilter.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellIterator.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellPredicateBase.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellPredicate.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellSetBase.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCellSet.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceCFMeshPair.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceDimensionalCellFilter.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceExplicitCellSet.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceImplicitCellSet.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceLabelCellPredicate.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceMaximalCellFilter.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundancePositionalCellPredicate.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceSideCellFilter.cpp\
$(top_srcdir)/src-std-mesh/CellFilter/SundanceSubsetCellFilter.cpp

SUNDANCE_MESHINTERFACE_H = \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellJacobianBatch.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellType.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReordererBase.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReordererImplemBase.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReorderer.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceGeomUtils.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceIncrementallyCreatableMesh.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceIdentityReorderer.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMesh.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshBase.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshTypeBase.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshType.hpp 



SUNDANCE_MESHINTERFACE = \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellJacobianBatch.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellType.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReordererImplemBase.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReorderer.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceGeomUtils.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceIdentityReorderer.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMesh.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshBase.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshType.cpp 

SUNDANCE_CONCRETEMESHES_H = \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicSimplicialMesh.hpp \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicSimplicialMeshType.hpp \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicVertexView.hpp 

SUNDANCE_CONCRETEMESHES = \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicSimplicialMesh.cpp \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicVertexView.cpp 

SUNDANCE_MESHSOURCES_H = \
$(top_srcdir)/src-std-mesh/Sources/SundanceExodusMeshReader.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceExodusNetCDFMeshReader.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshBuilder.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshReaderBase.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshSourceBase.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshSource.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundancePartitionedLineMesher.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundancePartitionedRectangleMesher.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceTriangleMeshReader.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceBamgMeshReader.hpp 

SUNDANCE_MESHSOURCES = \
$(top_srcdir)/src-std-mesh/Sources/SundanceExodusMeshReader.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceExodusNetCDFMeshReader.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshBuilder.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshReaderBase.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshSourceBase.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshSource.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundancePartitionedLineMesher.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundancePartitionedRectangleMesher.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceTriangleMeshReader.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceBamgMeshReader.cpp 

SUNDANCE_MESHTRANSFORMATIONS_H = \
$(top_srcdir)/src-std-mesh/Transformations/SundanceExtrusionMeshTransformation.hpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceFileIOChacoPartitioner.hpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceMeshTransformationBase.hpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceMeshTransformation.hpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceSerialPartitionerBase.hpp 

SUNDANCE_MESHTRANSFORMATIONS = \
$(top_srcdir)/src-std-mesh/Transformations/SundanceExtrusionMeshTransformation.cpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceFileIOChacoPartitioner.cpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceMeshTransformationBase.cpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceMeshTransformation.cpp  \
$(top_srcdir)/src-std-mesh/Transformations/SundanceSerialPartitionerBase.cpp 

SUNDANCE_MESHWRITERS_H = \
$(top_srcdir)/src-std-mesh/Writers/SundanceExodusWriter.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldBase.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldWriter.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldWriterBase.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceMatlabWriter.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceTriangleWriter.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceVerboseFieldWriter.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceVTKWriter.hpp 

SUNDANCE_MESHWRITERS = \
$(top_srcdir)/src-std-mesh/Writers/SundanceExodusWriter.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldWriter.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldWriterBase.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceMatlabWriter.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceTriangleWriter.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceVerboseFieldWriter.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceVTKWriter.cpp 

STD_MESH_SOURCES = \
$(SUNDANCE_MESHINTERFACE) \
$(SUNDANCE_CONCRETEMESHES) \
$(SUNDANCE_MESHTRANSFORMATIONS) \
$(SUNDANCE_MESHSOURCES) \
$(SUNDANCE_CELLFILTERS) \
$(SUNDANCE_MESHWRITERS) 

STD_MESH_SOURCES_H = \
$(SUNDANCE_MESHINTERFACE_H) \
$(SUNDANCE_CONCRETEMESHES_H) \
$(SUNDANCE_MESHTRANSFORMATIONS_H) \
$(SUNDANCE_MESHSOURCES_H)  \
$(SUNDANCE_CELLFILTERS_H) \
$(SUNDANCE_MESHWRITERS_H) 

include $(top_srcdir)/Makefile.common

#
#np#  Files which are included in select builds are typically have
#np#  names that start with EXTRA
#
#np#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#np#  set up by TAC_ARC_CONFIG_MPI 
#
EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif

AM_CPPFLAGS = $(STD_MESH_INCLUDES) $(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# Sundance library specifications
# ------------------------------------------------------------------------

STD_MESH_LIB = libsundanceStdMesh.a


#np# replace new_package with the name of the package being autotool'ed here
libsundanceStdMesh_a_SOURCES = $(STD_MESH_SOURCES)

#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libsundanceStdMesh_a_SOURCES = \
$(EXTRA_MPI) \
$(EXTRA_MPI_H) 

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libsundanceStdMesh_a_AR = $(ALTERNATE_AR)
else

libsundanceStdMesh_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libsundanceStdMesh_a_LIBADD = $(SUNDANCE_TEMPLATE_OBJS) $(XTRALDADD) 

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(STD_MESH_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(STD_MESH_SOURCES_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SundanceConfig.h was installed when the utiities were built, so skip it.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = 

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






