# @HEADER
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


SUNDANCE_MESHINTERFACE_H = \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellJacobianBatch.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellType.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReordererBase.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReordererImplemBase.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReorderer.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCreatableMesh.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceIdentityReorderer.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMesh.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshBase.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshTypeBase.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshType.hpp \
$(top_srcdir)/src-std-mesh/Interface/SundancePoint.hpp


SUNDANCE_MESHINTERFACE = \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellJacobianBatch.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellType.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReordererImplemBase.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceCellReorderer.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceIdentityReorderer.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMesh.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshBase.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundanceMeshType.cpp \
$(top_srcdir)/src-std-mesh/Interface/SundancePoint.cpp

SUNDANCE_CONCRETEMESHES_H = \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicSimplicialMesh.hpp \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicSimplicialMeshType.hpp \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicVertexSet.hpp 

SUNDANCE_CONCRETEMESHES = \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicSimplicialMesh.cpp \
$(top_srcdir)/src-std-mesh/Concrete/SundanceBasicVertexSet.cpp 

SUNDANCE_MESHSOURCES_H = \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshReaderBase.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshSourceBase.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshSource.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundancePartitionedLineMesher.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundancePartitionedRectangleMesher.hpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceTriangleMeshReader.hpp 

SUNDANCE_MESHSOURCES = \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshReaderBase.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshSourceBase.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceMeshSource.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundancePartitionedLineMesher.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundancePartitionedRectangleMesher.cpp \
$(top_srcdir)/src-std-mesh/Sources/SundanceTriangleMeshReader.cpp 

SUNDANCE_MESHTRANSFORMATIONS_H = \
$(top_srcdir)/src-std-mesh/Transformations/SundanceExtrusionMeshTransformation.hpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceMeshTransformationBase.hpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceMeshTransformation.hpp 

SUNDANCE_MESHTRANSFORMATIONS = \
$(top_srcdir)/src-std-mesh/Transformations/SundanceExtrusionMeshTransformation.cpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceMeshTransformationBase.cpp \
$(top_srcdir)/src-std-mesh/Transformations/SundanceMeshTransformation.cpp 

SUNDANCE_MESHWRITERS_H = \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldBase.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldWriter.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldWriterBase.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceTriangleWriter.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceVerboseFieldWriter.hpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceVTKWriter.hpp 

SUNDANCE_MESHWRITERS = \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldWriter.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceFieldWriterBase.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceTriangleWriter.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceVerboseFieldWriter.cpp \
$(top_srcdir)/src-std-mesh/Writers/SundanceVTKWriter.cpp 

STD_MESH_SOURCES = \
$(SUNDANCE_MESHINTERFACE) \
$(SUNDANCE_CONCRETEMESHES) \
$(SUNDANCE_MESHTRANSFORMATIONS) \
$(SUNDANCE_MESHSOURCES) \
$(SUNDANCE_MESHWRITERS) 

STD_MESH_SOURCES_H = \
$(SUNDANCE_MESHINTERFACE_H) \
$(SUNDANCE_CONCRETEMESHES_H) \
$(SUNDANCE_MESHTRANSFORMATIONS_H) \
$(SUNDANCE_MESHSOURCES_H)  \
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

AM_CPPFLAGS = $(STD_MESH_INCLUDES)

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






