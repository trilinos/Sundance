
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(TSF_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/../src-solvers/TSF_config.h
  )

#
# Cell Filtres files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/CellFilter)

APPEND_SET(HEADERS
  CellFilter/SundanceBinaryCellFilter.hpp
  CellFilter/SundanceBoundaryCellFilter.hpp
  CellFilter/SundanceBoundaryCellPredicate.hpp
  CellFilter/SundanceCellFilterBase.hpp
  CellFilter/SundanceCellFilter.hpp
  CellFilter/SundanceCellIterator.hpp
  CellFilter/SundanceCellPredicateBase.hpp
  CellFilter/SundanceCellPredicate.hpp
  CellFilter/SundanceCellSetBase.hpp
  CellFilter/SundanceCellSet.hpp
  CellFilter/SundanceCFMeshPair.hpp
  CellFilter/SundanceDimensionalCellFilter.hpp
  CellFilter/SundanceExplicitCellIterator.hpp
  CellFilter/SundanceExplicitCellSet.hpp
  CellFilter/SundanceImplicitCellSet.hpp
  CellFilter/SundanceLabelCellPredicate.hpp
  CellFilter/SundanceMaximalCellFilter.hpp
  CellFilter/SundanceNullCellFilter.hpp
  CellFilter/SundancePositionalCellPredicate.hpp
  CellFilter/SundanceSideCellFilter.hpp
  CellFilter/SundanceSubsetCellFilter.hpp
  )

APPEND_SET(SOURCES
  CellFilter/SundanceBinaryCellFilter.cpp
  CellFilter/SundanceBoundaryCellPredicate.cpp
  CellFilter/SundanceCellFilterBase.cpp
  CellFilter/SundanceCellFilter.cpp
  CellFilter/SundanceCellIterator.cpp
  CellFilter/SundanceCellPredicateBase.cpp
  CellFilter/SundanceCellPredicate.cpp
  CellFilter/SundanceCellSetBase.cpp
  CellFilter/SundanceCellSet.cpp
  CellFilter/SundanceCFMeshPair.cpp
  CellFilter/SundanceDimensionalCellFilter.cpp
  CellFilter/SundanceExplicitCellSet.cpp
  CellFilter/SundanceImplicitCellSet.cpp
  CellFilter/SundanceLabelCellPredicate.cpp
  CellFilter/SundanceMaximalCellFilter.cpp
  CellFilter/SundanceNullCellFilter.cpp
  CellFilter/SundancePositionalCellPredicate.cpp
  CellFilter/SundanceSideCellFilter.cpp
  CellFilter/SundanceSubsetCellFilter.cpp
  )

#
# Mesh Interface files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Interface)

APPEND_SET(HEADERS
  Interface/SundanceCellJacobianBatch.hpp
  Interface/SundanceCellType.hpp
  Interface/SundanceCellReordererBase.hpp
  Interface/SundanceCellReordererImplemBase.hpp
  Interface/SundanceCellReorderer.hpp
  Interface/SundanceGeomUtils.hpp
  Interface/SundanceIncrementallyCreatableMesh.hpp
  Interface/SundanceIdentityReorderer.hpp
  Interface/SundanceMaximalCofacetBatch.hpp
  Interface/SundanceMesh.hpp
  Interface/SundanceMeshBase.hpp
  Interface/SundanceMeshTypeBase.hpp
  Interface/SundanceMeshType.hpp
  )

APPEND_SET(SOURCES
  Interface/SundanceCellJacobianBatch.cpp
  Interface/SundanceCellType.cpp
  Interface/SundanceCellReordererImplemBase.cpp
  Interface/SundanceCellReorderer.cpp
  Interface/SundanceGeomUtils.cpp
  Interface/SundanceIdentityReorderer.cpp
  Interface/SundanceMaximalCofacetBatch.cpp
  Interface/SundanceMesh.cpp
  Interface/SundanceMeshBase.cpp
  Interface/SundanceMeshType.cpp
  )

#
# Concrete Mesh files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Concrete)

APPEND_SET(HEADERS
  Concrete/SundanceBasicSimplicialMesh.hpp
  Concrete/SundanceBasicSimplicialMeshType.hpp
  Concrete/SundanceBasicVertexView.hpp
  )

APPEND_SET(SOURCES
  Concrete/SundanceBasicSimplicialMesh.cpp
  Concrete/SundanceBasicVertexView.cpp
  )

#
# Mesh Source files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Sources)

APPEND_SET(HEADERS
  Sources/SundanceExodusMeshReader.hpp
  Sources/SundanceExodusNetCDFMeshReader.hpp
  Sources/SundanceMeshBuilder.hpp
  Sources/SundanceMeshReaderBase.hpp
  Sources/SundanceMeshSourceBase.hpp
  Sources/SundanceMeshSource.hpp
  Sources/SundancePartitionedLineMesher.hpp
  Sources/SundancePartitionedRectangleMesher.hpp
  Sources/SundanceTriangleMeshReader.hpp
  Sources/SundanceBamgMeshReader.hpp 
  )

APPEND_SET(SOURCES
  Sources/SundanceExodusMeshReader.cpp
  Sources/SundanceExodusNetCDFMeshReader.cpp
  Sources/SundanceMeshBuilder.cpp
  Sources/SundanceMeshReaderBase.cpp
  Sources/SundanceMeshSourceBase.cpp
  Sources/SundanceMeshSource.cpp
  Sources/SundancePartitionedLineMesher.cpp
  Sources/SundancePartitionedRectangleMesher.cpp
  Sources/SundanceTriangleMeshReader.cpp
  Sources/SundanceBamgMeshReader.cpp
  )

#
# Mesh Transformation files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Transformations)

APPEND_SET(HEADERS
  Transformations/SundanceExtrusionMeshTransformation.hpp
  Transformations/SundanceFileIOChacoPartitioner.hpp
  Transformations/SundanceMeshTransformationBase.hpp
  Transformations/SundanceMeshTransformation.hpp
  Transformations/SundanceRivaraEdge.hpp
  Transformations/SundanceRivaraElement.hpp
  Transformations/SundanceRivaraElementIterator.hpp
  Transformations/SundanceRivaraFace.hpp
  Transformations/SundanceRivaraMesh.hpp
  Transformations/SundanceRivaraNode.hpp
  Transformations/SundanceRivaraTreeNode.hpp
  Transformations/SundanceSerialPartitionerBase.hpp
  )

APPEND_SET(SOURCES
  Transformations/SundanceExtrusionMeshTransformation.cpp
  Transformations/SundanceFileIOChacoPartitioner.cpp
  Transformations/SundanceMeshTransformationBase.cpp
  Transformations/SundanceMeshTransformation.cpp
  Transformations/SundanceRivaraEdge.cpp
  Transformations/SundanceRivaraElement.cpp
  Transformations/SundanceRivaraElementIterator.cpp
  Transformations/SundanceRivaraMesh.cpp
  Transformations/SundanceRivaraNode.cpp
  Transformations/SundanceRivaraTreeNode.cpp
  Transformations/SundanceSerialPartitionerBase.cpp
  )

#
# Mesh Writers files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Writers)

APPEND_SET(HEADERS
  Writers/SundanceExodusWriter.hpp
  Writers/SundanceFieldBase.hpp
  Writers/SundanceFieldWriter.hpp
  Writers/SundanceFieldWriterBase.hpp
  Writers/SundanceMatlabWriter.hpp
  Writers/SundanceTriangleWriter.hpp
  Writers/SundanceVerboseFieldWriter.hpp
  Writers/SundanceVTKWriter.hpp
  )

APPEND_SET(SOURCES
  Writers/SundanceExodusWriter.cpp
  Writers/SundanceFieldBase.cpp
  Writers/SundanceFieldWriter.cpp
  Writers/SundanceFieldWriterBase.cpp
  Writers/SundanceMatlabWriter.cpp
  Writers/SundanceTriangleWriter.cpp
  Writers/SundanceVerboseFieldWriter.cpp
  Writers/SundanceVTKWriter.cpp
  )


IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  sundanceStdMesh
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
