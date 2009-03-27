
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
# Elements files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-solvers/Vectors)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-solvers/Operators)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-core/UserExprs)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-core/InternalExprs)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-core/FrameworkInterface)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Elements)

APPEND_SET(HEADERS
  Elements/SundanceBasisDOFTopologyBase.hpp
  Elements/SundanceBasisFamilyBase.hpp
  Elements/SundanceBasisFamily.hpp
  Elements/SundanceBasisReferenceEvaluationBase.hpp
  Elements/SundanceBubble.hpp
  Elements/SundanceClosedNewtonCotes.hpp
  Elements/SundanceClosedNewtonCotesType.hpp
  Elements/SundanceDiscreteFunction.hpp
  Elements/SundanceDiscreteFunctionData.hpp
  Elements/SundanceDiscreteSpace.hpp
  Elements/SundanceFIATExpansion.hpp
  Elements/SundanceFIATLagrange.hpp
  Elements/SundanceFIATP1NC.hpp
  Elements/SundanceFuncWithBasis.hpp
  Elements/SundanceGauss1D.hpp
  Elements/SundanceGaussianQuadrature.hpp
  Elements/SundanceGaussianQuadratureType.hpp
  Elements/SundanceLagrange.hpp
  Elements/SundanceNedelec.hpp
  Elements/SundanceP1NC.hpp
  Elements/SundanceQuadratureFamily.hpp
  Elements/SundanceQuadratureFamilyBase.hpp
  Elements/SundanceQuadratureType.hpp
  Elements/SundanceQuadratureTypeBase.hpp
  Elements/SundanceTestFunction.hpp
  Elements/SundanceTestFunctionData.hpp
  Elements/SundanceTetQuadrature.hpp
  Elements/SundanceTriangleQuadrature.hpp
  Elements/SundanceUnknownFunction.hpp
  Elements/SundanceUnknownFunctionData.hpp
  Elements/SundanceVectorBasisComponent.hpp
  )

APPEND_SET(SOURCES
  Elements/SundanceBasisFamily.cpp
  Elements/SundanceBasisFamilyBase.cpp
  Elements/SundanceBubble.cpp
  Elements/SundanceClosedNewtonCotes.cpp
  Elements/SundanceClosedNewtonCotesType.cpp
  Elements/SundanceDiscreteFunction.cpp
  Elements/SundanceDiscreteFunctionData.cpp
  Elements/SundanceDiscreteSpace.cpp
  Elements/SundanceFIATExpansion.cpp
  Elements/SundanceFIATLagrange.cpp
  Elements/SundanceFIATP1NC.cpp
  Elements/SundanceGauss1D.cpp
  Elements/SundanceGaussianQuadrature.cpp
  Elements/SundanceGaussianQuadratureType.cpp
  Elements/SundanceLagrange.cpp
  Elements/SundanceNedelec.cpp
  Elements/SundanceP1NC.cpp
  Elements/SundanceQuadratureFamily.cpp
  Elements/SundanceQuadratureFamilyBase.cpp
  Elements/SundanceQuadratureType.cpp
  Elements/SundanceTestFunction.cpp
  Elements/SundanceTestFunctionData.cpp
  Elements/SundanceTetQuadrature.cpp
  Elements/SundanceTriangleQuadrature.cpp
  Elements/SundanceUnknownFunction.cpp
  Elements/SundanceUnknownFunctionData.cpp
  Elements/SundanceVectorBasisComponent.cpp
  )

#
# DOF Mapping files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/DOFMapping)

APPEND_SET(HEADERS
  DOFMapping/SundanceBlock.hpp
  DOFMapping/SundanceDOFMapBase.hpp
  DOFMapping/SundanceDOFMapBuilder.hpp
  DOFMapping/SundanceInhomogeneousNodalDOFMap.hpp
  DOFMapping/SundanceMapStructure.hpp
  DOFMapping/SundanceMixedDOFMap.hpp
  DOFMapping/SundanceNodalDOFMap.hpp
  DOFMapping/SundancePartialElementDOFMap.hpp
  DOFMapping/SundanceSpatiallyHomogeneousDOFMapBase.hpp
  )

APPEND_SET(SOURCES
  DOFMapping/SundanceDOFMapBase.cpp
  DOFMapping/SundanceDOFMapBuilder.cpp
  DOFMapping/SundanceInhomogeneousNodalDOFMap.cpp
  DOFMapping/SundanceMapStructure.cpp
  DOFMapping/SundanceMixedDOFMap.cpp
  DOFMapping/SundanceNodalDOFMap.cpp
  DOFMapping/SundancePartialElementDOFMap.cpp
  DOFMapping/SundanceSpatiallyHomogeneousDOFMapBase.cpp
  )

#
# Assembly files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Assembly)

APPEND_SET(HEADERS
  Assembly/SundanceAssembler.hpp
  Assembly/SundanceAssemblyKernelBase.hpp
  Assembly/SundanceElementIntegral.hpp
  Assembly/SundanceFunctionalAssemblyKernel.hpp
  Assembly/SundanceFunctionalGradientAssemblyKernel.hpp
  Assembly/SundanceFunctionalEvaluator.hpp
  Assembly/SundanceGrouperBase.hpp
  Assembly/SundanceIntegralGroup.hpp
  Assembly/SundanceIntegrationCellSpecifier.hpp
  Assembly/SundanceLocalDOFMap.hpp
  Assembly/SundanceLocalMatrixContainer.hpp
  Assembly/SundanceMapBundle.hpp
  Assembly/SundanceMatrixVectorAssemblyKernel.hpp
  Assembly/SundanceQuadratureEvalMediator.hpp
  Assembly/SundanceQuadratureIntegral.hpp
  Assembly/SundanceRefIntegral.hpp
  Assembly/SundanceStdFwkEvalMediator.hpp
  Assembly/SundanceTrivialGrouper.hpp
  Assembly/SundanceVectorAssemblyKernel.hpp
  Assembly/SundanceVectorFillingAssemblyKernel.hpp
  )

APPEND_SET(SOURCES
  Assembly/SundanceAssembler.cpp
  Assembly/SundanceElementIntegral.cpp
  Assembly/SundanceFunctionalAssemblyKernel.cpp
  Assembly/SundanceFunctionalEvaluator.cpp
  Assembly/SundanceGrouperBase.cpp
  Assembly/SundanceIntegralGroup.cpp
  Assembly/SundanceLocalDOFMap.cpp
  Assembly/SundanceLocalMatrixContainer.cpp
  Assembly/SundanceMapBundle.cpp
  Assembly/SundanceMatrixVectorAssemblyKernel.cpp
  Assembly/SundanceQuadratureEvalMediator.cpp
  Assembly/SundanceQuadratureIntegral.cpp
  Assembly/SundanceRefIntegral.cpp
  Assembly/SundanceStdFwkEvalMediator.cpp
  Assembly/SundanceTrivialGrouper.cpp
  Assembly/SundanceVectorAssemblyKernel.cpp
  Assembly/SundanceVectorFillingAssemblyKernel.cpp
  )

#
# Problem files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Problem)

APPEND_SET(HEADERS
  Problem/Sundance.hpp
  Problem/SundanceAToCDensitySampler.hpp
  Problem/SundanceAToCPointLocator.hpp
  Problem/SundanceCToAInterpolator.hpp
  Problem/SundanceCoordinateSystem.hpp
  Problem/SundanceEigensolution.hpp
  Problem/SundanceExprFieldWrapper.hpp
  Problem/SundanceFunctional.hpp
  Problem/SundanceL2Projector.hpp
  Problem/SundanceLinearEigenproblem.hpp
  Problem/SundanceLinearProblem.hpp
  Problem/SundanceModelEvaluatorBase.hpp
  Problem/SundanceNLPModelEvaluator.hpp
  Problem/SundanceRivaraDriver.hpp
  Problem/SundanceNonlinearProblem.hpp 
  )

APPEND_SET(SOURCES
  Problem/Sundance.cpp
  Problem/SundanceAToCDensitySampler.cpp
  Problem/SundanceAToCPointLocator.cpp
  Problem/SundanceCToAInterpolator.cpp
  Problem/SundanceCoordinateSystem.cpp
  Problem/SundanceExprFieldWrapper.cpp
  Problem/SundanceFunctional.cpp
  Problem/SundanceL2Projector.cpp
  Problem/SundanceLinearEigenproblem.cpp
  Problem/SundanceLinearProblem.cpp
  Problem/SundanceModelEvaluatorBase.cpp
  Problem/SundanceNLPModelEvaluator.cpp
  Problem/SundanceRivaraDriver.cpp
  Problem/SundanceNonlinearProblem.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  sundanceStdFwk
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
