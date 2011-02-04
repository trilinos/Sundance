
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(PlayaConfig.h)

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
  ${CMAKE_CURRENT_BINARY_DIR}/PlayaConfig.h
  )

#
# Utils files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  EpetraPlayaOperator.hpp
  NOX_Playa_Group.hpp
  NOX_Playa.hpp
  NOX_Playa_StatusTestBuilder.hpp
  NOX_Playa_Vector.hpp
  NOX_StatusTest_SafeCombo.hpp
  PlayaAmesosSolver.hpp
  PlayaAnasaziAdapter.hpp
  PlayaAnasaziEigensolverDecl.hpp
  PlayaAnasaziEigensolverImpl.hpp
  PlayaAztecSolver.hpp
  PlayaBelosSolver.hpp
  PlayaBICGSTABSolverDecl.hpp
  PlayaBICGSTABSolverImpl.hpp
  PlayaBlockIteratorDecl.hpp
  PlayaBlockIteratorImpl.hpp
  PlayaBlockTriangularSolverDecl.hpp
  PlayaBlockTriangularSolver.hpp
  PlayaBlockTriangularSolverImpl.hpp
  PlayaBlockVectorBaseDecl.hpp
  PlayaBlockVectorBaseImpl.hpp
  PlayaBlockVectorSpaceDecl.hpp
  PlayaBlockVectorSpaceImpl.hpp
  PlayaCollectivelyConfigurableMatrixFactory.hpp
  PlayaCompoundTester.hpp
  PlayaDataChunk.hpp
  PlayaDebug.hpp
  PlayaDefaultBlockVectorDecl.hpp
  PlayaDefaultBlockVectorImpl.hpp
  PlayaDefaultBlockVectorSpaceDecl.hpp
  PlayaDefaultBlockVectorSpaceImpl.hpp
  PlayaDefs.hpp
  PlayaDenseLUSolver.hpp
  PlayaDenseSerialMatrixFactory.hpp
  PlayaDenseSerialMatrix.hpp
  PlayaEigensolverBase.hpp
  PlayaEigensolver.hpp
  PlayaEpetraGhostImporter.hpp
  PlayaEpetraGhostView.hpp
  PlayaEpetraMatrixFactory.hpp
  PlayaEpetraMatrix.hpp
  PlayaEpetraMatrixMatrixProduct.hpp
  PlayaEpetraMatrixMatrixSum.hpp
  PlayaEpetraMatrixOps.hpp
  PlayaEpetraVector.hpp
  PlayaEpetraVectorSpace.hpp
  PlayaEpetraVectorType.hpp
  PlayaExceptions.hpp
  PlayaGeneralizedIndex.hpp
  PlayaGenericLeftPreconditioner.hpp
  PlayaGenericRightPreconditioner.hpp
  PlayaGhostImporter.hpp
  PlayaGhostView.hpp
  PlayaGlobalAnd.hpp
  PlayaHandleable.hpp
  PlayaHandle.hpp
  PlayaHeatOperator1D.hpp
  PlayaIfpackOperator.hpp
  PlayaILUFactorizableOp.hpp
  PlayaILUKPreconditionerFactory.hpp
  PlayaIncrementallyConfigurableMatrixFactory.hpp
  PlayaInverseOperatorDecl.hpp
  PlayaInverseOperatorImpl.hpp
  PlayaIterativeSolver.hpp
  PlayaKrylovSolver.hpp
  PlayaLeftPreconditioner.hpp
  PlayaLinearCombinationDecl.hpp
  PlayaLinearCombinationImpl.hpp
  PlayaLinearCombinationTester.hpp
  PlayaLinearOperatorBaseDecl.hpp
  PlayaLinearOperatorDecl.hpp
  PlayaLinearOperatorImpl.hpp
  PlayaLinearOpWithSpacesDecl.hpp
  PlayaLinearSolverBaseDecl.hpp
  PlayaLinearSolverBaseImpl.hpp
  PlayaLinearSolverBuilder.hpp
  PlayaLinearSolverDecl.hpp
  PlayaLinearSolverImpl.hpp
  PlayaLoadableMatrix.hpp
  PlayaLoadableVector.hpp
  PlayaMatrixFactory.hpp
  PlayaMatrixLaplacian1D.hpp
  PlayaMatrixMatrixTester.hpp
  PlayaMLOperator.hpp
  PlayaMultiVectorOperatorDecl.hpp
  PlayaMultiVectorOperatorImpl.hpp
  PlayaNonlinearSolverBase.hpp
  PlayaNOXSolver.hpp
  PlayaObjectWithVerbosity.hpp
  PlayaOperatorBuilder.hpp
  PlayaOut.hpp
  PlayaParameterListPreconditionerFactory.hpp
  PlayaPoissonBoltzmannJacobian.hpp
  PlayaPoissonBoltzmannOp.hpp
  PlayaPreconditionerBase.hpp
  PlayaPreconditionerFactoryBase.hpp
  PlayaPreconditionerFactory.hpp
  PlayaPreconditioner.hpp
  PlayaPrintable.hpp
  PlayaRand.hpp
  PlayaRandomBlockMatrixBuilderDecl.hpp
  PlayaRandomBlockMatrixBuilderImpl.hpp
  PlayaRandomSparseMatrixBuilderDecl.hpp
  PlayaRandomSparseMatrixBuilderImpl.hpp
  PlayaRowAccessibleOp.hpp
  PlayaSerialGhostImporter.hpp
  PlayaSerialGhostView.hpp
  PlayaSerialVector.hpp
  PlayaSerialVectorSpace.hpp
  PlayaSerialVectorType.hpp
  PlayaSimpleAddedOpDecl.hpp
  PlayaSimpleAddedOpImpl.hpp
  PlayaSimpleBlockOpDecl.hpp
  PlayaSimpleBlockOpImpl.hpp
  PlayaSimpleComposedOpDecl.hpp
  PlayaSimpleComposedOpImpl.hpp
  PlayaSimpleDiagonalOpDecl.hpp
  PlayaSimpleDiagonalOpImpl.hpp
  PlayaSimpleIdentityOpDecl.hpp
  PlayaSimpleIdentityOpImpl.hpp
  PlayaSimpleScaledOpDecl.hpp
  PlayaSimpleScaledOpImpl.hpp
  PlayaSimpleTransposedOpDecl.hpp
  PlayaSimpleTransposedOpImpl.hpp
  PlayaSimpleZeroOpDecl.hpp
  PlayaSimpleZeroOpImpl.hpp
  PlayaSingleChunkVector.hpp
  PlayaSolverState.hpp
  PlayaTabs.hpp
  PlayaTesterBase.hpp
  PlayaTestSpecifier.hpp
  PlayaVectorBaseDecl.hpp
  PlayaVectorDecl.hpp
  PlayaVectorFunctorsDecl.hpp
  PlayaVectorFunctorsImpl.hpp
  PlayaVectorImpl.hpp
  PlayaVectorSpaceBaseDecl.hpp
  PlayaVectorSpaceBaseImpl.hpp
  PlayaVectorSpaceDecl.hpp
  PlayaVectorSpaceImpl.hpp
  PlayaVectorTester.hpp
  PlayaVectorTypeBase.hpp
  PlayaVectorType.hpp

  )

APPEND_SET(SOURCES
  EpetraPlayaOperator.cpp
  NOX_Playa_Group.cpp
  NOX_Playa_StatusTestBuilder.cpp
  NOX_Playa_Vector.cpp
  NOX_StatusTest_SafeCombo.cpp
  PlayaAmesosSolver.cpp
  PlayaAnasaziEigensolver.cpp
  PlayaAztecSolver.cpp
  PlayaBelosSolver.cpp
  PlayaBICGSTABSolver.cpp
  PlayaBlockIterator.cpp
  PlayaBlockTriangularSolver.cpp
  PlayaBlockVectorBase.cpp
  PlayaBlockVectorSpace.cpp
  PlayaDebug.cpp
  PlayaDefaultBlockVector.cpp
  PlayaDefaultBlockVectorSpace.cpp
  PlayaDenseLUSolver.cpp
  PlayaDenseSerialMatrix.cpp
  PlayaDenseSerialMatrixFactory.cpp
  PlayaEpetraGhostImporter.cpp
  PlayaEpetraGhostView.cpp
  PlayaEpetraMatrix.cpp
  PlayaEpetraMatrixFactory.cpp
  PlayaEpetraMatrixMatrixProduct.cpp
  PlayaEpetraMatrixMatrixSum.cpp
  PlayaEpetraMatrixOps.cpp
  PlayaEpetraVector.cpp
  PlayaEpetraVectorSpace.cpp
  PlayaEpetraVectorType.cpp
  PlayaExceptions.cpp
  PlayaGlobalAnd.cpp
  PlayaHeatOperator1D.cpp
  PlayaIfpackOperator.cpp
  PlayaInverseOperator.cpp
  PlayaLinearCombination.cpp
  PlayaLinearOperator.cpp
  PlayaLinearSolverBase.cpp
  PlayaLinearSolverBuilder.cpp
  PlayaLinearSolver.cpp
  PlayaMatrixLaplacian1D.cpp
  PlayaMLOperator.cpp
  PlayaMultiVectorOperator.cpp
  PlayaNOXSolver.cpp
  PlayaOut.cpp
  PlayaParameterListPreconditionerFactory.cpp
  PlayaPoissonBoltzmannJacobian.cpp
  PlayaPoissonBoltzmannOp.cpp
  PlayaRand.cpp
  PlayaRandomBlockMatrixBuilder.cpp
  PlayaRandomSparseMatrixBuilder.cpp
  PlayaSerialGhostImporter.cpp
  PlayaSerialVector.cpp
  PlayaSerialVectorSpace.cpp
  PlayaSerialVectorType.cpp
  PlayaSimpleOps.cpp
  PlayaTabs.cpp
  PlayaVector.cpp
  PlayaVectorSpaceBase.cpp
  PlayaVectorSpace.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  playa
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
