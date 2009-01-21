
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(TSF_config.h)
PACKAGE_CONFIGURE_FILE(SundanceConfig.h)

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
  ${CMAKE_CURRENT_BINARY_DIR}/SundanceConfig.h
  )

#
# Utils files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities)

APPEND_SET(HEADERS
  Utilities/SundanceADReal.hpp
  Utilities/SundanceArrayOfTuples.hpp
  Utilities/SundanceCollectiveExceptionCheck.hpp
  Utilities/SundanceCombinatorialUtils.hpp
  Utilities/SundanceDebug.hpp
  Utilities/SundanceDefs.hpp
  Utilities/SundanceExceptions.hpp
  Utilities/SundanceFunctorDomain.hpp
  Utilities/SundanceIntHashSet.hpp
  Utilities/SundanceMap.hpp
  Utilities/SundanceMultiSet.hpp
  Utilities/SundanceNoncopyable.hpp
  Utilities/SundanceObjectWithInstanceID.hpp
  Utilities/SundanceOrderedHandle.hpp
  Utilities/SundanceOrderedTuple.hpp
  Utilities/SundanceOut.hpp
  Utilities/SundancePoint.hpp
  Utilities/SundanceSet.hpp
  Utilities/SundanceStdMathFunctors.hpp
  Utilities/SundanceTabs.hpp
  Utilities/SundanceUnaryFunctor.hpp
  Utilities/SundanceVersionString.hpp
  )

APPEND_SET(SOURCES
  Utilities/SundanceADReal.cpp
  Utilities/SundanceCollectiveExceptionCheck.cpp
  Utilities/SundanceCombinatorialUtils.cpp
  Utilities/SundanceExceptions.cpp
  Utilities/SundanceFunctorDomain.cpp
  Utilities/SundanceIntHashSet.cpp
  Utilities/SundancePoint.cpp
  Utilities/SundanceStdMathFunctors.cpp
  Utilities/SundanceTabs.cpp
  Utilities/SundanceUnaryFunctor.cpp
  )


IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  sundanceUtils
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
