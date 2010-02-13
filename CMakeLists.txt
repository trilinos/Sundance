
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Sundance)

#
# B) Set up package-specific options
#



ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Python
  HAVE_SUNDANCE_PYTHON
  "Enable python functionality."
  OFF )


ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_BROKEN_CODE
  HAVE_SUNDANCE_BROKEN_CODE
  "Enable broken code. For expert developer use only."
  OFF )


ADD_OPTION_AND_DEFINE(TPL_ENABLE_ExodusII
  HAVE_SUNDANCE_EXODUS
  "Enable ExodusII reader in Sundance (requires ExodusII TPL)."
  OFF )

ADD_OPTION_AND_DEFINE(Sundance_ENABLE_EXTERNAL_CHACO
  HAVE_SUNDANCE_CHACO
  "Enable file-based chaco partitioning in Sundance (requires chaco executable in path)."
  OFF )



# ========== Peano configuration ===============
IF(TPL_ENABLE_PEANO)
  ADD_OPTION_AND_DEFINE(TPL_ENABLE_PEANO
     HAVE_SUNDANCE_PEANO
    "Enable Peano Mesh in Sundance (requires Peano TPL)."
    OFF )  
  INCLUDE_DIRECTORIES(${PEANO_INCLUDE_DIRS}) 
  LINK_DIRECTORIES(${PEANO_LIBRARY_DIRS}/libpeano2D.a)
  LINK_DIRECTORIES(${PEANO_LIBRARY_DIRS}/libpeano3D.a)
ENDIF()
# ============ end of Peano configuration ================

IF(${Proj}_VERBOSE_CONFIGURE)
  PRINT_VAR(HAVE_SUNDANCE_EXODUS)
  PRINT_VAR(HAVE_SUNDANCE_PYTHON)
  PRINT_VAR(HAVE_SUNDANCE_CHACO)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src-utils)
PACKAGE_ADD_TEST_DIRECTORIES(tests-utils)

ADD_SUBDIRECTORY(src-solvers)
PACKAGE_ADD_TEST_DIRECTORIES(tests-solvers)

ADD_SUBDIRECTORY(src-interop)

ADD_SUBDIRECTORY(src-core)
PACKAGE_ADD_TEST_DIRECTORIES(tests-core)

ADD_SUBDIRECTORY(src-std-mesh)
PACKAGE_ADD_TEST_DIRECTORIES(tests-std-mesh)

ADD_SUBDIRECTORY(src-std-framework)
PACKAGE_ADD_TEST_DIRECTORIES(tests-std-framework)

IF(Sundance_ENABLE_Python)
  ADD_SUBDIRECTORY(python)
ENDIF()

PACKAGE_ADD_TEST_DIRECTORIES(examples-tutorial)


#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
