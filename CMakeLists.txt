
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(WriteExportMakefile)

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



PRINT_VAR(HAVE_SUNDANCE_EXODUS)


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src-solvers)
ADD_SUBDIRECTORY(src-utils)
ADD_SUBDIRECTORY(src-interop)
ADD_SUBDIRECTORY(src-core)
ADD_SUBDIRECTORY(src-std-mesh)
ADD_SUBDIRECTORY(src-std-framework)

IF(Sundance_ENABLE_Python)
  ADD_SUBDIRECTORY(python)
ENDIF()

PACKAGE_ADD_TEST_DIRECTORIES(tests-utils)
PACKAGE_ADD_TEST_DIRECTORIES(tests-solvers)
PACKAGE_ADD_TEST_DIRECTORIES(tests-core)
PACKAGE_ADD_TEST_DIRECTORIES(tests-std-mesh)
PACKAGE_ADD_TEST_DIRECTORIES(tests-std-framework)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)



WRITE_EXPORT_MAKEFILE(Trilinos Makefile.export.cmake)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
