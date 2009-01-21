
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

#
# C) Add the libraries, tests, and examples
#

#ADD_SUBDIRECTORY(src-solvers)
ADD_SUBDIRECTORY(src-utils)
ADD_SUBDIRECTORY(src-interop)
ADD_SUBDIRECTORY(src-core)
ADD_SUBDIRECTORY(src-std-mesh)
#ADD_SUBDIRECTORY(src-std-framework)

IF(Sundance_ENABLE_Python)
  ADD_SUBDIRECTORY(python)
ENDIF()

#PACKAGE_ADD_TEST_DIRECTORIES(test)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
