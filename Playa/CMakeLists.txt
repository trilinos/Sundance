
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

# dummy change to test secondary git repo
# another dummy change

#
# A) Define the package
#

#PACKAGE(Playa  DISABLE_CIRCULAR_REF_DETECTION_FAILURE)

#
# B) Set up package-specific options
#

# ------- nothing to do here

#
# C) Add the libraries, tests, and examples
#

MESSAGE("adding Playa/src")
ADD_SUBDIRECTORY(src)
MESSAGE("adding Playa/tests")
PACKAGE_ADD_TEST_DIRECTORIES(tests)


#
# Exclude files for source package.
#

#PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

#PACKAGE_POSTPROCESS()
