# CMake tests specification for Playa/tests/Vectors


INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
        BasicVectorTest
        SOURCES BasicVectorTest.cpp
        COMM serial mpi
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
        BlockVectorTest
        SOURCES BlockVectorTest.cpp
        COMM serial mpi
)
