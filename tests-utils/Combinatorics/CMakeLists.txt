# CMake tests specification 


INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
        CombTest
        SOURCES CombTest.cpp
        COMM serial
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
        VertexOrderTest
        SOURCES VertexOrderTest.cpp
        COMM serial mpi
        NUM_MPI_PROCS 1
)
