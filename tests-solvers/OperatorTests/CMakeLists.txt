# CMake tests specification for Sundance/tests-solvers/OperatorTests


INCLUDE(AddTestBatch)



SET(MPITests 
        LinearCombinationTest 
        CompoundOpTest
        MultiVectorOpTest
        BlockOpTest)

ADD_TEST_BATCH(
        MPITests
        COMM serial mpi
)
