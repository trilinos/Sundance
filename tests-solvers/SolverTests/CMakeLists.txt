# CMake tests specification for Sundance/tests-solvers/SolverTests


INCLUDE(AddTestBatch)

SET(MPITests 
             PoissonTest 
             BlockTriangularTest 
             PoissonBoltzmannTest 
             PartitionedPoissonTest)

SET(SerialOnlyTests 
             UserDefPrecondTest)


SET(SolverParamPath ${PACKAGE_SOURCE_DIR}/etc/SolverParameters)

ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        DATA_FILES nox.xml poissonParams.xml
        DATA_SOURCE_DIR ${SolverParamPath}
        )

ADD_TEST_BATCH(SerialOnlyTests 
        COMM serial
        DATA_FILES userPrecParams.xml
        DATA_SOURCE_DIR ${SolverParamPath}
        )      



