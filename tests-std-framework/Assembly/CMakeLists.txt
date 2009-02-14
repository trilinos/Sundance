# CMake tests specification for Sundance/tests-solvers/SolverTests


INCLUDE(AddTestBatch)

SET(SerialTests BasisCheck TransformedIntegral2D  QuadratureTest)


ADD_TEST_BATCH(SerialTests 
        COMM serial
        )      



