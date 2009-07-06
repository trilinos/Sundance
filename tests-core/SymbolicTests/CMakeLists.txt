# CMake tests specification 


INCLUDE(AddTestBatch)

# deactivated:       SensitivityTest 

SET(Tests 
        DerivSpecTest
        PreprocTest
        ChainRule 
        SimpleEval 
        SensitivityTest 
        SensEqnTest 
        TestEval 
        PolynomialTest          
        StdTransformation 
        VariationTest 
        testExpr
)
    
ADD_TEST_BATCH(Tests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )



