# CMake tests specification 


INCLUDE(AddTestBatch)

# deactivated:       SensitivityTest 

SET(Tests 
        SetTest 
        ChainRule 
        SimpleEval 
        PolynomialTest          
        StdTransformation 
        VariationTest 
        testExpr
)
    
ADD_TEST_BATCH(Tests 
        COMM serial
        )



