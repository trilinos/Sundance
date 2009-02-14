# CMake tests specification 


INCLUDE(AddTestBatch)

SET(Tests 
        SetTest 
        ChainRule 
        SimpleEval 
        SensitivityTest 
        PolynomialTest          
        StdTransformation 
        VariationTest 
        testExpr
)
    
ADD_TEST_BATCH(Tests 
        COMM serial
        )



