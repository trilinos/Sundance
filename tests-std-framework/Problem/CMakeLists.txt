

INCLUDE(AddTestBatch)

SET(MPITests 
             Poisson1D
   )


SET(SolverParamPath ${PACKAGE_SOURCE_DIR}/etc/SolverParameters)


IF (${HAVE_SUNDANCE_BROKEN_CODE})

ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        DATA_FILES nox.xml poissonParams.xml aztec-ifpack.xml aztec-ml.xml
        DATA_SOURCE_DIR ${SolverParamPath}
        )

ENDIF()