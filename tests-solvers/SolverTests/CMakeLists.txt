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
        )

ADD_TEST_BATCH(SerialOnlyTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )      


PACKAGE_COPY_FILES_TO_BINARY_DIR(SolverTestCopyFiles
  DEST_FILES 
                   amesos.xml 
                   anasazi-ml.xml 
                   aztec-ifpack.xml 
                   aztec-ml.xml 
                   aztec-native.xml
                   aztec.xml
                   bicgstab.xml
                   gmres.xml
                   nox.xml            
                   poissonParams.xml
                   userPrecParams.xml
  SOURCE_DIR ${SolverParamPath}
  EXEDEPS ${MPITests} ${SerialTests}
  )
