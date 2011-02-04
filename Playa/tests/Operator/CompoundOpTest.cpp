//@HEADER

//@HEADER


#include <cstdlib>
#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaRandomSparseMatrixBuilderDecl.hpp"
#include "PlayaCompoundTester.hpp"
#include "PlayaMatrixMatrixTester.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#endif

STREAM_OUT(Vector<double>)

using namespace Playa;
using namespace PlayaOps;

int main(int argc, char *argv[]) 
{
  int stat = 0;
  try
    {
      GlobalMPISession session(&argc, &argv);
 

      VectorType<double> type = new EpetraVectorType();

      int nLocalRows = 4;
      
      double onProcDensity = 0.5;
      double offProcDensity = 0.1;

      RandomSparseMatrixBuilder<double> ABuilder(nLocalRows, nLocalRows, 
        onProcDensity, offProcDensity, type);
      RandomSparseMatrixBuilder<double> BBuilder(nLocalRows, nLocalRows, 
        onProcDensity, offProcDensity, type);

      /* Build some rectangular matrices to test products */
      RandomSparseMatrixBuilder<double> CBuilder(2*nLocalRows, nLocalRows, 
        onProcDensity, offProcDensity, type);
      RandomSparseMatrixBuilder<double> DBuilder(3*nLocalRows, 2*nLocalRows, 
        onProcDensity, offProcDensity, type);

      LinearOperator<double> A = ABuilder.getOp();
      LinearOperator<double> B = BBuilder.getOp();

      LinearOperator<double> C = CBuilder.getOp();
      LinearOperator<double> D = DBuilder.getOp();

      Out::root() << "A = " << std::endl;
      Out::os() << A << std::endl;
      Out::root() << "B = " << std::endl;
      Out::os() << B << std::endl;

      Out::root() << "C = " << std::endl;
      Out::os() << C << std::endl;
      Out::root() << "D = " << std::endl;
      Out::os() << D << std::endl;
      
      CompoundTester<double> tester(A, B, 
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

      bool allPass =  tester.runAllTests();

      Out::root() << std::endl << std::endl 
                  << "testing multiplication of square matrices " 
                  << std::endl << std::endl;

      MatrixMatrixTester<double> mmTester(A, B, 
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

      allPass = mmTester.runAllTests() && allPass;

      Out::root() << std::endl << std::endl 
                  << "testing multiplication of rectangular matrices " 
                  << std::endl << std::endl;

      MatrixMatrixTester<double> rectMMTester(C, D, 
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

      allPass = rectMMTester.runAllTests() && allPass;

     if (!allPass) stat = -1;
    }
  catch(std::exception& e)
    {
      stat = 0;
      std::cerr << "Caught exception: " << e.what() << std::endl;
    }
  return stat;
}



