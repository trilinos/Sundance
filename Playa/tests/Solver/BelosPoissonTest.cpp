//@HEADER@

//@HEADER@

#include "BelosBlockGmresSolMgr.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaPreconditioner.hpp"
#include "PlayaPreconditionerFactory.hpp"
#include "PlayaParameterListPreconditionerFactory.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaLoadableVector.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaMatrixLaplacian1D.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaInverseOperatorDecl.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaInverseOperatorImpl.hpp"

#endif

using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;


int main(int argc, char *argv[]) 
{
  typedef Teuchos::ScalarTraits<double> ST;

  try
  {
    GlobalMPISession session(&argc, &argv);

    MPIComm::world().synchronize();

    VectorType<double> type = new EpetraVectorType();


    ParameterXMLFileReader reader("belos-ml.xml");
    ParameterList solverParams = reader.getParameters();

    /* create the range space  */
    int nLocalRows = 20;

    MatrixLaplacian1D builder(nLocalRows, type);

    LinearOperator<double> A = builder.getOp();

    Vector<double> x = A.domain().createMember();
    int myRank = MPIComm::world().getRank();
    int nProcs = MPIComm::world().getNProc();
      

    x.randomize();

    if (myRank==0) loadable(x)->setElement(0, 0);
    if (myRank==nProcs-1) loadable(x)->setElement(nProcs * nLocalRows - 1, 0.0);

    cout << "input is " << std::endl;
    x.print(cout);
    Vector<double> b = A*x;

    cout << "rhs is " << std::endl;
    b.print(cout);




    LinearSolver<double> solver = LinearSolverBuilder::createSolver(solverParams);
    LinearOperator<double> AInv = inverse(A, solver);
    Vector<double> ans = AInv * b;

    cout << "answer is " << std::endl;
    ans.print(cout);
      
    double err = (x-ans).norm2()/((double) nProcs * nLocalRows);
    cout << "error norm = " << err << std::endl;

    if (err <= 1.0e-8)
    {
      cout << "Belos poisson solve test PASSED" << std::endl;
      return 0;
    }
    else
    {
      cout << "Belos poisson solve test FAILED" << std::endl;
      return 1;
    }
  }
  catch(std::exception& e)
  {
    cout << "Caught exception: " << e.what() << std::endl;
    return -1;
  }
}

