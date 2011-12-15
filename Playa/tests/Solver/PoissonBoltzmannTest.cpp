//@HEADER@

//@HEADER@

#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaPoissonBoltzmannOp.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaNOXSolver.hpp"
#include "PlayaLinearCombinationImpl.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"

#endif

using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;




int main(int argc, char *argv[]) 
{
  try
  {
    GlobalMPISession session(&argc, &argv);


    MPIComm::world().synchronize();

    /* create the nonlinear operator */
    VectorType<double> type = new EpetraVectorType();
    int nProc = MPIComm::world().getNProc();
    int nLocalRows = 128/nProc;
    PoissonBoltzmannOp* prob = new PoissonBoltzmannOp(nLocalRows, type);
    NonlinearOperator<double> F = prob;

    /* create the nox solver */
    ParameterXMLFileReader reader("nox.xml");
    ParameterList noxParams = reader.getParameters();

    Out::root() << "solver params = " << noxParams << std::endl;

    NOXSolver solver(noxParams);

    Vector<double> soln;
    SolverState<double> stat = solver.solve(F, soln);
    TEUCHOS_TEST_FOR_EXCEPTION(stat.finalState() != SolveConverged,
      runtime_error, "solve failed");

    Out::root() << "numerical solution = " << std::endl;
    Out::os() << soln << std::endl;

    Vector<double> exact = prob->exactSoln();

    Out::root() << "exact solution = " << std::endl;
    Out::os() << exact << std::endl;

//bvbw reddish port hack
    double temp_val = nLocalRows*nProc;
    double err = (exact-soln).norm2()/sqrt(temp_val);
    Out::root() << "error norm = " << err << std::endl;
      

    double tol = 1.0e-6;
    if (err > tol)
    {
      Out::root() << "NOX Poisson-Boltzmann test FAILED" << std::endl;
      return 1;
    }
    else
    {
      Out::root() << "NOX Poisson-Boltzmann test PASSED" << std::endl;
      return 0;
    }
  }
  catch(std::exception& e)
  {
    Out::root() << "Caught exception: " << e.what() << std::endl;
    return -1;
  }
}

