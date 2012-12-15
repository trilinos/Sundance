#include <iostream>
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaMPISession.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaPCGSolver.hpp"
#include "PlayaICCPreconditionerFactory.hpp"
#include "PlayaMatrixMarketIO.hpp"

using std::cout;
using std::endl;
using std::setw;
using std::exception;
using namespace Playa;

int main(int argc, char** argv)
{
  int rtn = 0;
  try
    {
      MPISession::init(&argc, &argv);     
      Tabs::showDepth() = false;
      

      LinearOperator<double> A = readMatrixMarketFile("nut-stiffness-n-2-p-1.mtx");
      VectorSpace<double> spc = A.domain();

      Vector<double> ans = spc.createMember();
      ans.randomize();

      Vector<double> b = A*ans;

      ParameterList cgParams("Linear Solver");
      cgParams.set("Verbosity", 2);
      cgParams.set("Max Iterations", 100);
      cgParams.set("Tolerance", 1.0e-6);

      PreconditionerFactory<double> prec
      	= new ICCPreconditionerFactory<double>();
      LinearSolver<double> solver = new PCGSolver(cgParams, prec);

      Vector<double> x = b.copy();
      SolverState<double> state = solver.solve(A, b, x);

      if (state.finalState()==SolveConverged)
	{
	  Out::root() << "solve succeeded!" << endl;
	  Out::root() << "|error|=" << norm2(x-ans) << endl;
	}
      else
	{
	  Out::root() << "solve FAILED! Message=" << state.finalMsg() << endl;
	  rtn = -1;
	}
    }
  catch(std::exception& e)
    {
      cerr << "exception detected: " << e.what() << endl;
      rtn = -1;
    }

  MPISession::finalize();
  return rtn;
}
