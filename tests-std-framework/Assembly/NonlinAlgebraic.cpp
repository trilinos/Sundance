#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceBruteForceEvaluator.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceBruteForceEvaluator.hpp"
#include "SundanceBasicInserter.hpp"
#include "SundanceIntegrator.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "TSFVectorType.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFBICGSTABSolver.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFLinearCombination.hpp"

using namespace TSFExtended;
using namespace TSFExtendedOps;
using namespace Teuchos;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;

static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

      int np = MPIComm::world().getNProc();

      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10*np, meshType);

      Mesh mesh = mesher.getMesh();

      CellFilter interior = new MaximalCellFilter();

      VectorType<double> vecType = new EpetraVectorType();

      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

      DiscreteFunction* discFunc = dynamic_cast<DiscreteFunction*>(u0.ptr().get());

      
      Vector<double> x0 = discFunc->vector();

      QuadratureFamily quad = new GaussianQuadrature(2);
      
      Expr eqn = Integral(interior, v*(u*u-2.0), quad);
      Expr bc;

      RefCountPtr<EquationSet> eqnSet 
        = rcp(new EquationSet(eqn, bc, v, u, u0, 
                              rcp(new BruteForceEvaluatorFactory())));

      Assembler assembler(mesh, eqnSet, vecType); 

      LinearOperator<double> A;
      Vector<double> b;
      Vector<double> solnVec;

      ParameterList solverParams;

      solverParams.set(LinearSolverBase<double>::verbosityParam(), 0);
      solverParams.set(IterativeSolver<double>::maxitersParam(), 100);
      solverParams.set(IterativeSolver<double>::tolParam(), 1.0e-12);

      LinearSolver<double> solver = new BICGSTABSolver<double>(solverParams);

      bool converged = false;
      for (int i=0; i<20; i++)
        {
          assembler.assemble(A, b);
          SolverState<double> state = solver.solve(A, b, solnVec);

          //          cerr << "solver state = " << endl << state << endl;

          x0 = x0 - solnVec;
          //          cerr << "updated soln = " << x0 << endl;
          cerr << "step norm = " << solnVec.norm2() << endl;

          discFunc->setVector(x0);

          if (solnVec.norm2() < 1.0e-14) 
            {
              cerr << "Newton's method converged!" << endl;
              converged = true;
              break;
            }
        }
      
      if (!converged) 
        {
          cerr << "FAILED TO CONVERGE!" << endl;
        }

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
