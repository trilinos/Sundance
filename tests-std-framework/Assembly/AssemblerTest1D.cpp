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
#include "TSFVectorType.hpp"
#include "TSFEpetraVectorType.hpp"

using namespace TSFExtended;
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

bool leftPointTest(const Point& x) {return fabs(x[0]) < 1.0e-10;}

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

        int np = MPIComm::world().getNProc();

      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10*np, meshType);

      //      mesher.ptr()->verbosity() = VerbExtreme;

      mesher.serializeLocal();

      Mesh mesh = mesher.getMesh();

      Array<int> funcs = tuple(0);

      //verbosity<CellFilter>() = VerbExtreme;
      //verbosity<CellSet>() = VerbExtreme;

      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellPredicate leftPointFunc = new PositionalCellPredicate(leftPointTest);
      CellFilter leftPoint = points.subset(leftPointFunc);
      
      Expr x = new CoordExpr(0);
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      //      Expr u0 = new DiscreteFunction(new Lagrange(1), "u0");
      Expr u0 = new ZeroExpr();
      Expr v = new TestFunction(new Lagrange(1), "v");
      Expr dx = new Derivative(0);
      

      verbosity<Assembler>() = VerbExtreme;
      verbosity<EquationSet>() = VerbExtreme;

      QuadratureFamily quad = new GaussianQuadrature(2);
      Expr eqn = Integral(interior, (dx*v)*(dx*u) + v, quad);
      Expr bc = EssentialBC(leftPoint, v*u, quad);

      RefCountPtr<EquationSet> eqnSet 
        = rcp(new EquationSet(eqn, bc, v, u, u0, 
                              rcp(new BruteForceEvaluatorFactory())));

      //      verbosity<EvalVector>() = VerbExtreme;
      // verbosity<QuadratureEvalMediator>() = VerbExtreme;
      //verbosity<EvaluatableExpr>() = VerbExtreme;
      //verbosity<Evaluator>() = VerbExtreme;
      //      EvalVector::shadowOps() = true;

      EquationSet::classVerbosity() = VerbHigh;
      Assembler::classVerbosity() = VerbHigh;
      IntegralGroup::classVerbosity() = VerbHigh;
      Expr::showAllParens() = true;

      VectorType<double> vecType = new EpetraVectorType();

      Assembler assembler(mesh, eqnSet, vecType); 

      Array<Set<int> > graph;
      assembler.getGraph(graph);

      cerr << "graph" << endl;
      for (int i=0; i<graph.size(); i++) 
        {
          cerr << "row=" << i << " " << graph[i] << endl;
        }

      LinearOperator<double> A;
      Vector<double> b;

      assembler.assemble(A, b);

      cerr << "Matrix = " << endl << A << endl;
      cerr << "RHS = " << endl << b << endl;

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
