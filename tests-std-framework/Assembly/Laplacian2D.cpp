#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceTabs.hpp"
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
#include "SundanceBasicIntegrator.hpp"
#include "SundanceRefIntegral.hpp"
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


int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());

      int pMax = 2;
      int dim=2;

      verbosity<RefIntegral>() = VerbHigh;

      CellType cellType = TriangleCell;

      Point a = Point(1.0, 1.0);
      Point b = Point(1.2, 1.6);
      Point c = Point(0.8, 1.3);
      CellJacobianBatch JBatch;
      JBatch.resize(1, dim, dim);
      double* J = JBatch.jVals(0);
      J[0] = b[0] - a[0];
      J[1] = c[0] - a[0];
      J[2] = b[1] - a[1];
      J[3] = c[1] - a[1];
      
      Array<double> coeff = tuple(1.0);
      RefCountPtr<Array<double> > A = rcp(new Array<double>());
          
      for (int p=1; p<=pMax; p++)
        {
          Tabs tab;
          BasisFamily P = new Lagrange(p);
          RefIntegral ref(dim, cellType, P, 1, P, 1);
          Array<int> alpha = tuple(0, 1);
          Array<int> beta = tuple(0, 1);
          Array<double> coeff = tuple(1.0, 1.0);
          
          ref.transformTwoForm(JBatch, 
                               alpha, beta, coeff, A);
          cerr << tab << "transformed element" << endl;
          cerr << tab << "{";
          for (int r=0; r<ref.nNodesTest(); r++)
            {
              if (r!=0) cerr << ", ";
              cerr << "{";
              for (int c=0; c<ref.nNodesUnk(); c++)
                {
                  if (c!=0) cerr << ", ";
                  cerr << (*A)[r + ref.nNodesTest()*c];
                }
              cerr << "}";
            }
          cerr << "}" << endl;
        }

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
