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

      int pMax = 1;
      int maxDim=2;
      Array<int> alpha(1);
      Array<int> beta(1);

      verbosity<RefIntegral>() = VerbMedium;

      for (int dim=1; dim<=maxDim; dim++)
        {
          CellType cellType;
          if (dim==1) cellType=LineCell;
          if (dim==2) cellType=TriangleCell;

          /* do one-forms */
          for (int p=0; p<=pMax; p++)
            {
              BasisFamily P = new Lagrange(p);
              for (int d=0; d<=1; d++)
                {
                  RefIntegral ref(dim, cellType, P, alpha, d);
                }
            }
          
          /* do two-forms */
          for (int p=0; p<=pMax; p++)
            {
              BasisFamily P = new Lagrange(p);
              for (int q=0; q<=pMax; q++)
                {
                  BasisFamily Q = new Lagrange(q);
                  for (int dp=0; dp<=1; dp++)
                    {
                      for (int dq=0; dq<=1; dq++)
                        {
                          RefIntegral ref(dim, cellType, P, alpha, dp, 
                                          Q, beta, dq);
                        }
                    }
                }
            }
        }

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
