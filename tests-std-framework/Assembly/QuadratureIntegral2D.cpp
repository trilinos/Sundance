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
#include "SundanceFIATLagrange.hpp"
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
#include "SundanceQuadratureIntegral.hpp"
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


double chop(double x) {if (::fabs(x) < 1.0e-14) return 0.0; return x;}

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());

      int pMax = 1;
      int dim=2;

      verbosity<RefIntegral>() = VerbMedium;

      CellType cellType = TriangleCell;

      Point a = Point(0.0, 0.0);
      Point b = Point(1.0, 0.0);
      Point c = Point(0.0, 1.0);
      CellJacobianBatch JBatch;
      JBatch.resize(1, 2, 2);
      double* J = JBatch.jVals(0);
      J[0] = b[0] - a[0];
      J[1] = c[0] - a[0];
      J[2] = b[1] - a[1];
      J[3] = c[1] - a[1];
      
      RefCountPtr<Array<double> > A = rcp(new Array<double>());
          
      QuadratureFamily quad = new GaussianQuadrature(4);
      Array<double> quadWeights;
      Array<Point> quadPts;
      quad.getPoints(cellType, quadPts, quadWeights);
      int nQuad = quadPts.size();

      Array<double> coeff(nQuad);
      for (int i=0; i<nQuad; i++) 
        {
          double s = quadPts[i][0];
          double t = quadPts[i][1];
          double x = a[0] + J[0]*s + J[1]*t;
          double y = a[1] + J[2]*s + J[3]*t;
          coeff[i] = x*y;
        }
      const double* const f = &(coeff[0]);

      cerr << endl << endl 
           << "---------------- One-forms --------------------" 
           << endl << endl;
      for (int p=1; p<=pMax; p++)
        {
          BasisFamily P = new FIATLagrange(p);
          for (int dp=0; dp<=1; dp++)
            {
              if (dp > p) continue;

              int numTestDir = 1;
              if (dp==1) numTestDir = dim;
              for (int t=0; t<numTestDir; t++)
                {
                  Array<int> alpha = tuple(t);
                  Tabs tab;
                  QuadratureIntegral ref(dim, cellType, P, alpha, dp, quad);
                  ref.transformOneForm(JBatch, f, A);
                  cerr << tab << "transformed element" << endl;
                  cerr << tab << "t=" << t << endl;
                  cerr << tab << "{";
                  for (int r=0; r<ref.nNodesTest(); r++)
                    {
                      if (r!=0) cerr << ", ";
                      cerr << (*A)[r];
                    }
                  cerr << "}";
                }
            }
        }

      cerr << endl << endl 
           << "---------------- Two-forms --------------------" 
           << endl << endl;
      for (int p=1; p<=pMax; p++)
        {
          BasisFamily P = new FIATLagrange(p);
          for (int q=1; q<=pMax; q++)
            {
              BasisFamily Q = new FIATLagrange(q);
              for (int dp=0; dp<=1; dp++)
                {
                  if (dp > p) continue;
                  for (int dq=0; dq<=1; dq++)
                    {
                      if (dq > q) continue;

                      int numTestDir = 1;
                      if (dp==1) numTestDir = dim;
                      for (int t=0; t<numTestDir; t++)
                        {
                          Array<int> alpha = tuple(t);
                          int numUnkDir = 1;
                          if (dq==1) numUnkDir = dim;
                          for (int u=0; u<numUnkDir; u++)
                            {
                              Tabs tab;
                              Array<int> beta = tuple(u);
                              QuadratureIntegral ref(dim, cellType, P, alpha, 
                                                     dp, Q, beta, dq, quad);
                              ref.transformTwoForm(JBatch, f, A);
                              cerr << tab << "transformed element" << endl;
                              cerr << tab << "t=" << t << ", u=" << u << endl;
                              cerr << tab << "{";
                              for (int r=0; r<ref.nNodesTest(); r++)
                                {
                                  if (r!=0) cerr << ", ";
                                  cerr << "{";
                                  for (int c=0; c<ref.nNodesUnk(); c++)
                                    {
                                      if (c!=0) cerr << ", ";
                                      cerr << chop((*A)[r + ref.nNodesTest()*c]);
                                    }
                                  cerr << "}";
                                }
                              cerr << "}" << endl;
                            }
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
