/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBruteForceIntegrator.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceQuadratureEvalMediator.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;



BruteForceIntegrator::BruteForceIntegrator(const Mesh& mesh, 
                                 const Expr& expr,
                                 const DerivSet& nonzeroDerivs,
                                 const RegionQuadCombo& rqc,
                                 const RefCountPtr<EvalManager>& evalMgr)
  : IntegratorBase(mesh, expr, nonzeroDerivs, rqc, evalMgr)
{

  for (int i=0; i<sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& deriv = sparsity()->deriv(i);
      int testID;
      int unkID;
      MultiIndex miTest;
      MultiIndex miUnk;
      BasisFamily testBasis;
      BasisFamily unkBasis;
      bool isTwoForm = getWeakForm(deriv, testID, unkID, miTest, miUnk,
                                   testBasis, unkBasis);
      isTwoForm_.append(isTwoForm);
      testID_.append(testID);
      unkID_.append(unkID);      
      miTest_.append(miTest);
      miUnk_.append(miUnk);
      testBasis_.append(testBasis);
      unkBasis_.append(unkBasis);
    }
}



RefCountPtr<StdFwkEvalMediator> 
BruteForceIntegrator::createEvalMediator(const Mesh& mesh, const RegionQuadCombo& rqc) const
{
  return rcp(new QuadratureEvalMediator(mesh, cellDim(), 
                                        QuadratureFamily(rqc.quad())));
}

void BruteForceIntegrator
::innerIntegrate(const RefCountPtr<Array<int> >& workSet,
                 RefCountPtr<LocalMatrixBatch>& localMat) const
{
  CellJacobianBatch J;
  getJacobians(workSet, J);
  RefCountPtr<EvalVectorArray> results;
  evaluate(results);

  


  
}


void BruteForceIntegrator::init()
{
}


                                  
