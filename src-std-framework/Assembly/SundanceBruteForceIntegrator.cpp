/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBruteForceIntegrator.hpp"
#include "SundanceIPow.hpp"
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
      bool isConstant = sparsity()->isConstant(i);
      int testID;
      int unkID;
      MultiIndex miTest;
      MultiIndex miUnk;
      BasisFamily testBasis;
      BasisFamily unkBasis;
      bool isTwoForm = getWeakForm(deriv, testID, unkID, miTest, miUnk,
                                   testBasis, unkBasis);
      isTwoForm_.append(isTwoForm);
      isConstant_.append(isConstant);
      testID_.append(testID);
      unkID_.append(unkID);      
      miTest_.append(miTest);
      miUnk_.append(miUnk);
      testBasis_.append(testBasis);
      unkBasis_.append(unkBasis);
      
      if (isConstant)
        {
          if (isTwoForm)
            {
              integrateReferenceTwoForm(testBasis, miTest, unkBasis, miUnk, W2[0]);
            }
          else
            {
              integrateReferenceOneForm(testBasis, miTest, W1[0]);
            }
        }
      else
        {
          if (isTwoForm)
            {
              getReferenceTwoFormQuad(testBasis, miTest, unkBasis, miUnk, W2);
            }
          else
            {
              getReferenceOneFormQuad(testBasis, miTest, W1);
            }
        }
      
      oneFormCache_.append(W1);
      twoFormCache_.append(W2);
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
                 RefCountPtr<LocalMatrixContainer>& localMat) const
{
  CellJacobianBatch J;
  getJacobians(workSet, J);
  RefCountPtr<EvalVectorArray> results;
  evaluate(results);

  for (int i=0; i<sparsity()->numDerivs(); i++)
    {
      if (!isConstant_[i])
        {
          if (isTwoForm_[i])
            {
              sumTwoForms(twoFormQuad_, (*results)[i], twoFormCache_);
            }
          else
            {
              sumOneForms(oneFormQuad_, (*results)[i], oneFormCache_);
            }
        }
      if (isTwoForm_[i])
        {
          twoFormCache_[i]->transformToPhys(J, A);
        }
      else
        {
          oneFormCache_[i]->transformToPhys(J, b);
        }
    }
  
}


void BruteForceIntegrator::init()
{
}


                                  
