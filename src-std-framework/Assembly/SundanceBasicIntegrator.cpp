/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBasicIntegrator.hpp"
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



BasicIntegrator::BasicIntegrator(const Mesh& mesh, 
                                 const Expr& expr,
                                 const DerivSet& nonzeroDerivs,
                                 const RegionQuadCombo& rqc,
                                 const RefCountPtr<EvalManager>& evalMgr)
  : IntegratorBase(mesh, expr, nonzeroDerivs, rqc, evalMgr)
{
//   Set<BasisFamily> basisSet;
//   Set<MultiIndex> miSet;
  

//   for (int i=0; i<sparsity()->numDerivs(); i++)
//     {
//       const MultipleDeriv& deriv = sparsity()->deriv(i);
//       WeakForm wf = getWeakForm(deriv, i, sparsity()->isConstant(i));
//       if (wf.isOneForm())
//         {
//           oneForms_.append(wf);
//         }
//       else
//         {
//           twoForms_.append(wf);
//         }
//     }
}



RefCountPtr<StdFwkEvalMediator> 
BasicIntegrator::createEvalMediator(const Mesh& mesh, const RegionQuadCombo& rqc) const
{
  return rcp(new QuadratureEvalMediator(mesh, cellDim(), 
                                        QuadratureFamily(rqc.quad())));
}

void BasicIntegrator
::innerIntegrate(const RefCountPtr<Array<int> >& workSet,
                 RefCountPtr<LocalMatrixContainer>& localMat) const
{
  // CellJacobianBatch J;
//   getJacobians(workSet, J);
//   RefCountPtr<EvalVectorArray> results;
//   evaluate(results);

//    for (int i=0; i<twoForms_.size(); i++)
//     {
//       const WeakForm& wf = twoForms_[i];
//       if (wf.hasConstantCoeff())
//         {
//           constantCoeffTwoFormIntegration(wf, workSet, localMat);
//         }
//       else
//         {
//           variableCoeffOneFormIntegration(wf, workSet, localMat);
//         }
//     }
  
}


void BasicIntegrator::init()
{
}


                                  
