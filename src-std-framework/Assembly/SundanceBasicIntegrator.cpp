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
{;}



RefCountPtr<StdFwkEvalMediator> 
BasicIntegrator::createEvalMediator(const Mesh& mesh, const RegionQuadCombo& rqc) const
{
  return rcp(new QuadratureEvalMediator(mesh, cellDim(), 
                                        QuadratureFamily(rqc.quad())));
}

void BasicIntegrator
::innerIntegrate(const RefCountPtr<Array<int> >& workSet,
                 RefCountPtr<LocalMatrixBatch>& localMat) const
{
  
}


void BasicIntegrator::init()
{
}
