/* @HEADER@ */
/* @HEADER@ */

#include "SundanceIntegrator.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceCellJacobianBatch.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;


Integrator::Integrator(const Mesh& mesh, 
                       const Expr& expr,
                       const DerivSet& nonzeroDerivs,
                       const RegionQuadCombo& rqc,
                       const RefCountPtr<EvalManager>& evalMgr)
  : cellDim_(CellFilter(rqc.domain()).dimension(mesh)),
    cellType_(mesh.cellType(cellDim_)),
    mesh_(mesh),
    expr_(EvaluatableExpr::getEvalExpr(expr)),
    evalMgr_(evalMgr),
    nonzeroDerivs_(nonzeroDerivs),
    rqc_(rqc),
    mediator_(),
    sparsity_(expr_->sparsity(expr_->getDerivSetIndex(rqc)))
{
  mediator_ = rcp(new QuadratureEvalMediator(mesh, cellDim(), 
                                             QuadratureFamily(rqc.quad())));
  mediator_->setCellType(cellType_);

  // grouper_->findGroups(sparsity_, alpha_, beta_, constCoeff_, 
//                        elemIntegrals_);
}

void Integrator::integrate(const RefCountPtr<Array<int> >& workSet,
                           RefCountPtr<LocalMatrixContainer>& localMat) const
{
  RefCountPtr<EvalVectorArray> results;
  RefCountPtr<CellJacobianBatch> J;

  evalMgr_->setMediator(mediator_);
  evalMgr_->setRegion(rqc_);
  mediator_->setCellBatch(workSet);


  mesh().getJacobians(cellDim(), *workSet, *J);

  expr()->flushResultCache();

  expr()->evaluate(*evalMgr_, results);

  
}





