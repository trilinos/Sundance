/* @HEADER@ */
/* @HEADER@ */

#include "SundanceIntegratorBase.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceEvaluatableExpr.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;


IntegratorBase::IntegratorBase(const Mesh& mesh, 
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
    sparsity_(expr_->sparsity(expr_->getDerivSetIndex(rqc))),
    needsInit_(true)
{}

void IntegratorBase::integrate(const RefCountPtr<Array<int> >& workSet,
                               RefCountPtr<LocalMatrixContainer>& localMat) const
{
  if (needsInit_) 
    {
      mediator_ = createEvalMediator(mesh_, rqc_);
      mediator_->setCellType(cellType_);
      const_cast<IntegratorBase&>(*this).init();
      needsInit_ = false;
    }

  evalMgr_->setMediator(mediator_);
  evalMgr_->setRegion(rqc_);
  mediator_->setCellBatch(workSet);
  innerIntegrate(workSet, localMat);
}

void IntegratorBase::getJacobians(const Array<int>& workSet,
                                  CellJacobianBatch& J) const 
{
  mesh().getJacobians(cellDim(), workSet, J);
}

void IntegratorBase::evaluate(RefCountPtr<EvalVectorArray>& results) const
{
  expr()->flushResultCache();
  expr()->evaluate(*evalMgr_, results);
}



