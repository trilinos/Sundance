/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFunctionalEvaluator.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellSet.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;


FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh, 
                                         const CellFilter& filter,
                                         const Expr& expr,
                                         const QuadratureFamily& quad,
                                         const VerbositySetting& verb)
  : mesh_(mesh), 
    filter_(filter), 
    expr_(expr), 
    evalExpr_(EvaluatableExpr::getEvalExpr(expr)), 
    quad_(quad),
    rqc_(rcp(new RegionQuadCombo(filter.ptr(), quad.ptr()))),
    mediator_(),
    evalMgr_(rcp(new EvalManager()))
{
  verbosity() = verb;

  int cellDim = filter_.dimension(mesh_);
  mediator_ = rcp(new QuadratureEvalMediator(mesh, cellDim, quad));
}

double FunctionalEvaluator::evaluate() const
{
  /* specify the mediator for this RQC */
  evalMgr_->setMediator(mediators_);
  evalMgr_->setRegion(rqc_);

  RefCountPtr<EvalVectorArray> values;

  CellSet cells = filter_.getCells(mesh_);
  int cellDim = filter_.dimension(mesh_);
  CellType cellType = mesh_.cellType(cellDim);

  mediator_->setCellType(cellType);  

  double localSum = 0.0;

  /* do the cells in batches of the work set size */

  CellIterator iter=cells.begin();
  int workSetCounter = 0;

  while (iter != cells.end())
    {
      Tabs tab1;
      /* build up the work set */
      workSet->resize(0);
      for (int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
        {
          workSet->append(*iter);
        }
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   tab1 << "doing work set " << workSetCounter
                   << " consisting of " << workSet->size() << " cells");
      workSetCounter++;

      mediator_->setCellBatch(workSet);

      evalExpr_->flushResultCache();
      evalExpr_->evaluateAndSum(*evalMgr_, values);

      for (int c=0; c<workSet->size(); c++)
        {
          for (int q=0; q<mediator_->quadWeights().size(); q++)
            {
              localSum += detJ*quadWeights[q] * (*localValues)[c*nQuad + q];
            }
        }
    }

  
}
