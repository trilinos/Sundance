/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFunctionalEvaluator.hpp"
#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellSet.hpp"
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

static Time& functionalEvalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("functional evaluation"); 
  return *rtn;
}

static Time& functionalEvalSetupTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("functional evaluation setup"); 
  return *rtn;
}

namespace SundanceStdFwk
{
  double evaluateIntegral(const Mesh& mesh, const Expr& expr)
  {
    FunctionalEvaluator eval(mesh, expr);
    return eval.evaluate();
  }
}

FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh,
                                         const Expr& integral)
  : mesh_(mesh), 
    rqc_(),
    context_(),
    rqcExprs_(),
    mediators_(),
    evalExprs_(),
    evalMgr_(rcp(new EvalManager()))
{
  TimeMonitor timer(functionalEvalSetupTimer());
  verbosity() = classVerbosity();;


  /* begin with a sanity check to ensure that the input equation set 
   * exists and is integral form */
  SUNDANCE_OUT(verbosity() > VerbLow, 
               "checking existence of input eqn set...");

  TEST_FOR_EXCEPTION(integral.ptr().get()==0, RuntimeError,
                     "EquationSet ctor detected empty equation set input");

  const SumOfIntegrals* integralSum
    = dynamic_cast<const SumOfIntegrals*>(integral.ptr().get());

  TEST_FOR_EXCEPTION(integralSum==0, RuntimeError,
                     "EquationSet ctor detected an input equation set that "
                     "is not in integral form");
  SUNDANCE_OUT(verbosity() > VerbLow, 
               "...input eqn set is OK");

  Set<RegionQuadCombo> rqcSet;

  int contextID = EvalContext::nextID();

  for (int d=0; d<integralSum->numRegions(); d++)
    {
      /* make sure we have no test functions */
      OrderedHandle<CellFilterStub> reg = integralSum->region(d);

      TEST_FOR_EXCEPTION(integralSum->testsOnRegion(d).size() != 0,
                         RuntimeError,
                         "Expr with test function detected in functional "
                         "ctor in region " << reg);

      for (int t=0; t<integralSum->numTerms(d); t++)
        {
          int maxDiffOrder = 0;
          RegionQuadCombo rqc(reg.ptr(), integralSum->quad(d,t));
          EvalContext context(rqc, maxDiffOrder, contextID);
          context_.append(context);
          Expr term = integralSum->expr(d,t);
          rqcSet.put(rqc);
          rqcExprs_.put(rqc, term);
          DerivSet nonzeros = SymbPreprocessor::setupExpr(term, 
                                                          context);

        }
    }


  rqc_ = rqcSet.elements();

  for (int r=0; r<rqc_.size(); r++)
    {
      
      int cellDim = CellFilter(rqc_[r].domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      QuadratureFamily quad(rqc_[r].quad());
      const Expr& expr = rqcExprs_.get(rqc_[r]);
      const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
      evalExprs_.append(ee);
      mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
                                                       quad)));
    }

  

}

double FunctionalEvaluator::evaluate() const
{
  TimeMonitor timer(functionalEvalTimer());

  Array<RefCountPtr<EvalVector> > vectorValues;
  Array<double> constantValues;

  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  RefCountPtr<CellJacobianBatch> J = rcp(new CellJacobianBatch());

  double localSum = 0.0;
      
  for (int r=0; r<rqc_.size(); r++)
    {
      /* specify the mediator for this RQC */
      evalMgr_->setMediator(mediators_[r]);
      evalMgr_->setRegion(context_[r]);


      /* get the cells for the current domain */
      CellFilter filter = rqc_[r].domain();
      CellSet cells = filter.getCells(mesh_);
      int cellDim = filter.dimension(mesh_);
      CellType cellType = mesh_.cellType(cellDim);

      mediators_[r]->setCellType(cellType);  

      const Array<double>& quadWgts = mediators_[r]->quadWgts();
      int nQuad = quadWgts.size();



      const Evaluator* evaluator 
        = evalExprs_[r]->evaluator(context_[r]).get();

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

          mediators_[r]->setCellBatch(workSet, J);

          evaluator->resetNumCalls();
          evalExprs_[r]->evaluate(*evalMgr_, constantValues, vectorValues);
          
          //          mesh_.getJacobians(cellDim, *workSet, *J);
          const Array<double>& detJ = J->detJ();

          if (vectorValues.size()==0)
            {
              double f = constantValues[0];
              for (int c=0; c<workSet->size(); c++)
                {
                  localSum += fabs(detJ[c])*f;
                }
            }
          else
            {
              const double* const vals = vectorValues[0]->start();
              
              for (int c=0; c<workSet->size(); c++)
                {
                  double dJ = fabs(detJ[c]);
                  double cellSum = 0.0;
                  const double* const p = &(vals[c*nQuad]);
                  for (int q=0; q<nQuad; q++)
                    {
                      cellSum += quadWgts[q]*p[q];
                    }
                  localSum += dJ*cellSum;
                }
            }
        }
    }
  
  double globalSum = localSum;

  mesh_.comm().allReduce((void*) &localSum, (void*) &globalSum, 1, 
                         MPIComm::DOUBLE, MPIComm::SUM);
  return globalSum ; 
}
