/* @HEADER@ */
/* @HEADER@ */

#include "SundanceAssembler.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
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


Assembler::Assembler(const Mesh& mesh, 
                     const RefCountPtr<EquationSet>& eqn)
  : mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    bcRows_(),
    rqc_(),
    isBCRqc_(),
    rqcExprs_(),
    rqcDerivSet_(),
    rqcEval_(),
    weakForms_(),
    evalMgr_(rcp(new EvalManager())),
    rqcEvaluatableExpr_()
{
  verbosity() = Assembler::classVerbosity();

  DOFMapBuilder mapBuilder(mesh, eqn);
  rowMap_ = mapBuilder.rowMap();
  colMap_ = mapBuilder.colMap();
  bcRows_ = mapBuilder.bcRows();

  for (int r=0; r<eqn_->regionQuadCombos().size(); r++)
    {
      rqc_.append(eqn_->regionQuadCombos()[r]);
      isBCRqc_.append(false);
      rqcExprs_.append(eqn_->expr(eqn_->regionQuadCombos()[r]));
      rqcDerivSet_.append(eqn_->nonzeroFunctionalDerivs(eqn_->regionQuadCombos()[r]));
      addToWeakFormBatch(eqn_->nonzeroFunctionalDerivs(eqn_->regionQuadCombos()[r]));
    }
  for (int r=0; r<eqn_->bcRegionQuadCombos().size(); r++)
    {
      rqc_.append(eqn_->bcRegionQuadCombos()[r]);
      isBCRqc_.append(true);
      rqcExprs_.append(eqn_->bcExpr(eqn_->bcRegionQuadCombos()[r]));
      rqcDerivSet_.append(eqn_->nonzeroBCFunctionalDerivs(eqn_->bcRegionQuadCombos()[r]));
      addToWeakFormBatch(eqn_->nonzeroBCFunctionalDerivs(eqn_->bcRegionQuadCombos()[r]));
    }

  for (int r=0; r<rqc_.size(); r++)
    {
      SUNDANCE_OUT(verbosity() > VerbLow, "rqc[" << r << "] = " << rqc_[r]);
      CellFilter filter = rqc_[r].domain();
      QuadratureFamily quad = rqc_[r].quad();
      int cellDim = filter.dimension(mesh);
      rqcEval_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, quad)));
      const EvaluatableExpr* ee 
        = dynamic_cast<const EvaluatableExpr*>(rqcExprs_[r][0].ptr().get());
      TEST_FOR_EXCEPTION(ee==0, RuntimeError, "expression " << rqcExprs_[r]
                         << " could not be cast to evaluatable expr");
      rqcEvaluatableExpr_.append(ee);
    }
}

void Assembler::addToWeakFormBatch(const DerivSet& derivs) 
{
  Array<RefCountPtr<WeakFormBatch> > w;
  if (derivs.size() > 0) 
    {
      DerivSet::const_iterator i;

      int pos = 0;
      for (i=derivs.begin(); i != derivs.end(); i++, pos++)
        {
          const MultipleDeriv& d = *i;
          if (d.order()==0) continue;
          bool foundMatch = false;
          for (int j=0; j<w.size(); j++)
            {
              if (w[j]->tryToAdd(d, pos))
                {
                  foundMatch = true;
                  break;
                }
            }
          if (!foundMatch) 
            {
              w.append(rcp(new WeakFormBatch(d, pos)));
            }
        }
    }
  weakForms_.append(w);
}

                                 


void Assembler::print(ostream& os) const 
{
  for (int r=0; r<rqc_.size(); r++)
    {
      os << "Region/Quad combination: " << rqc_[r] << endl;
      {
        Tabs tab;
        os << tab << "isBC = " << Teuchos::toString(isBCRqc_[r]) << endl;
        os << tab << "nonzero derivs = " << rqcDerivSet_[r] << endl;
        os << tab << "Weak forms: " << endl;
        for (int w=0; w<weakForms_[r].size(); w++)
          {
            Tabs tab1 ;
            os << tab1;
            weakForms_[r][w]->print(os);
            os << endl;
          }
      }
    }
}

void Assembler::assemble() const 
{
  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());
  RefCountPtr<EvalVectorArray> results = rcp(new EvalVectorArray());

  for (int r=0; r<rqc_.size(); r++)
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "assembling rqc=" << rqc_[r]);
      evalMgr_->setMediator(rqcEval_[r]);
      evalMgr_->setRegion(rqc_[r]);
      CellFilter filter = rqc_[r].domain();
      CellSet cells = filter.getCells(mesh_);
      int cellDim = filter.dimension(mesh_);
      CellType cellType = mesh_.cellType(cellDim);

      CellIterator iter=cells.begin();
      rqcEval_[r]->setCellType(cellType);
      int numWorkSets = 0;
      while (iter != cells.end())
        {
          workSet->resize(0);
          for (int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       "doing work set=" << numWorkSets);
          numWorkSets++;
          rqcEval_[r]->setCellBatch(workSet);
          rqcEvaluatableExpr_[r]->flushResultCache();
          rqcEvaluatableExpr_[r]->evaluate(*evalMgr_, results);
          if (verbosity() > VerbHigh) dumpResults(rqcEval_[r], results,
                                                  rqcDerivSet_[r]);
        }
      
    }
}

void Assembler::dumpResults(const RefCountPtr<StdFwkEvalMediator>& eval,
                            const RefCountPtr<EvalVectorArray>& results,
                            const DerivSet& derivs) const
{
  eval->print(cerr);
  results->print(cerr, derivs);
}
