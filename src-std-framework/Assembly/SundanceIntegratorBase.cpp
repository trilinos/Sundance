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

void IntegratorBase::integrateReferenceTwoForm(const BasisFamily& testBasis,
                                               const MultiIndex& miTest,
                                               const BasisFamily& unkBasis,
                                               const MultiIndex& miUnk,
                                               ReferenceTwoForm& W) const 
{
  int nRefDerivTest = ipow(cellDim(), miTest.order());
  int nRefDerivUnk = ipow(cellDim(), miUnk.order());
      
  int nNodesTest = testBasis.nNodes(cellType());
  int nNodesUnk = unkBasis.nNodes(cellType());
  
  Array<Array<Array<double> > > testBasisVals(nRefDerivTest);
  Array<Array<Array<double> > > unkBasisVals(nRefDerivUnk);

  W.setSize(nRefDerivTest, nNodesTest, nRefDerivUnk, nNodesUnk);
  
  for (int r=0; r<nRefDerivTest; r++)
    {
      MultiIndex mi;
      if (miTest.order()==1) mi[r] = 1;
      testBasis.ptr()->refEval(cellType(), quadPts, mi, testBasisVals[r]);
    }
  for (int r=0; r<nRefDerivUnk; r++)
    {
      MultiIndex mi;
      if (miUnk.order()==1) mi[r] = 1;
      unkBasis.ptr()->refEval(cellType(), quadPts, mi, unkBasisVals[r]);
    }
  
  for (int i=0; i<W.size(); i++) W[i] = 0.0;
  for (int q=0; q<nQuad; q++)
    {
      for (int t=0; t<nRefDerivTest; t++)
        {
          for (int nt=0; nt<nNodesTest; nt++)
            {
              for (int u=0; u<nRefDerivTest; u++)
                {
                  for (int nu=0; nu<nNodesUnk; nu++)
                    {
                      W(t, nt, u, nu) 
                        += quadWgts[q] * testBasisVals[t][q][nt] 
                        * unkBasisVals[u][q][nu];
                    }
                }
            }
        }
    }
}

