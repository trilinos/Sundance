/* @HEADER@ */
/* @HEADER@ */

#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;


QuadratureEvalMediator
::QuadratureEvalMediator(const Mesh& mesh, 
                         int cellDim,
                         const QuadratureFamily& quad)
  : StdFwkEvalMediator(mesh, cellDim),
    quad_(quad),
    refQuadPts_(),
    refQuadWeights_(),
    physQuadPts_(),
    refBasisVals_(2)
{}

void QuadratureEvalMediator::setCellType(const CellType& cellType) 
{
  StdFwkEvalMediator::setCellType(cellType);
  
  if (refQuadPts_.containsKey(cellType)) return;

  
  RefCountPtr<Array<Point> > pts = rcp(new Array<Point>());
  RefCountPtr<Array<double> > wgts = rcp(new Array<double>());

  quad_.getPoints(cellType, *pts, *wgts);
  refQuadPts_.put(cellType, pts);
  refQuadWeights_.put(cellType, wgts);
  
}

void QuadratureEvalMediator::evalCoordExpr(const CoordExpr* expr,
                                           RefCountPtr<EvalVector>& vec) const 
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs 
                       << "QuadratureEvalMediator evaluating coord expr " 
                       << expr->toString());
  
  computePhysQuadPts();
  int nQuad = physQuadPts_.length();
  int d = expr->dir();
  
  SUNDANCE_VERB_HIGH(tabs << "number of quad pts=" << nQuad);

  vec->resize(nQuad);
  double * const xx = vec->start();
  for (int q=0; q<nQuad; q++) 
    {
      xx[q] = physQuadPts_[q][d];
    }
}

RefCountPtr<Array<Array<Array<double> > > > QuadratureEvalMediator
::getRefBasisVals(const BasisFamily& basis, int diffOrder) const
{
  Tabs tab;
  RefCountPtr<Array<Array<Array<double> > > > rtn ;

  typedef OrderedPair<BasisFamily, CellType> key;

  if (!refBasisVals_[diffOrder].containsKey(key(basis, cellType())))
    {
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   tab << "computing basis values on quad pts");
      rtn = rcp(new Array<Array<Array<double> > >());
      if (diffOrder==0)
        {
          rtn->resize(1);
          basis.ptr()->refEval(cellType(), *(refQuadPts_.get(cellType())), 
                               MultiIndex(), (*rtn)[0]);
        }
      else
        {
          rtn->resize(cellDim());
          for (int r=0; r<cellDim(); r++)
            {
              MultiIndex mi;
              mi[r]=1;
              basis.ptr()->refEval(cellType(), *(refQuadPts_.get(cellType())), 
                                   mi, (*rtn)[r]);
            }
        }
      refBasisVals_[diffOrder].put(key(basis, cellType()), rtn);
    }
  else
    {
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   tab << "reusing basis values on quad pts");
      rtn = refBasisVals_[diffOrder].get(key(basis, cellType()));
    }
  return rtn;
}


void QuadratureEvalMediator
::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                          const Array<MultiIndex>& multiIndices,
                          Array<RefCountPtr<EvalVector> >& vec) const
{
  const DiscreteFunction* f = dynamic_cast<const DiscreteFunction*>(expr->master());
  TEST_FOR_EXCEPTION(f==0, InternalError,
                     "QuadratureEvalMediator::evalDiscreteFuncElement() called "
                     "with expr that is not a discrete function");

  for (int i=0; i<multiIndices.size(); i++)
    {
      const MultiIndex& mi = multiIndices[i];
      vec[i]->resize(cellLID()->size() * quadWgts().size());
  
      if (mi.order() == 0)
        {
          if (!fCache().containsKey(f) || !fCacheIsValid()[f])
            {
              fillFunctionCache(f, mi);
            }
          const RefCountPtr<Array<double> >& cacheVals 
            = fCache()[f];
          int nFuncs = f->discreteSpace().nFunc();
          int nPts = cellLID()->size() * quadWgts().size();
          int myIndex = expr->myIndex();
          const double* cachePtr = &((*cacheVals)[0]);
          double* vecPtr = vec[i]->start();
          for (int i=0; i<nPts; i++) 
            {
              vecPtr[i] = cachePtr[i*nFuncs + myIndex];
            }
        }
      else
        {
          if (!dfCache().containsKey(f) || !dfCacheIsValid()[f])
            {
              fillFunctionCache(f, mi);
            }
          const RefCountPtr<Array<double> >& cacheVals 
            = dfCache()[f];
          int nFuncs = f->discreteSpace().nFunc();
          int dim = cellDim();
          int nPts = cellLID()->size() * quadWgts().size();
          int pDir = mi.firstOrderDirection();
          int myIndex = expr->myIndex();
          /* offset to the first entry of the pDir'th derivative */
          const double* cachePtr = &((*cacheVals)[pDir*nPts*nFuncs]);
          double* vecPtr = vec[i]->start();
          for (int i=0; i<nPts; i++) 
            {
              vecPtr[i] = cachePtr[i*nFuncs + myIndex];
            }
        }
    }
}

void QuadratureEvalMediator::fillFunctionCache(const DiscreteFunction* f,
                                               const MultiIndex& mi) const 
{
  int nFuncs = f->discreteSpace().nFunc();
  int diffOrder = mi.order();

  RefCountPtr<Array<double> > cacheVals;
  if (mi.order()==0)
    {
      if (fCache().containsKey(f))
        {
          cacheVals = fCache().get(f);
        }
      else
        {
          cacheVals = rcp(new Array<double>());
          fCache().put(f, cacheVals);
        }
      fCacheIsValid().put(f, true);
    }
  else
    {
      if (dfCache().containsKey(f))
        {
          cacheVals = dfCache().get(f);
        }
      else
        {
          cacheVals = rcp(new Array<double>());
          dfCache().put(f, cacheVals);
        }
      dfCacheIsValid().put(f, true);
    }

  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());
  f->getLocalValues(cellDim(), *cellLID(), *localValues);

  Array<RefCountPtr<Array<Array<Array<double> > > > > refBasisValues(nFuncs) ;

  int nTotalNodes = 0;
  Array<int> nNodes(nFuncs);
  for (int i=0; i<nFuncs; i++)
    {
      const BasisFamily& basis = f->basis()[i];
      nTotalNodes += basis.nNodes(cellType());
      nNodes[i] = basis.nNodes(cellType());
      refBasisValues[i] = getRefBasisVals(basis, mi.order());
    }

  RefCountPtr<CellJacobianBatch> J = rcp(new CellJacobianBatch());
  if (mi.order() != 0) mesh().getJacobians(cellDim(), *cellLID(), *J);
  
  int nQuad = quadWgts().size();

  int nCells = cellLID()->size();
  int nDir;

  if (mi.order()==1)
    {
      nDir = cellDim();
      cacheVals->resize(cellLID()->size() * nQuad * cellDim() * nFuncs);
    }
  else
    {
      nDir = 1;
      cacheVals->resize(cellLID()->size() * nQuad * nFuncs);
    }


  for (int p=0; p<nDir; p++)
    {
      for (int c=0; c<nCells; c++)
        {
          const double* ptr = &((*localValues)[c*nTotalNodes]);
          double* vPtr = &((*cacheVals)[p*nCells*nFuncs*nQuad + c*nFuncs*nQuad]);
          int valsPerCell = nFuncs*nQuad;
          /* initialize to zero */
          for (int q=0; q<valsPerCell; q++) 
            {
              vPtr[q] = 0.0;
            }
          if (mi.order()==0)
            {
              for (int fid=0; fid<nFuncs; fid++)
                {
                  for (int q=0; q<nQuad; q++)
                    {
                      double& sum = vPtr[q*nFuncs + fid];
                      for (int i=0; i<nNodes[fid]; i++)
                        {
                          double coeff = ptr[nFuncs*i + fid];
                          double basisVals = (*(refBasisValues[fid]))[0][q][i];
                          sum += coeff * basisVals;
                        }
                    }
                }
            }
          else
            {
              static Array<double> invJ;
              J->getInvJ(c, invJ);
              const double* invJPtr = &(invJ[0]);
              int dim = cellDim();

              for (int fid=0; fid<nFuncs; fid++)
                {
                  for (int q=0; q<nQuad; q++)
                    {
                      double& sum = vPtr[q*nFuncs + fid];
                      for (int i=0; i<nNodes[fid]; i++)
                        {
                          double g = ptr[nFuncs*i + fid];
                          for (int r=0; r<dim; r++)
                            {
                              sum += g*invJPtr[p + r*dim]*(*(refBasisValues[fid]))[r][q][i];
                            }
                        }
                    }
                }
            }
        }
    }
}

void QuadratureEvalMediator::computePhysQuadPts() const 
{
  if (cacheIsValid()) 
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "reusing cached phys quad points");
    }
  else
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "computing phys quad points");
      const Array<Point>& refPts = *(refQuadPts_.get(cellType()));
      mesh().pushForward(cellDim(), *cellLID(), 
                         refPts, physQuadPts_); 
      cacheIsValid() = true;
    }
  SUNDANCE_OUT(verbosity() > VerbMedium, 
               "phys quad: " << physQuadPts_);
}

void QuadratureEvalMediator::print(ostream& os) const 
{
  if (cacheIsValid())
    {
      Tabs tab0;
      os << tab0 << "Physical quadrature points" << endl;
      int i=0;
      for (int c=0; c<cellLID()->size(); c++)
        {
          Tabs tab1;
          os << tab1 << "cell " << c << endl;
          for (int q=0; q<physQuadPts_.size()/cellLID()->size(); q++, i++)
            {
              Tabs tab2;
              os << tab2 << "q=" << q << " " << physQuadPts_[i] << endl;
            }
        }
    }
}
