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
                                           SundanceCore::Internal::LoadableVector* const vec) const
{
  SUNDANCE_OUT(verbosity() > VerbSilent, "evaluating coord expr " << expr->toXML().toString());
  
  computePhysQuadPts();
  int nQuad = physQuadPts_.length();
  int d = expr->dir();
  cerr << "num phys quad=" << nQuad << endl;
  vec->resize(nQuad);
  for (int q=0; q<nQuad; q++) 
    {
      vec->setElement(q, physQuadPts_[q][d]);
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
                          const MultiIndex& mi,
                          SundanceCore::Internal::LoadableVector* const vec) const
{
  const DiscreteFunction* f = dynamic_cast<const DiscreteFunction*>(expr->master());
  TEST_FOR_EXCEPTION(f==0, InternalError,
                     "QuadratureEvalMediator::evalDiscreteFuncElement() called "
                     "with expr that is not a discrete function");

  int myIndex = expr->myIndex();


  const BasisFamily& basis = f->basis()[myIndex];

  const Vector<double>& fValues = f->vector();
  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());
  f->getLocalValues(cellDim(), *cellLID(), *localValues);

  RefCountPtr<Array<Array<Array<double> > > > refBasisValues 
    = getRefBasisVals(basis, mi.order());

  RefCountPtr<CellJacobianBatch> J = rcp(new CellJacobianBatch());
  if (mi.order() != 0) mesh().getJacobians(cellDim(), *cellLID(), *J);
  
  int nQuad = quadWgts().size();
  int nNodes = basis.nNodes(cellType());
  int nFuncs = f->discreteSpace().nFunc();

  vec->resize(cellLID()->size() * nQuad);

  for (int c=0; c<cellLID()->size(); c++)
    {
      /* initialize to zero */
      for (int q=0; q<nQuad; q++) 
        {
          vec->start()[c*nQuad + q] = 0.0;
        }
      if (mi.order()==0)
        {
          for (int q=0; q<nQuad; q++)
            {
              double& sum = vec->start()[c*nQuad + q];
              for (int i=0; i<nNodes; i++)
                {
                  double coeff = (*localValues)[c*nNodes*nFuncs + nFuncs*i + myIndex];
                  double basisVals = (*refBasisValues)[0][q][i];
                  sum += coeff * basisVals;
                }
            }
        }
      else
        {
          Array<double> invJ;
          J->getInvJ(c, invJ);
          int pDir = mi.firstOrderDirection();
          
          /* initialize to zero */
          for (int q=0; q<nQuad; q++) 
            {
              vec->start()[c*nQuad + q] = 0.0;
            }
          for (int q=0; q<nQuad; q++)
            {
              double& sum = vec->start()[c*nQuad + q];
              for (int i=0; i<nNodes; i++)
                {
                  double g = (*localValues)[c*nNodes*nFuncs + nFuncs*i + myIndex];
                  for (int r=0; r<cellDim(); r++)
                    {
                      sum += g*invJ[pDir + r*cellDim()]*(*refBasisValues)[r][q][i];
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
