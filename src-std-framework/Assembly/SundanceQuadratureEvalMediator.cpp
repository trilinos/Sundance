/* @HEADER@ */
/* @HEADER@ */

#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
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
    physQuadPts_()
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
                                           LoadableVector* const vec) const
{
  SUNDANCE_OUT(verbosity() > VerbSilent, "evaluating coord expr " << expr->toXML().toString());
  
  computePhysQuadPts();
  int nQuad = physQuadPts_.length();
  int d = expr->dir();
  vec->resize(nQuad);
  for (int q=0; q<nQuad; q++) 
    {
      vec->setElement(q, physQuadPts_[q][d]);
    }
}

void QuadratureEvalMediator
::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                          const MultiIndex& mi,
                          LoadableVector* const vec) const
{
  SUNDANCE_ERROR("QuadratureEvalMediator::evalDiscreteFuncElement() unimplemented");
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
