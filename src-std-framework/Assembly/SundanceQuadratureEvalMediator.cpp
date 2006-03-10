/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"
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

void QuadratureEvalMediator::evalCellDiameterExpr(const CellDiameterExpr* expr,
                                                  RefCountPtr<EvalVector>& vec) const 
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs 
                       << "QuadratureEvalMediator evaluating cell diameter expr " 
                       << expr->toString());

  int nQuad = quadWgts().size();
  int nCells = cellLID()->size();

  SUNDANCE_VERB_HIGH(tabs << "number of quad pts=" << nQuad);
  Array<double> diameters;
  mesh().getCellDiameters(cellDim(), *cellLID(), diameters);

  vec->resize(nQuad*nCells);
  double * const xx = vec->start();
  int k=0;
  for (int c=0; c<nCells; c++)
    {
      double h = diameters[c];
      for (int q=0; q<nQuad; q++, k++) 
        {
          xx[k] = h;
        }
    }
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

RefCountPtr<Array<double> > QuadratureEvalMediator
::getRefBasisVals(const BasisFamily& basis, int diffOrder) const
{
  Tabs tab;
  RefCountPtr<Array<double> > rtn ;

  typedef OrderedPair<BasisFamily, CellType> key;

  if (!refBasisVals_[diffOrder].containsKey(key(basis, cellType())))
    {
      SUNDANCE_OUT(this->verbosity() > VerbMedium,
                   tab << "computing basis values on quad pts");
      rtn = rcp(new Array<double>());

      if (diffOrder==0)
        {
          Array<Array<double> > tmp;
          basis.ptr()->refEval(mesh().spatialDim(), cellType(), 
                               *(refQuadPts_.get(cellType())), 
                               MultiIndex(), tmp);
          /* the tmp array contains values indexed as [quad][node]. 
           * We need to put this into fortran order with quad index running
           * fastest */
          int nQuad = tmp.size();
          int nNodes = tmp[0].size();
          int nTot = tmp.size() * tmp[0].size();
          rtn->resize(nTot);
          
          for (int q=0; q<nQuad; q++)
            {
              for (int n=0; n<nNodes; n++)
                {
                  (*rtn)[n*nQuad + q] = tmp[q][n];
                }
            }
        }
      else
        {
          Array<Array<Array<double> > > tmp(cellDim());          
          for (int r=0; r<cellDim(); r++)
            {
              MultiIndex mi;
              mi[r]=1;
              basis.ptr()->refEval(mesh().spatialDim(), cellType(), 
                                   *(refQuadPts_.get(cellType())), 
                                   mi, tmp[r]);
            }
          /* the tmp array contains values indexed as [quad][node]. 
           * We need to put this into fortran order with quad index running
           * fastest */
          int dim = cellDim();
          int nQuad = tmp[0].size();
          int nNodes = tmp[0][0].size();
          int nTot = dim * nQuad * nNodes;
          rtn->resize(nTot);
          
          for (int r=0; r<dim; r++)
            {
              for (int q=0; q<nQuad; q++)
                {
                  for (int n=0; n<nNodes; n++)
                    {
                      (*rtn)[(n*nQuad + q)*dim + r] = tmp[r][q][n];
                    }
                }
            }
        }
      refBasisVals_[diffOrder].put(key(basis, cellType()), rtn);
    }
  else
    {
      SUNDANCE_OUT(this->verbosity() > VerbMedium,
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
  const DiscreteFunctionData* f = DiscreteFunctionData::getData(expr);
  TEST_FOR_EXCEPTION(f==0, InternalError,
                     "QuadratureEvalMediator::evalDiscreteFuncElement() called "
                     "with expr that is not a discrete function");
  Tabs tab;
   if (Evaluator::classVerbosity() > VerbHigh)
     {
       cerr << tab << "evaluting DF " << expr->name() << endl;
     }

   const RefCountPtr<DOFMapBase>& dofMap = f->map(); 

   int nQuad = quadWgts().size();
   int myIndex = expr->myIndex();
   int chunk = dofMap->chunkForFuncID(myIndex);
   int funcIndex = dofMap->indexForFuncID(myIndex);
   int nFuncs = dofMap->nFuncs(chunk);

   for (unsigned int i=0; i<multiIndices.size(); i++)
    {
      Tabs tab1;
      const MultiIndex& mi = multiIndices[i];
      if (Evaluator::classVerbosity() > VerbHigh)
        {
          cerr << tab1 << "evaluting DF for multiindex " << mi << endl;
          cerr << tab1 << "num cells = " << cellLID()->size() << endl;
          cerr << tab1 << "num quad points = " << quadWgts().size() << endl;
          cerr << tab1 << "my index = " << expr->myIndex() << endl;
          cerr << tab1 << "num funcs = " << f->discreteSpace().nFunc() << endl;
        }
      vec[i]->resize(cellLID()->size() * quadWgts().size());
  
      if (mi.order() == 0)
        {
          if (!fCache().containsKey(f) || !fCacheIsValid()[f])
            {
              fillFunctionCache(f, mi);
            }
          else
            {
            }

          const RefCountPtr<Array<Array<double> > >& cacheVals 
            = fCache()[f];

          const double* cachePtr = &((*cacheVals)[chunk][0]);
          double* vecPtr = vec[i]->start();
          
          int cellSize = nQuad*nFuncs;
          int offset = funcIndex*nQuad;
          int k = 0;
          for (unsigned int c=0; c<cellLID()->size(); c++)
            {
              for (int q=0; q<nQuad; q++, k++)
                {
                  vecPtr[k] = cachePtr[c*cellSize + offset + q];
                }
            }
        }
      else
        {
          if (!dfCache().containsKey(f) || !dfCacheIsValid()[f])
            {
              fillFunctionCache(f, mi);
            }
          else
            {
            }

          const RefCountPtr<Array<Array<double> > >& cacheVals 
            = dfCache()[f];

          int dim = cellDim();
          int pDir = mi.firstOrderDirection();
          const double* cachePtr = &((*cacheVals)[chunk][0]);
          double* vecPtr = vec[i]->start();

          int cellSize = nQuad*nFuncs*dim;
          int offset = myIndex * nQuad * dim;
          int k = 0;

          for (unsigned int c=0; c<cellLID()->size(); c++)
            {
              for (int q=0; q<nQuad; q++, k++)
                {
                  vecPtr[k] = cachePtr[c*cellSize + offset + q*dim + pDir];
                }
            }
        }
    }
}

void QuadratureEvalMediator::fillFunctionCache(const DiscreteFunctionData* f,
                                               const MultiIndex& mi) const 
{
  const RefCountPtr<DOFMapBase>& dofMap = f->discreteSpace().map();
  int diffOrder = mi.order();

  int flops = 0;
  double jFlops = CellJacobianBatch::totalFlops();

  RefCountPtr<Array<Array<double> > > cacheVals;
  if (mi.order()==0)
    {
      if (fCache().containsKey(f))
        {
          cacheVals = fCache().get(f);
        }
      else
        {
          cacheVals = rcp(new Array<Array<double> >(dofMap->nChunks()));
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
          cacheVals = rcp(new Array<Array<double> >(dofMap->nChunks()));
          dfCache().put(f, cacheVals);
        }
      dfCacheIsValid().put(f, true);
    }

  RefCountPtr<Array<Array<double> > > localValues;
  if (!localValueCacheIsValid().containsKey(f) 
      || !localValueCacheIsValid().get(f))
    {
      localValues = rcp(new Array<Array<double> >());
      f->getLocalValues(cellDim(), *cellLID(), *localValues);
      localValueCache().put(f, localValues);
      localValueCacheIsValid().put(f, true);
    }
  else
    {
      localValues = localValueCache().get(f);
    }

  
  for (int chunk=0; chunk<dofMap->nChunks(); chunk++)
    {
      const BasisFamily& basis = dofMap->basis(chunk);
      int nFuncs = dofMap->nFuncs(chunk);

      RefCountPtr<Array<double> > refBasisValues 
        = getRefBasisVals(basis, diffOrder);

      Array<double>& cache = (*cacheVals)[chunk];

      int nQuad = quadWgts().size();
      int nCells = cellLID()->size();
      int nNodes = basis.nNodes(mesh().spatialDim(), cellType());
      int nDir;

      if (mi.order()==1)
        {
          nDir = cellDim();
          cache.resize(cellLID()->size() * nQuad * cellDim() * nFuncs);
        }
      else
        {
          nDir = 1;
          cache.resize(cellLID()->size() * nQuad * nFuncs);
        }

      /* 
       * Sum over nodal values, which we can do with a matrix-matrix multiply
       * between the ref basis values and the local function values.
       */
      int nRowsA = nQuad*nDir;
      int nColsA = nNodes;
      int nColsB = nFuncs*nCells; 
      int lda = nRowsA;
      int ldb = nNodes;
      int ldc = lda;
      double alpha = 1.0;
      double beta = 0.0;
      double* A = &((*refBasisValues)[0]);
      double* B = &((*localValues)[chunk][0]);
      double* C = &((*cacheVals)[chunk][0]);
      
      dgemm_("n", "n", &nRowsA, &nColsB, &nColsA, &alpha, A, &lda, 
             B, &ldb, &beta, C, &ldc);
      
      
      /* Transform derivatives to physical coordinates */
      const CellJacobianBatch& J = JTrans();
      if (mi.order()==1)
        {
          int nRhs = nQuad * nFuncs;
          for (unsigned int c=0; c<cellLID()->size(); c++)
            {
              double* rhsPtr = &(C[(nRhs * nDir)*c]);
              J.applyInvJ(c, 0, rhsPtr, nRhs, false);
            }
        }
    }

  jFlops = CellJacobianBatch::totalFlops() - jFlops;
  addFlops(flops + jFlops);
}

void QuadratureEvalMediator::computePhysQuadPts() const 
{
  if (cacheIsValid()) 
    {
      SUNDANCE_OUT(this->verbosity() > VerbLow, 
                   "reusing cached phys quad points");
    }
  else
    {
      double jFlops = CellJacobianBatch::totalFlops();
      SUNDANCE_OUT(this->verbosity() > VerbLow, 
                   "computing phys quad points");
      const Array<Point>& refPts = *(refQuadPts_.get(cellType()));
      mesh().pushForward(cellDim(), *cellLID(), 
                         refPts, physQuadPts_); 

      addFlops(CellJacobianBatch::totalFlops() - jFlops);
      cacheIsValid() = true;
    }
  SUNDANCE_OUT(this->verbosity() > VerbMedium, 
               "phys quad: " << physQuadPts_);
}


void QuadratureEvalMediator::print(ostream& os) const 
{
  if (cacheIsValid())
    {
      Tabs tab0;
      os << tab0 << "Physical quadrature points" << endl;
      int i=0;
      for (unsigned int c=0; c<cellLID()->size(); c++)
        {
          Tabs tab1;
          os << tab1 << "cell " << c << endl;
          for (unsigned int q=0; q<physQuadPts_.size()/cellLID()->size(); q++, i++)
            {
              Tabs tab2;
              os << tab2 << "q=" << q << " " << physQuadPts_[i] << endl;
            }
        }
    }
}


