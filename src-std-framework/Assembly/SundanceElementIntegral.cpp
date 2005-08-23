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

#include "SundanceElementIntegral.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;

static Time& transCreationTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("building integral transformation matrices"); 
  return *rtn;
}

ElementIntegral::ElementIntegral(int spatialDim,
                                 int dim, 
                                 const CellType& cellType)
  : spatialDim_(spatialDim),
    dim_(dim),
    testDerivOrder_(-1), 
    nRefDerivTest_(-1),
    nNodesTest_(-1),
    unkDerivOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    nNodes_(-1),
    order_(0),
    alpha_(),
    beta_()
{;}

ElementIntegral::ElementIntegral(int spatialDim,
                                 int dim, 
                                 const CellType& cellType,
                                 const BasisFamily& testBasis,
                                 int alpha,
                                 int testDerivOrder)
  : spatialDim_(spatialDim),
    dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(spatialDim, cellType)),
    unkDerivOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    nNodes_(nNodesTest_),
    order_(1),
    alpha_(alpha),
    beta_(-1)
{;}



ElementIntegral::ElementIntegral(int spatialDim,
                                 int dim,
                                 const CellType& cellType,
                                 const BasisFamily& testBasis,
                                 int alpha,
                                 int testDerivOrder,
                                 const BasisFamily& unkBasis,
                                 int beta,
                                 int unkDerivOrder)
  : spatialDim_(spatialDim),
    dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(spatialDim, cellType)), 
    unkDerivOrder_(unkDerivOrder), 
    nRefDerivUnk_(ipow(dim, unkDerivOrder)),
    nNodesUnk_(unkBasis.nNodes(spatialDim, cellType)), 
    nNodes_(nNodesTest_*nNodesUnk_),
    order_(2),
    alpha_(alpha),
    beta_(beta)
{;}


Array<double>& ElementIntegral::G(int alpha)
{
  static Array<Array<double> > rtn(3);

  return rtn[alpha];
}

Array<double>& ElementIntegral::G(int alpha, int beta)
{
  static Array<Array<double> > rtn(9);

  return rtn[3*alpha + beta];
}

int& ElementIntegral::transformationMatrixIsValid(int alpha, int beta)
{
  static Array<int> rtn(9, false);
  return rtn[3*alpha + beta];
}

int& ElementIntegral::transformationMatrixIsValid(int alpha)
{
  static Array<int> rtn(3, false);
  return rtn[alpha];
}

void ElementIntegral::invalidateTransformationMatrices()
{
  for (int i=0; i<3; i++)
    {
      transformationMatrixIsValid(i) = false;
      for (int j=0; j<3; j++)
        {
          transformationMatrixIsValid(i, j) = false;
        }
    }
}




int ElementIntegral::ipow(int base, int power) 
{
  int rtn = 1;
  for (int i=0; i<power; i++) rtn *= base;
  return rtn;
}

void ElementIntegral::createTwoFormTransformationMatrix(const CellJacobianBatch& J) const
{
  TimeMonitor timer(transCreationTimer());

  TEST_FOR_EXCEPTION(J.cellDim() != dim(), InternalError,
                     "Inconsistency between Jacobian dimension " << J.cellDim()
                     << " and cell dimension " << dim() 
                     << " in ElementIntegral::createTwoFormTransformationMatrix()");

  int flops = 0;

  if (testDerivOrder() == 1 && unkDerivOrder() == 1)
    {
      if (transformationMatrixIsValid(alpha(), beta())) return;
      transformationMatrixIsValid(alpha(), beta()) = true;

      G(alpha(), beta()).resize(J.numCells() * J.cellDim() * J.cellDim());

      SUNDANCE_OUT(this->verbosity() > VerbMedium, 
                   Tabs() << "both derivs are first order");

      double* GPtr = &(G(alpha(),beta())[0]);
      int k = 0;

      for (int c=0; c<J.numCells(); c++)
        {
          static Array<double> invJ;
          J.getInvJ(c, invJ);
          double detJ = fabs(J.detJ()[c]);
          for (int gamma=0; gamma<dim(); gamma++)
            {
              for (int delta=0; delta<dim(); delta++, k++)
                {
                  GPtr[k] =  detJ*invJ[alpha() + gamma*dim()]
                    * invJ[beta() + dim()*delta];
                }
            }
        }
      flops = 2 * J.numCells() * dim() * dim() + J.numCells();
    }

  else if (testDerivOrder() == 1 && unkDerivOrder() == 0)
    {
      if (transformationMatrixIsValid(alpha())) return;
      transformationMatrixIsValid(alpha()) = true;

      G(alpha()).resize(J.numCells() * J.cellDim());

      int k = 0;
      double* GPtr = &(G(alpha())[0]);

      for (int c=0; c<J.numCells(); c++)
        {
          static Array<double> invJ;
          J.getInvJ(c, invJ);
          double detJ = fabs(J.detJ()[c]);
          for (int gamma=0; gamma<dim(); gamma++,k++)
            {
              GPtr[k] = detJ*invJ[alpha() + dim() * gamma];
            }
        }
      flops = J.numCells() * dim() + J.numCells();
    }

  else 
    {
      if (transformationMatrixIsValid(beta())) return;
      transformationMatrixIsValid(beta()) = true;

      G(beta()).resize(J.numCells() * J.cellDim());

      int k = 0;
      double* GPtr = &(G(beta())[0]);

      for (int c=0; c<J.numCells(); c++)
        {
          static Array<double> invJ;
          J.getInvJ(c, invJ);
          double detJ = fabs(J.detJ()[c]);
          for (int gamma=0; gamma<dim(); gamma++,k++)
            {
              GPtr[k] = detJ*invJ[beta() + dim() * gamma];
            }
        }
      flops = J.numCells() * dim() + J.numCells();
    }

  addFlops(flops);
}


void ElementIntegral
::createOneFormTransformationMatrix(const CellJacobianBatch& J) const 
{
  TimeMonitor timer(transCreationTimer());

  if (transformationMatrixIsValid(alpha())) return;
  transformationMatrixIsValid(alpha()) = true;

  TEST_FOR_EXCEPTION(J.cellDim() != dim(), InternalError,
                     "Inconsistency between Jacobian dimension " << J.cellDim()
                     << " and cell dimension " << dim() 
                     << " in ElementIntegral::createOneFormTransformationMatrix()");

  int flops = J.numCells() * dim() + J.numCells();

  G(alpha()).resize(J.numCells() * J.cellDim());

  int k = 0;
  double* GPtr = &(G(alpha())[0]);

  for (int c=0; c<J.numCells(); c++)
    {
      Array<double> invJ;
      J.getInvJ(c, invJ);
      double detJ = fabs(J.detJ()[c]);
      for (int gamma=0; gamma<dim(); gamma++, k++)
        {
          GPtr[k] = detJ*invJ[alpha() + dim() * gamma]; 
        }
    }
  
  addFlops(flops);
}

