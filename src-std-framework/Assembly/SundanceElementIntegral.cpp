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

ElementIntegral::ElementIntegral(int dim, 
                                 const CellType& cellType)
  : dim_(dim),
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

ElementIntegral::ElementIntegral(int dim, 
                                 const CellType& cellType,
                                 const BasisFamily& testBasis,
                                 const Array<int>& alpha,
                                 int testDerivOrder)
  : dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(cellType)),
    unkDerivOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    nNodes_(nNodesTest_),
    order_(1),
    alpha_(alpha),
    beta_()
{;}



ElementIntegral::ElementIntegral(int dim,
                                 const CellType& cellType,
                                 const BasisFamily& testBasis,
                                 const Array<int>& alpha,
                                 int testDerivOrder,
                                 const BasisFamily& unkBasis,
                                 const Array<int>& beta,
                                 int unkDerivOrder)
  : dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(cellType)), 
    unkDerivOrder_(unkDerivOrder), 
    nRefDerivUnk_(ipow(dim, unkDerivOrder)),
    nNodesUnk_(unkBasis.nNodes(cellType)), 
    nNodes_(nNodesTest_*nNodesUnk_),
    order_(2),
    alpha_(alpha),
    beta_(beta)
{;}

int ElementIntegral::ipow(int base, int power) 
{
  int rtn = 1;
  for (int i=0; i<power; i++) rtn *= base;
  return rtn;
}

