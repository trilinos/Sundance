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

#include "SundanceStdFwkEvalMediator.hpp"
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




void StdFwkEvalMediator::setCellBatch(const RefCountPtr<Array<int> >& cellLID,
                                      RefCountPtr<CellJacobianBatch>& J) 
{
  cellLID_ = cellLID; 
  cacheIsValid() = false; 
  jCacheIsValid_=false;
  mesh_.getJacobians(cellDim(), *cellLID, *J);
  J_ = J;

  /* mark the function caches as invalid */
  Map<const DiscreteFunction*, bool>::iterator iter;
  for (iter = fCacheIsValid_.begin(); iter != fCacheIsValid_.end(); iter++)
    {
      iter->second = false;
    }
  for (iter = dfCacheIsValid_.begin(); iter != dfCacheIsValid_.end(); iter++)
    {
      iter->second = false;
    }
  for (iter = localValueCacheIsValid_.begin(); iter != localValueCacheIsValid_.end(); iter++)
    {
      iter->second = false;
    }
}
