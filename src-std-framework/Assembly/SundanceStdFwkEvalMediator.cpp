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




StdFwkEvalMediator::StdFwkEvalMediator(const Mesh& mesh, int cellDim)
  : mesh_(mesh),
    cellDim_(cellDim),
    cellType_(NullCell),
    cellLID_(),
    useMaximalCells_(false),
    JVol_(rcp(new CellJacobianBatch())),
    JTrans_(rcp(new CellJacobianBatch())),
    facetIndices_(rcp(new Array<int>())),
    maxCellLIDs_(rcp(new Array<int>())),
    cacheIsValid_(false),
    jCacheIsValid_(false),
    fCache_(),
    dfCache_(),
    localValueCache_(),
    fCacheIsValid_(),
    dfCacheIsValid_(),
    localValueCacheIsValid_()
{;}

void StdFwkEvalMediator::setCellType(const CellType& cellType,
                               const CellType& maxCellType) 
{
  cellType_=cellType; 
  maxCellType_ = maxCellType_;
  cacheIsValid() = false; 
  jCacheIsValid_=false;
}

void StdFwkEvalMediator::setCellBatch(bool useMaximalCells,
                                      const RefCountPtr<Array<int> >& cellLID) 
{
  cellLID_ = cellLID; 
  cacheIsValid() = false; 
  jCacheIsValid_=false;
  mesh_.getJacobians(cellDim(), *cellLID, *JVol_);
  useMaximalCells_ = useMaximalCells;

  if (useMaximalCells_)
    {
      const Array<int>& cells = *cellLID;
      facetIndices_->resize(cells.size());
      maxCellLIDs_->resize(cells.size());
      for (unsigned int c=0; c<cells.size(); c++)
        {
          (*maxCellLIDs_)[c] 
            = mesh_.cofacetLID(cellDim(), cells[c], 0, (*facetIndices_)[c]);
          cout << "max cellLID=" << (*maxCellLIDs_)[c] 
               << ", facetIndex=" << (*facetIndices_)[c] << endl;
        }
      mesh_.getJacobians(mesh_.spatialDim(), *maxCellLIDs_, *JTrans_);
    }

  /* mark the function caches as invalid */
  Map<const DiscreteFunctionData*, bool>::iterator iter;
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



const CellJacobianBatch& StdFwkEvalMediator::JTrans() const
{
  /* If we're integrating a derivative on a boundary, JVol and JTrans will be
   * different. Otherwise, they'll be the same, and we use JVol for both
   * volume computations and vector transformations */
  if (useMaximalCells_) return *JTrans_;
  return *JVol_;
}

const Array<int>& StdFwkEvalMediator::dofCellLIDs() const
{
  if (useMaximalCells_) return *maxCellLIDs_;
  return *cellLID_;
}

const CellType& StdFwkEvalMediator::dofCellType() const 
{
  if (useMaximalCells_) return maxCellType_;
  return cellType_;
}
