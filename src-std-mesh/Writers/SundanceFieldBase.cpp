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

#include "SundanceFieldBase.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceMesh.hpp"

using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace Teuchos;

using SundanceStdFwk::Internal::CellSet;
using SundanceStdFwk::Internal::CellIterator;
using SundanceStdFwk::CellFilter;
using SundanceStdFwk::MaximalCellFilter;


void FieldBase::getDataBatch(
  int cellDim, 
  const Array<int>& cellID,
  const Array<int>& funcElem,
  Array<double>& batch) const
{
  /* This is a silly default implementation */
  
  int nFunc = funcElem.size();
  batch.resize(cellID.size() * funcElem.size());

  for (unsigned int c=0; c<cellID.size(); c++)
  {
    for (unsigned int f=0; f<funcElem.size(); f++)
    {
      batch[c*nFunc + f] = getData(cellDim, cellID[c], f);
    }
  }
}

const CellFilter& FieldBase::domain() const 
{
  static CellFilter dum = new MaximalCellFilter();
  return dum;
}


namespace SundanceStdMesh
{

using SundanceStdFwk::Internal::CellSet;

CellSet connectedNodeSet(const CellFilter& f, const Mesh& mesh)
{
  CellSet cells = f.getCells(mesh);
  int dim = cells.dimension();
  if (dim==0) return cells;


  Array<int> cellLID;

  for (CellIterator i=cells.begin(); i!=cells.end(); i++)
  {
    cellLID.append(*i);
  }

  Array<int> nodes;
  Array<int> fo;

  mesh.getFacetLIDs(dim, cellLID, 0, nodes, fo);

  Set<int> nodeSet;

  for (unsigned int i=0; i<nodes.size(); i++)
  {
    nodeSet.put(nodes[i]);
  }
  
  return CellSet(mesh, 0, PointCell, nodeSet);
}


RefCountPtr<Array<int> > cellSetToLIDArray(const CellSet& cs)
{
  RefCountPtr<Array<int> > cellLID = rcp(new Array<int>());

  for (CellIterator i=cs.begin(); i!=cs.end(); i++)
  {
    cellLID->append(*i);
  }
  
  return cellLID;
}

}
