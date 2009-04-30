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

#include "SundanceRaviartThomas.hpp"
#include "SundancePoint.hpp"
#include "SundanceCellType.hpp"
#include "SundanceMultiIndex.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceTypeUtils.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace TSFExtended;

RaviartThomas::RaviartThomas(int spatialDim)
  : HDivVectorBasis(spatialDim)
{}

bool RaviartThomas::supportsCellTypePair(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case TriangleCell:
      switch(cellType)
      {
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    case TetCell:
      switch(cellType)
      {
        case TetCell:
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    default:
      return false;
  }
}

std::string RaviartThomas::description() const 
{
  return "RaviartThomas()";
}

int RaviartThomas::nReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(cellType)
  {
    case PointCell:
      return 0;
    case LineCell:
      return 2;
    case TriangleCell:
      return 3;
    case TetCell:
      return 4;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
        << cellType << " not implemented in RaviartThomas basis");
      return -1; // -Wall
  }
}

void RaviartThomas::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  switch(cellType)
  {
    case PointCell:
      dofs.resize(1);
      dofs[0] = tuple(Array<int>());
      return;
    case LineCell:
      dofs.resize(2);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Array<int> >(tuple(0));
      return;
    case TriangleCell:
      dofs.resize(3);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Array<int> >(tuple(0), tuple(1), tuple(2));
      dofs[2] = tuple(Array<int>());
      return;
    case TetCell:
      dofs.resize(4);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Array<int> >(tuple(0), tuple(1), tuple(2), tuple(3));
      dofs[2] = tuple(Array<int>());
      dofs[3] = tuple(Array<int>());
      return;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
        << cellType << " not implemented in RaviartThomas basis");
  }
}



bool RaviartThomas::lessThan(const BasisDOFTopologyBase* other) const 
{
  if (typeLessThan(this, other)) return true;
  if (typeLessThan(other, this)) return false;

  return false;
}


void RaviartThomas::refEval(
  const CellType& maximalCellType,
  const CellType& cellType,
  const Array<Point>& pts,
  const MultiIndex& deriv,
  Array<Array<Array<double> > >& result) const
{
  TEST_FOR_EXCEPTION(true, RuntimeError, "evaluation of RaviartThomas elements not yet supported");
}


void RaviartThomas::print(std::ostream& os) const 
{
  os << "RaviartThomas(" << dim() << ")";
}
