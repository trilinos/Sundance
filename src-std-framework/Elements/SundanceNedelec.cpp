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

#include "SundanceNedelec.hpp"
#include "SundanceADReal.hpp"
#include "SundanceExceptions.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace TSFExtended;

Nedelec::Nedelec()
{}

bool Nedelec::supportsCellTypePair(
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
      // tets not yet implemented
      return false;
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

void Nedelec::print(ostream& os) const 
{
  os << "Nedelec()";
}

int Nedelec::nReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(cellType)
    {
    case PointCell:
      return 0;
    case LineCell:
      return 1;
    case TriangleCell:
      return 3;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in Nedelec basis");
      return -1; // -Wall
    }
}

void Nedelec::getReferenceDOFs(
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
      dofs[1] = tuple(tuple(0));
      return;
    case TriangleCell:
      dofs.resize(3);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple(tuple(0), tuple(1), tuple(2));
      dofs[2] = tuple(Array<int>());
      return;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in Nedelec basis");
    }
}



void Nedelec::refEval(
  const CellType& maximalCellType,
  const CellType& cellType,
  const Array<Point>& pts,
  const MultiIndex& deriv,
  Array<Array<Array<double> > >& result) const
{
  int dim = dimension(maximalCellType);
  result.resize(dim);
  for (int i=0; i<dim; i++) result[i].resize(pts.length());

  switch(cellType)
    {
    case TriangleCell:
      for (int d=0; d<dim; d++)
      {
        for (int i=0; i<pts.length(); i++)
        {
          evalOnTriangle(d, pts[i], deriv, result[d][i]);
        }
      }
      return;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "Nedelec::refEval() unimplemented for cell type "
                         << cellType);

    }
}

/* ---------- evaluation on different cell types -------------- */

void Nedelec::evalOnTriangle(int dir,
  const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
	ADReal x = ADReal(pt[0], 0, 2);
	ADReal y = ADReal(pt[1], 1, 2);
	ADReal one(1.0, 2);

  Array<ADReal> tmp;

  SUNDANCE_OUT(this->verbosity() > VerbHigh, "x=" << x.value() << " y="
               << y.value());

  result.resize(3);
  tmp.resize(3);

  switch(dir)
  {
    case 0:
      tmp[0] = 1.0-y;
      tmp[1] = -y;
      tmp[2] = -y;
      break;
    case 1:
      tmp[0] = x;
      tmp[1] = x;
      tmp[2] = x-1.0;
      break;
    default:
      TEST_FOR_EXCEPTION(true, InternalError, "impossible direction="
        << dir << " in Nedelec::evalOnTriangle()");
  }

	for (int i=0; i<tmp.length(); i++)
  {
    SUNDANCE_OUT(this->verbosity() > VerbHigh,
      "tmp[" << i << "]=" << tmp[i].value() 
      << " grad=" << tmp[i].gradient());
    if (deriv.order()==0) result[i] = tmp[i].value();
    else 
      result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
  }
}

