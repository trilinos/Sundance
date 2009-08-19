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

#include "SundanceP1NC.hpp"
#include "SundanceADReal.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace Teuchos;



bool P1NC::supportsCellTypePair(
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


void P1NC::print(std::ostream& os) const 
{
  os << "P1NC()";
}

int P1NC::nReferenceDOFs(
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
                         << cellType << " not implemented in P1NC basis");
      return -1; // -Wall
    }
}

void P1NC::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > > &dofs ) const
{
  switch(cellType)
    {
    case PointCell:
      dofs.resize(1);
      dofs[0] = tuple(Array<int>());
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
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in P1NC basis");
    }
}



void P1NC::refEval(
    const CellType& maximalCellType,
    const CellType& cellType,
    const Array<Point>& pts,
    const MultiIndex& deriv,
    Array<Array<Array<double> > >& result) const
{
  result.resize(1);
  result[0].resize(pts.length());

  switch(cellType)
    {
    case LineCell:
      for (int i=0; i<pts.length(); i++)
        {
          evalOnLine(pts[i], deriv, result[0][i]);
        }
      return;
    case TriangleCell:
      for (int i=0; i<pts.length(); i++)
        {
          evalOnTriangle(pts[i], deriv, result[0][i]);
        }
      return;
    default:
#ifndef TRILINOS_7
      SUNDANCE_ERROR("P1NC::refEval() unimplemented for cell type "
                     << cellType);
#else
      SUNDANCE_ERROR7("P1NC::refEval() unimplemented for cell type "
                     << cellType);
#endif
    }
}

/* ---------- evaluation on different cell types -------------- */

void P1NC::evalOnLine(const Point& pt, 
                      const MultiIndex& deriv,
                      Array<double>& result) const
{
	ADReal x = ADReal(pt[0], 0, 1);
	ADReal one(1.0, 1);
	
	result.resize(1);
	Array<ADReal> tmp(result.length());

  tmp[0] = one;

	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else result[i] = tmp[i].gradient()[0];
		}
}

void P1NC::evalOnTriangle(const Point& pt, 
                          const MultiIndex& deriv,
                          Array<double>& result) const



{
	ADReal x = ADReal(pt[0], 0, 2);
	ADReal y = ADReal(pt[1], 1, 2);
	ADReal one(1.0, 2);

  Array<ADReal> tmp;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
               << y.value());

  result.resize(3);
  tmp.resize(3);

  tmp[0] = 1.0 - 2.0*y;
  tmp[1] = 2.0*(x + y) - 1.0;
  tmp[2] = 1.0 - 2.0*x;

	for (int i=0; i<tmp.length(); i++)
		{
      SUNDANCE_OUT(this->verb() > 3,
                   "tmp[" << i << "]=" << tmp[i].value() 
                   << " grad=" << tmp[i].gradient());
			if (deriv.order()==0) 
        result[i] = tmp[i].value();
			else 
        result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}



bool P1NC::lessThan(const BasisDOFTopologyBase* other) const 
{
  if (typeLessThan(this, other)) return true;
  if (typeLessThan(other, this)) return false;

  const P1NC* lag = dynamic_cast<const P1NC*>(other);
  return order() < lag->order();
}
