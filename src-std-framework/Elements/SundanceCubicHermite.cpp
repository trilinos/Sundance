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

#include "SundanceCubicHermite.hpp"
#include "SundanceADReal.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Teuchos;



bool CubicHermite::supportsCellTypePair(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case LineCell:
      switch(cellType)
      {
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
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
//     case TetCell:
//       switch(cellType)
//       {
//         case TetCell:
//         case TriangleCell:
//         case LineCell:
//         case PointCell:
//           return true;
//         default:
//           return false;
//       }
    default:
      return false;
  }
}

void CubicHermite::print(std::ostream& os) const 
{
  os << "CubicHermite";
}

int CubicHermite::nReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case LineCell:
      switch(cellType)
      {
        case PointCell:
          return 2;
        case LineCell:
          return 0;
        default:
          TEST_FOR_EXCEPTION( true , std::invalid_argument , "illegal combination of cell type and maximal cell type" );
          return -1;
      }
      break;
    case TriangleCell:
      switch(cellType)
      {
        case PointCell:
          return 3;
        case LineCell:
          return 0;
        case TriangleCell:
          return 1;
        default:
          TEST_FOR_EXCEPTION( true , std::invalid_argument , "illegal combination of cell type and maximal cell type" );
          return -1;
      }
      break;
    default:
      TEST_FOR_EXCEPTION( true , std::invalid_argument , "illegal combination of cell type and maximal cell type" );
      return -1;
  }
  
}

void CubicHermite::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  typedef Array<int> Aint;
  switch(cellType)
  {
    case PointCell:
    {
      std::cout << "woohoo" << std::endl;
      std::cout << maximalCellType << std::endl;
      dofs.resize(1);
      if (maximalCellType==LineCell)
        dofs[0] = tuple<Aint>(tuple(0,1));
      else if (maximalCellType==TriangleCell)
        dofs[0] = tuple<Aint>(tuple(0,1,2));
      else TEST_FOR_EXCEPT(1);
      return;
    }
    break;
    case LineCell:
    {
      dofs.resize(2);
      dofs[0].resize(2);
      if (maximalCellType==LineCell)
      {
        dofs[0][0].resize(2);
        dofs[0][0][0] = 0;
        dofs[0][0][1] = 1;
        dofs[0][1].resize(2);
        dofs[0][1][0] = 2;
        dofs[0][1][1] = 3;
      }
      else if (maximalCellType==TriangleCell)
      {
        dofs[0][0].resize(3);
        dofs[0][0][0] = 0;
        dofs[0][0][1] = 1;
        dofs[0][0][2] = 2;
        dofs[0][1].resize(3);
        dofs[0][1][0] = 3;
        dofs[0][1][1] = 4;
        dofs[0][1][1] = 5;
      }
      else
      {
        TEST_FOR_EXCEPT(1);
      }
      dofs[1].resize(1);
      dofs[1][0].resize(0);
      return;
    }
    break;
    case TriangleCell:
    {
      dofs.resize(3);
      dofs[0].resize(3);
      dofs[0][0] = tuple(0,1,2);
      dofs[0][1] = tuple(3,4,5);
      dofs[0][2] = tuple(6,7,8);
      dofs[1].resize(3);
      dofs[1][0].resize(0);
      dofs[1][1].resize(0);
      dofs[1][2].resize(0);
      dofs[2].resize(1);
      dofs[2][0].resize(1);
      dofs[2][0][0] = 9;
    }
    break;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
        << cellType << " not implemented in CubicHermite basis");
  }
}


void CubicHermite::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    RuntimeError,
    "cannot evaluate spatial derivative " << sds << " on CubicHermite basis");
  const MultiIndex& deriv = sds.mi();
  typedef Array<double> Adouble;
  result.resize(1);
  result[0].resize(pts.length());

  switch(cellType)
  {
    case PointCell:
      result[0] = tuple<Adouble>(tuple(1.0));
      return;
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
      TEST_FOR_EXCEPTION(true, RuntimeError,
        "CubicHermite::refEval() unimplemented for cell type "
        << cellType);

  }
}

/* ---------- evaluation on different cell types -------------- */

void CubicHermite::evalOnLine(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
  result.resize(4);
  ADReal x = ADReal(pt[0],0,1);
  Array<ADReal> tmp(4);

  tmp[0] = 1 + x * x * ( -3 + 2 * x );
  tmp[1] = x * ( 1 + (x - 2 ) * x );
  tmp[2] = ( 3 - 2*x ) * x * x;
  tmp[3] = (-1+x)*x*x;

  for (int i=0; i<tmp.length(); i++)
  {
    if (deriv.order()==0) 
    {
      result[i] = tmp[i].value();
    }
    else 
    {
      result[i] = tmp[i].gradient()[0];
    }
  }  
  return;
}

void CubicHermite::evalOnTriangle(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const

{
  result.resize(10);
  ADReal x = ADReal(pt[0], 0, 2);
  ADReal y = ADReal(pt[1], 1, 2);
  ADReal one(1.0, 2);
  
  Array<ADReal> tmp(10);

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
    << y.value());
 
  tmp[0] = (-1 + x + y) * (2 * x * x + (-1 + y) * (1 + 2 * y) + x*  (-1 + 11 *y));
  tmp[1] = x * (-1 + x + y) * (-1 + x + 2 * y);
  tmp[2] = y * (-1 + x + y) * (-1 + 2 * x + y);
  tmp[3] = x * (7 * (-1 + y) * y + x * (3 - 2 * x + 7 * y));
  tmp[4] = x * (x * (-1 + x - 2 * y) - 2 * (-1 + y) * y);
  tmp[5] = x * y * (-1 + 2 * x + y);
  tmp[6] = y * ((3 - 2 * y) * y + 7 * x * (-1 + x + y));
  tmp[7] = x * y * (-1 + x + 2 * y);
  tmp[8] = y * (-2 * x * x - 2 * x * (-1 + y) + (-1 + y) * y);
  tmp[9] = -27 * x * y * (-1 + x + y);

  for (int i=0; i<tmp.length(); i++)
  {
    if (deriv.order()==0) result[i] = tmp[i].value();
    else 
      result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
  }
}

void CubicHermite::evalOnTet(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
  ADReal x = ADReal(pt[0], 0, 3);
  ADReal y = ADReal(pt[1], 1, 3);
  ADReal z = ADReal(pt[2], 2, 3);
  ADReal one(1.0, 3);
  
  TEST_FOR_EXCEPT(true);
}




