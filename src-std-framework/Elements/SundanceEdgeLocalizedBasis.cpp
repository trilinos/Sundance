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

#include "SundanceEdgeLocalizedBasis.hpp"
#include "SundanceADReal.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace Teuchos;


EdgeLocalizedBasis::EdgeLocalizedBasis(int order)
  : order_(order)
{
TEST_FOR_EXCEPTION(order < 0, RuntimeError,
                     "invalid polynomial order=" << order
                     << " in EdgeLocalizedBasis ctor");
}

bool EdgeLocalizedBasis::supportsCellTypePair(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case TriangleCell:
    case TetCell:
      switch(cellType)
      {
        case LineCell:
          return true;
        default:
          return false;
      }
    default:
      return false;
  }
}

void EdgeLocalizedBasis::print(std::ostream& os) const 
{
  os << "EdgeLocalizedBasis(" << order_ << ")";
}

int EdgeLocalizedBasis::nReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(cellType)
  {
    case LineCell:
      return 1 + order_;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
        << cellType << " not implemented in EdgeLocalizedBasis basis");
      return -1; // -Wall
  }
}

void EdgeLocalizedBasis::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  typedef Array<int> Aint;
  switch(cellType)
  {
    case LineCell:
      dofs.resize(2);
      dofs[0] = Array<Array<int> >();
      dofs[1] = tuple<Aint>(makeRange(0, 1+order()));
      return;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
        << cellType << " not implemented in EdgeLocalizedBasis basis");
  }
}



Array<int> EdgeLocalizedBasis::makeRange(int low, int high)
{
  if (high < low) return Array<int>();

  Array<int> rtn(high-low+1);
  for (int i=0; i<rtn.length(); i++) rtn[i] = low+i;
  return rtn;
}

void EdgeLocalizedBasis::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    RuntimeError,
    "cannot evaluate spatial derivative " << sds << " on EdgeLocalizedBasis basis");
  const MultiIndex& deriv = sds.mi();
  typedef Array<double> Adouble;
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
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "EdgeLocalizedBasis::refEval() unimplemented for cell type "
                         << cellType);

    }
}

/* ---------- evaluation on different cell types -------------- */

void EdgeLocalizedBasis::evalOnLine(const Point& pt, 
													const MultiIndex& deriv,
													Array<double>& result) const
{
	ADReal x = ADReal(pt[0], 0, 1);
	ADReal one(1.0, 1);
	
	result.resize(order()+1);
	Array<ADReal> tmp(result.length());
  Array<double> x0(order()+1);

  if (order_ == 0)
    {
      tmp[0] = one;
    }
  else
    {
      x0[0] = 0.0;
      x0[1] = 1.0;
      for (int i=0; i<order_-1; i++)
        {
          x0[i+2] = (i+1.0)/order_;
        }

      for (int i=0; i<=order_; i++)
        {
          tmp[i] = one;
          for (int j=0; j<=order_; j++)
            {
              if (i==j) continue;
              tmp[i] *= (x - x0[j])/(x0[i]-x0[j]);
            }
        }
    }

	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else result[i] = tmp[i].gradient()[0];
		}
}


