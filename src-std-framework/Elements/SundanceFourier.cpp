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

#include "SundanceFourier.hpp"
#include "SundanceADReal.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Teuchos;
using std::endl;

// ctor
Fourier::Fourier(int N) : N_(N)
{}

bool Fourier::supportsCellTypePair(
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
    default:
      return false;
  }
}

void Fourier::print(std::ostream& os) const
{
  os << "Fourier(" << N_ << ")";
}

int Fourier::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(cellType)
  {
    case PointCell:
      return 0;
    case LineCell:
      return 2*N_+1;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument , "illegal combination of cell type and maximal cell type" );
      return -1;
  }
}

Array<int> Fourier::makeRange(int low, int high)
{
  if (high < low) return Array<int>();

  Array<int> rtn(high-low+1);
  for (int i=0; i<rtn.length(); i++) rtn[i] = low+i;
  return rtn;
}

void Fourier::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  typedef Array<int> Aint;

  switch(cellType)
  {
    case LineCell:
    {
      dofs.resize(2);
      dofs[0] = Array<Array<int> >();
      dofs[1] = tuple(makeRange(0, 2*N_));
    }
    break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
        << cellType << " not implemented in Fourier basis");
  }
}


void Fourier::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    std::runtime_error,
    "cannot evaluate spatial derivative " << sds << " on Fourier basis");
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
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Fourier::refEval() unimplemented for cell type "
        << cellType);
  }
}

/* ---------- evaluation on different cell types -------------- */

void Fourier::evalOnLine(const Point& pt,
  const MultiIndex& deriv,
  Array<double>& result) const
{
  result.resize(2*N_+1);
  double x = pt[0];
  const double pi = 4.0*atan(1.0);

  if (deriv.order()==0)
  {
    result[0] = 1.0;
    for (int n=1; n<=N_; n++)
    {
      result[2*n-1]=::cos(2.0*pi*n*x);
      result[2*n]=::sin(2.0*pi*n*x);
    }
  }
  else if (deriv.order()==1)
  {
    result[0] = 0.0;
    for (int n=1; n<=N_; n++)
    {
      result[2*n-1]=-2.0*pi*n*::sin(2.0*pi*n*x);
      result[2*n]=2.0*n*pi*::cos(2.0*pi*n*x);
    }
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  Out::os() << "quad point=" << x << endl;
  Out::os() << "bas: " << result << endl;
}

