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

#include "SundanceLegendre.hpp"
#include "SundanceADReal.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Teuchos;

// ctor
Legendre::Legendre(int order):order_(order)
{

	// set the nr DOFs

	if (order_ > 1) nrDOF_edge_ = order_ - 1;
	else nrDOF_edge_ = 0;

	if (order_ > 3) nrDOF_face_ = ((order_-2)*(order_-3))/2;
	else nrDOF_face_ = 0;

	nrDOF_brick_ = 0;
}

bool Legendre::supportsCellTypePair(
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
    case QuadCell:
      switch(cellType)
      {
        case QuadCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
     case BrickCell:
       switch(cellType)
       {
         case BrickCell:
         case QuadCell:
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

void Legendre::print(std::ostream& os) const
{
  os << "Legendre";
}

int Legendre::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
      switch(cellType)
      {
        case PointCell:
          return 1;
        case LineCell:
          return nrDOF_edge_;
        case QuadCell:
            return nrDOF_face_;
        case BrickCell:
            return nrDOF_brick_;
        default:
            TEST_FOR_EXCEPTION( true , std::invalid_argument , "illegal combination of cell type and maximal cell type" );
            return -1;
      }
}

Array<int> Legendre::makeRange(int low, int high)
{
  if (high < low) return Array<int>();

  Array<int> rtn(high-low+1);
  for (int i=0; i<rtn.length(); i++) rtn[i] = low+i;
  return rtn;
}

void Legendre::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{

  typedef Array<int> Aint;

  //switch(order_)

  switch(cellType)
  {
    case PointCell:
    {
        dofs.resize(1);
        dofs[0] = tuple<Aint>(tuple(0));
        return;
    }
    break;
    case LineCell:
    {
        dofs.resize(2);
        dofs[0] = tuple<Aint>(tuple(0), tuple(1));
        if (nrDOF_edge_>0)
        {
        	dofs[1] = tuple<Aint>(makeRange(2, 2+nrDOF_edge_-1));
        }
        else
        {
        	dofs[1] = tuple(Array<int>());
        }
      return;
    }
    break;
    case QuadCell:
    {
    	int offs = 0;
        dofs.resize(3);
        // dofs[0] are DOFs at Points
        dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2), tuple(3));
        offs = 4;
        if (nrDOF_edge_>0)
        {
        	dofs[1].resize(4);
        	dofs[1][0] = makeRange(offs, offs+nrDOF_edge_-1);  offs += nrDOF_edge_;
        	dofs[1][1] = makeRange(offs, offs+nrDOF_edge_-1);  offs += nrDOF_edge_;
        	dofs[1][2] = makeRange(offs, offs+nrDOF_edge_-1);  offs += nrDOF_edge_;
        	dofs[1][3] = makeRange(offs, offs+nrDOF_edge_-1);  offs += nrDOF_edge_;
        }
        else
        {
        	dofs[1] = tuple(Array<int>(), Array<int>(),
	                 Array<int>(), Array<int>());
        }

        if (nrDOF_edge_>0)
        {
        	dofs[2].resize(1);
        	dofs[2][0] = makeRange(offs, offs+nrDOF_face_-1);  offs += nrDOF_face_;
        }
        else
        {
        	dofs[2] = tuple(Array<int>());
        }
        //SUNDANCE_OUT( true , "Legendre::getReferenceDOFs offs:" << offs );
    }
    break;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
        << cellType << " not implemented in Legendre basis");
  }
}


void Legendre::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    RuntimeError,
    "cannot evaluate spatial derivative " << sds << " on Legendre basis");
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
    case QuadCell:
      for (int i=0; i<pts.length(); i++)
      {
        evalOnQuad(pts[i], deriv, result[0][i]);
      }
      return;
    case BrickCell:
      for (int i=0; i<pts.length(); i++)
      {
        evalOnBrick(pts[i], deriv, result[0][i]);
      }
      return;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
        "Legendre::refEval() unimplemented for cell type "
        << cellType);

  }
}

/* ---------- evaluation on different cell types -------------- */

void Legendre::evalOnLine(const Point& pt,
  const MultiIndex& deriv,
  Array<double>& result) const
{
  result.resize(2+nrDOF_edge_);
  ADReal x = ADReal(pt[0],0,1);

  Array<ADReal> refAll(7);

  refAll[0] = 1-x;
  refAll[1] = x;
  refAll[2] = 2.44948974278318 * ( (2*x-1)*(2*x-1) - 1 ) / 4;
  refAll[3] = 3.16227766016838 * ( (2*x-1)*(2*x-1) - 1 ) * (2*x-1) / 4;
  refAll[4] = 3.74165738677394 * ( 5*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) - 6*(2*x-1)*(2*x-1) + 1) / 16;
  refAll[5] = 4.24264068711929 * (2*x-1) * (7*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) - 10*(2*x-1)*(2*x-1) + 3) / 16;
  refAll[6] = 4.69041575982343 * (21*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) -
		                          35*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) + 15*(2*x-1)*(2*x-1) - 1) / 32;


  for (int i=0; i<result.length(); i++)
  {
    if (deriv.order()==0) 
    {
      result[i] = refAll[i].value();
    }
    else 
    {
      result[i] = refAll[i].gradient()[0];
    }
  }  
  //SUNDANCE_OUT( true , "Legendre::evalOnLine result.length():" << result.length() );
  return;
}

void Legendre::evalOnQuad(const Point& pt,
  const MultiIndex& deriv,
  Array<double>& result) const

{
  result.resize( 4 + 4*nrDOF_edge_ + nrDOF_face_);
  ADReal x = ADReal(pt[0], 0, 2);
  ADReal y = ADReal(pt[1], 1, 2);
  ADReal one(1.0, 2);
  
  Array<ADReal> refAllx(7);
  Array<ADReal> refAlly(7);

  refAllx[0] = 1-x;
  refAllx[1] = x;
  refAllx[2] = 2.44948974278318 * ( (2*x-1)*(2*x-1) - 1 ) / 4;
  refAllx[3] = 3.16227766016838 * ( (2*x-1)*(2*x-1) - 1 ) * (2*x-1) / 4;
  refAllx[4] = 3.74165738677394 * ( 5*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) - 6*(2*x-1)*(2*x-1) + 1) / 16;
  refAllx[5] = 4.24264068711929 * (2*x-1) * (7*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) - 10*(2*x-1)*(2*x-1) + 3) / 16;
  refAllx[6] = 4.69041575982343 * (21*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) -
		                          35*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) + 15*(2*x-1)*(2*x-1) - 1) / 32;

  refAlly[0] = 1-y;
  refAlly[1] = y;
  refAlly[2] = 2.44948974278318 * ( (2*y-1)*(2*y-1) - 1 ) / 4;
  refAlly[3] = 3.16227766016838 * ( (2*y-1)*(2*y-1) - 1 ) * (2*y-1) / 4;
  refAlly[4] = 3.74165738677394 * ( 5*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1) - 6*(2*y-1)*(2*y-1) + 1) / 16;
  refAlly[5] = 4.24264068711929 * (2*y-1) * (7*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1) - 10*(2*y-1)*(2*y-1) + 3) / 16;
  refAlly[6] = 4.69041575982343 * (21*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1) -
		                          35*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1) + 15*(2*y-1)*(2*y-1) - 1) / 32;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
    << y.value());

  int sideIndex[4][2] = { {0,0} , {1,0} , {0,1} , {1,1}};
  int edgeIndex[4]    = { 0 , 1 , 1 , 0};
  int edgeMultI[4]    = { 0 , 0 , 1 , 1};
  int offs = 0;
  Array<ADReal> tmp(4 + 4*nrDOF_edge_ + nrDOF_face_);

  // loop over vertexes
  for (int i=0; i < 4 ; i++, offs++){
	  tmp[offs] = refAllx[sideIndex[i][0]] * refAlly[sideIndex[i][1]];
  }

  // loop over edges
  for (int i=0; i < 4 ; i++){
	  // loop over each DOF on the edge
	  if (edgeIndex[i] == 0){
		  for (int j = 0 ; j < nrDOF_edge_ ; j++ , offs++){
			  tmp[offs] = refAllx[2+j] * refAlly[edgeMultI[i]];
		  }
	  }
	  else
	  {
		  for (int j = 0 ; j < nrDOF_edge_ ; j++ , offs++){
			  tmp[offs] = refAllx[edgeMultI[i]] * refAlly[2+j];
		  }
	  }
  }

  // loop over all internal DOFs
  for (int i=0 ; i < nrDOF_face_ ; i++ , offs++){
	  tmp[offs] = refAllx[2+i] * refAlly[2+(nrDOF_face_-1-i)];
  }

  // compute the results
  for (int i=0; i<result.length(); i++)
  {
    if (deriv.order()==0) result[i] = tmp[i].value();
    else 
      result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
  }
  //SUNDANCE_OUT( true , "Legendre::evalOnQuad result.length():" << result.length() );
}

void Legendre::evalOnBrick(const Point& pt,
  const MultiIndex& deriv,
  Array<double>& result) const
{
  ADReal x = ADReal(pt[0], 0, 3);
  ADReal y = ADReal(pt[1], 1, 3);
  ADReal z = ADReal(pt[2], 2, 3);
  ADReal one(1.0, 3);
  
  TEST_FOR_EXCEPT(true);
}