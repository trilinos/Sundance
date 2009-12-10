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

#include "SundanceLagrange.hpp"
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


Lagrange::Lagrange(int order)
  : order_(order)
{
TEST_FOR_EXCEPTION(order < 0, RuntimeError,
                     "invalid polynomial order=" << order
                     << " in Lagrange ctor");
}

bool Lagrange::supportsCellTypePair(
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
    case QuadCell:
        switch(cellType){
          case QuadCell:
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

void Lagrange::print(std::ostream& os) const 
{
  os << "Lagrange(" << order_ << ")";
}

int Lagrange::nReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(cellType)
    {
    case PointCell:
      return 1;
    case LineCell:
      return 1 + order_;
    case TriangleCell:
      {
        switch(order_)
          {
          case 0:
            return 1;
          case 1:
            return 3;
          case 2:
            return 6;
          case 3:
            return 10;
          default:
            TEST_FOR_EXCEPTION(true, RuntimeError, "order=" << order_ 
                               << " not implemented in Lagrange basis");
            return -1; // -Wall
          }
      }
    case QuadCell:{
        switch(order_) {
          case 0:
            return 1;
          case 1:
            return 4;
          case 2:
            return 8; // was 9 before
                      // !!! ---- change here also if one changes the number of unknowns ----
          default:
            TEST_FOR_EXCEPTION(true, RuntimeError, "order=" << order_
                               << " not implemented in Lagrange basis");
            return -1; // -Wall
          }
      }
    case TetCell:
      {switch(order_)
          {
          case 0:
            return 1;
          case 1:
            return 4;
          case 2:
            return 10;
          default:
            TEST_FOR_EXCEPTION(true, RuntimeError, "order=" << order_ 
                               << " not implemented in Lagrange basis");
            return -1; // -Wall
          }

      }
    case BrickCell:
      {switch(order_)
          {
          case 0:
            return 1;
          case 1:
            return 8;
          case 2:
            return 20;
          default:
            TEST_FOR_EXCEPTION(true, RuntimeError, "order=" << order_
                               << " not implemented in Lagrange basis");
            return -1; // -Wall
          }

      }
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in Lagrange basis");
      return -1; // -Wall
    }
}

void Lagrange::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  typedef Array<int> Aint;
  switch(cellType)
    {
    case PointCell:
      dofs.resize(1);
      dofs[0] = tuple<Aint>(tuple(0));
      return;
    case LineCell:
      dofs.resize(2);
      dofs[0] = tuple<Aint>(tuple(0), tuple(1));
      dofs[1] = tuple<Aint>(makeRange(2, order()));
      return;
    case TriangleCell:
      {
        int n = order()-1;
        dofs.resize(3);
        dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2));
        dofs[1] = tuple<Aint>(makeRange(3,2+n), 
                        makeRange(3+n, 2+2*n),
                        makeRange(3+2*n, 2+3*n));
        if (order()<3)
          {
            dofs[2] = tuple(Array<int>());
          }
        else 
          {
            dofs[2] = tuple<Aint>(makeRange(3+3*n, 3+3*n));
          }
        return;
      }
    case QuadCell:{ // TODO: ask Kevin about Q0 (P0)
         dofs.resize(3);
    	 // dofs[0] are DOFs at Points
    	 dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2), tuple(3));
    	 if (order() == 2)
    	 {
    	   // dofs[1] are the DOFs on the line
    	   dofs[1] = tuple<Aint>(tuple(4), tuple(5), tuple(6), tuple(7));
    	   //dofs[2] are DOFs inside the Cell

    	   dofs[2] = tuple(Array<int>());
           //dofs[2] = tuple<Aint>(tuple(8));
    	 } else {
    	   dofs[1] = tuple(Array<int>(), Array<int>(),
    		                 Array<int>(), Array<int>());
           dofs[2] = tuple(Array<int>());
    	 }
         return;
      }
    case TetCell:
      {
        dofs.resize(4);
        dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2), tuple(3));
        if (order() == 2)
          {
            dofs[1] = tuple<Aint>(tuple(4), tuple(5), tuple(6), 
                            tuple(7), tuple(8), tuple(9));
          }
        else
          {
            dofs[1] = tuple(Array<int>(), Array<int>(),
                            Array<int>(), Array<int>(), 
                            Array<int>(), Array<int>());
          }
        dofs[2] = tuple(Array<int>(), Array<int>(), 
                        Array<int>(), Array<int>());
        dofs[3] = tuple(Array<int>());
        return;
      }
    case BrickCell:{
         dofs.resize(4);
    	 // dofs[0] are DOFs at Points
    	 dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2), tuple(3),
    			               tuple(4), tuple(5), tuple(6), tuple(7));
    	 if (order() == 2)
    	 {
    	   // dofs[1] are the DOFs on the line
    	   dofs[1] = tuple<Aint>( tuple(8), tuple(9) , tuple(10), tuple(11),
    			                 tuple(12), tuple(13), tuple(14), tuple(15),
    			                 tuple(16), tuple(17), tuple(18), tuple(19));
    	   //dofs[2] are DOFs on the faces
    	 } else {
    	   dofs[1] = tuple(Array<int>(), Array<int>(), Array<int>(), Array<int>(),
    			           Array<int>(), Array<int>(), Array<int>(), Array<int>(),
    			           Array<int>(), Array<int>(), Array<int>(), Array<int>());
    	 }
         dofs[2] = tuple(Array<int>(), Array<int>(), Array<int>(), Array<int>(),
                         Array<int>(), Array<int>() );
         dofs[3] = tuple(Array<int>());
         return;
      }
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in Lagrange basis");
    }
}



Array<int> Lagrange::makeRange(int low, int high)
{
  if (high < low) return Array<int>();

  Array<int> rtn(high-low+1);
  for (int i=0; i<rtn.length(); i++) rtn[i] = low+i;
  return rtn;
}

void Lagrange::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    RuntimeError,
    "cannot evaluate spatial derivative " << sds << " on Lagrange basis");
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
    case QuadCell:
      for (int i=0; i<pts.length(); i++)
       {
          evalOnquad(pts[i], deriv, result[0][i]);
       }
       return;
    case TetCell:
      for (int i=0; i<pts.length(); i++)
        {
          evalOnTet(pts[i], deriv, result[0][i]);
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
                         "Lagrange::refEval() unimplemented for cell type "
                         << cellType);

    }
}

/* ---------- evaluation on different cell types -------------- */

void Lagrange::evalOnLine(const Point& pt, 
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
  /*
	switch(order_)
		{
		case 0:
			tmp[0] = one;
			break;
		case 1:
			tmp[0] = 1.0 - x;
			tmp[1] = x;
			break;
		case 2:
			tmp[0] = 2.0*(1.0-x)*(0.5-x);
			tmp[1] = 2.0*x*(x-0.5);
			tmp[2] = 4.0*x*(1.0-x);
			break;
    case 3:
      tmp[0] = 9.0/2.0 * (1.0 - x) * (1.0/3.0 - x) * (2.0/3.0 - x);
      tmp[1] = 9.0/2.0 * x * (x - 1.0/3.0) * (x - 2.0/3.0);
      tmp[2] = 27.0/2.0 * x * (1.0 - x) * (2.0/3.0 - x);
      tmp[3] = 27.0/2.0 * x * (1.0 - x) * (x - 1.0/3.0);

      break;
		default:
			SUNDANCE_ERROR("Lagrange::evalOnLine polynomial order > 2 has not been"
                     " implemented");
		}
  */
	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else result[i] = tmp[i].gradient()[0];
		}
}

void Lagrange::evalOnTriangle(const Point& pt, 
															const MultiIndex& deriv,
															Array<double>& result) const



{
	ADReal x = ADReal(pt[0], 0, 2);
	ADReal y = ADReal(pt[1], 1, 2);
	ADReal one(1.0, 2);

  Array<ADReal> tmp;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
               << y.value());

	switch(order_)
		{
		case 0:
      result.resize(1);
      tmp.resize(1);
			tmp[0] = one;
			break;
		case 1:
      result.resize(3);
      tmp.resize(3);
			tmp[0] = 1.0 - x - y;
			tmp[1] = x;
			tmp[2] = y;
			break;
		case 2:
      result.resize(6);
      tmp.resize(6);
			tmp[0] = (1.0-x-y)*(1.0-2.0*x-2.0*y);
			tmp[1] = 2.0*x*(x-0.5);
			tmp[2] = 2.0*y*(y-0.5);
			tmp[3] = 4.0*x*(1.0-x-y); 
			tmp[4] = 4.0*x*y;
			tmp[5] = 4.0*y*(1.0-x-y);
			break;
    case 3:
      result.resize(10);
      tmp.resize(10);
      {
        ADReal q1 = 1.0 - x - y;
        ADReal q2 = x;
        ADReal q3 = y;
        tmp[0] = 9.0/2.0 * q1 * (q1 - 2.0/3.0) * (q1 - 1.0/3.0);
        tmp[1] = 9.0/2.0 * q2 * (q2 - 2.0/3.0) * (q2 - 1.0/3.0);
        tmp[2] = 9.0/2.0 * q3 * (q3 - 2.0/3.0) * (q3 - 1.0/3.0);
        tmp[3] = 27.0/2.0 * q1 * q2   * (q1 - 1.0/3.0);
        tmp[4] = 27.0/2.0 * q1 * q2   * (q2 - 1.0/3.0);
        tmp[5] = 27.0/2.0 * q2 * q3   * (q2 - 1.0/3.0);
        tmp[6] = 27.0/2.0 * q2 * q3   * (q3 - 1.0/3.0);
        tmp[7] = 27.0/2.0 * q3 * q1   * (q3 - 1.0/3.0);
        tmp[8] = 27.0/2.0 * q3 * q1   * (q1 - 1.0/3.0);
        tmp[9] = 27.0 * q1 * q2 * q3;
      }
      break;
		default:
#ifndef TRILINOS_7
			SUNDANCE_ERROR("Lagrange::evalOnTriangle poly order > 2");
#else
			SUNDANCE_ERROR7("Lagrange::evalOnTriangle poly order > 2");
#endif
		}

	for (int i=0; i<tmp.length(); i++)
		{
      SUNDANCE_OUT(this->verb() > 3,
                   "tmp[" << i << "]=" << tmp[i].value() 
                   << " grad=" << tmp[i].gradient());
			if (deriv.order()==0) result[i] = tmp[i].value();
			else 
          result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}

void Lagrange::evalOnquad(const Point& pt,
	    const MultiIndex& deriv,
	    Array<double>& result) const {

	ADReal x = ADReal(pt[0], 0, 2);
	ADReal y = ADReal(pt[1], 1, 2);
	ADReal one(1.0, 2);

  Array<ADReal> tmp;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
               << y.value());

	switch(order_)
		{
		case 0:
            result.resize(1);
            tmp.resize(1);
			tmp[0] = one;
			break;
		case 1:
            result.resize(4);
            tmp.resize(4);
			tmp[0] = (1.0 - x)*(1.0 - y);
			tmp[1] = x*(1.0 - y);
			tmp[2] = (1.0 - x)*y;
			tmp[3] = x*y; //These are the 4 basis functions for the first order
			break;
		case 2: //The basis functions calculated with Octave / Matlab
			/*
            result.resize(9);
            tmp.resize(9);
                   //9  -9  -9   9  -6   0  -6   0   4
			tmp[0] = 9.0 - 9.0*x - 9.0*y + 9.0*x*y - 6.0*(x*x -1)*(y-1)
					- 6.0*(x-1)*(y*y-1)                     + 4.0*(x*x-1)*(y*y-1);
			       //6  -3  -6   3  -6   0  -3   1   4
			tmp[1] = 6.0 - 3.0*x - 6.0*y + 3.0*x*y - 6.0*(x*x -1)*(y-1)
                    - 3.0*(x-1)*(y*y-1) + 1.0*(x+1)*(y*y-1) + 4.0*(x*x-1)*(y*y-1);
			       //6  -6  -3   3  -3   1  -6   0   4
			tmp[2] = 6.0 - 6.0*x - 3.0*y + 3.0*x*y - 3.0*(x*x -1)*(y-1) + 1.0*(x*x -1)*(y+1)
					-6.0*(x-1)*(y*y-1) +                    + 4.0*(x*x-1)*(y*y-1);
			       //4  -2  -2   1  -3   1  -3   1   4
			tmp[3] = 4.0 - 2.0*x - 2.0*y + 1.0*x*y - 3.0*(x*x -1)*(y-1) + 1.0*(x*x -1)*(y+1)
					- 3.0*(x-1)*(y*y-1) + 1.0*(x+1)*(y*y-1) + 4.0*(x*x-1)*(y*y-1);
			       //-12   12   12  -12   12    0    8    0   -8
			tmp[4] = -12.0 + 12.0*x + 12.0*y - 12.0*x*y + 12.0*(x*x -1)*(y-1)
					+ 8.0*(x-1)*(y*y-1)                     - 8.0*(x*x-1)*(y*y-1);
			       //-12   12   12  -12    8    0   12    0   -8
			tmp[5] = -12.0 + 12.0*x + 12.0*y - 12.0*x*y + 8.0*(x*x -1)*(y-1)
					+ 12.0*(x-1)*(y*y-1)                    - 8.0*(x*x-1)*(y*y-1);
			       //-8   4   8  -4   8   0   6  -2  -8
			tmp[6] = -8.0 + 4.0*x + 8.0*y - 4.0*x*y + 8.0*(x*x -1)*(y-1)
					+ 6.0*(x-1)*(y*y-1) - 2.0*(x+1)*(y*y-1) - 8.0*(x*x-1)*(y*y-1);
			       //-8   8   4  -4   6  -2   8   0  -8
			tmp[7] = -8.0 + 8.0*x + 4.0*y - 4.0*x*y + 6.0*(x*x -1)*(y-1) - 2.0*(x*x -1)*(y+1)
					+ 8.0*(x-1)*(y*y-1)                     - 8.0*(x*x-1)*(y*y-1);
			       //16  -16  -16   16  -16   -0  -16   -0   16
			tmp[8] = 16.0 - 16.0*x - 16.0*y + 16.0*x*y - 16.0*(x*x -1)*(y-1)
					- 16.0*(x-1)*(y*y-1)                    + 16.0*(x*x-1)*(y*y-1);*/
	        result.resize(8);
	        tmp.resize(8);
	               //5  -5  -5   5  -2   0  -2   0
			tmp[0] = 5.0 - 5.0*x - 5.0*y + 5.0*x*y - 2.0*(x*x -1)*(y-1)
					- 2.0*(x-1)*(y*y-1);
			       //2   1  -2  -1  -2   0   1   1
			tmp[1] = 2.0 + 1.0*x - 2.0*y - 1.0*x*y - 2.0*(x*x -1)*(y-1)
	                + 1.0*(x-1)*(y*y-1) + 1.0*(x+1)*(y*y-1);
			       //2  -2   1  -1   1   1  -2   0
			tmp[2] = 2.0 - 2.0*x + 1.0*y - 1.0*x*y + 1.0*(x*x -1)*(y-1) + 1.0*(x*x -1)*(y+1)
					- 2.0*(x-1)*(y*y-1);
			       //0   2   2  -3   1   1   1   1
			tmp[3] =      2.0*x + 2.0*y - 3.0*x*y + 1.0*(x*x -1)*(y-1) + 1.0*(x*x -1)*(y+1)
					+ 1.0*(x-1)*(y*y-1)  + 1.0*(x+1)*(y*y-1);
			       //-4   4   4  -4   4   0   0   0
			tmp[4] = -4.0 + 4.0*x + 4.0*y - 4.0*x*y + 4.0*(x*x -1)*(y-1);

			       //-4   4   4  -4  -0   0   4   0
			tmp[5] = -4.0 + 4.0*x + 4.0*y - 4.0*x*y
					+ 4.0*(x-1)*(y*y-1);
			       //0  -4   0   4  -0  -0  -2  -2
			tmp[6] =      - 4.0*x         + 4.0*x*y
					- 2.0*(x-1)*(y*y-1)  - 2.0*(x+1)*(y*y-1);
			       //0   0  -4   4  -2  -2   0   0
			tmp[7] =             - 4.0*y + 4.0*x*y - 2.0*(x*x -1)*(y-1) - 2.0*(x*x -1)*(y+1);
			break;
		default:{}
		}

	 for (int i=0; i<tmp.length(); i++) {
		 SUNDANCE_OUT(this->verb() > 3 ,
                   "tmp[" << i << "]=" << tmp[i].value()
                   << " grad=" << tmp[i].gradient());
		if (deriv.order()==0)
				result[i] = tmp[i].value();
		else
                result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
	}
}

void Lagrange::evalOnTet(const Point& pt, 
												 const MultiIndex& deriv,
												 Array<double>& result) const
{
	ADReal x = ADReal(pt[0], 0, 3);
	ADReal y = ADReal(pt[1], 1, 3);
	ADReal z = ADReal(pt[2], 2, 3);
	ADReal one(1.0, 3);


	Array<ADReal> tmp(result.length());

	switch(order_)
		{
		case 0:
      tmp.resize(1);
      result.resize(1);
			tmp[0] = one;
			break;
		case 1:
      result.resize(4);
      tmp.resize(4);
			tmp[0] = 1.0 - x - y - z;
			tmp[1] = x;
			tmp[2] = y;
			tmp[3] = z;
			break;
		case 2:
      result.resize(10);
      tmp.resize(10);
			tmp[0] = (1.0-x-y-z)*(1.0-2.0*x-2.0*y-2.0*z);
			tmp[1] = 2.0*x*(x-0.5);
			tmp[2] = 2.0*y*(y-0.5);
			tmp[3] = 2.0*z*(z-0.5);
			tmp[4] = 4.0*x*(1.0-x-y-z);
			tmp[5] = 4.0*x*y;
			tmp[6] = 4.0*y*(1.0-x-y-z);
			tmp[7] = 4.0*z*(1.0-x-y-z);
			tmp[8] = 4.0*x*z;
			tmp[9] = 4.0*y*z;
			break;
		default:
#ifndef TRILINOS_7
			SUNDANCE_ERROR("Lagrange::evalOnTet poly order > 2");
#else
			SUNDANCE_ERROR7("Lagrange::evalOnTet poly order > 2");
#endif
		}

	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else 
				result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}

void Lagrange::evalOnBrick(const Point& pt,
												 const MultiIndex& deriv,
												 Array<double>& result) const
{
	ADReal x = ADReal(pt[0], 0, 3);
	ADReal y = ADReal(pt[1], 1, 3);
	ADReal z = ADReal(pt[2], 2, 3);
	ADReal one(1.0, 3);


	Array<ADReal> tmp(result.length());

	switch(order_)
		{
		case 0:
            tmp.resize(1);
            result.resize(1);
			tmp[0] = one;
			break;
		case 1: // tri linear basis functions for bricks
            result.resize(8);
            tmp.resize(8);
			tmp[0] = (1.0 - x)*(1.0 - y)*(1.0 - z);
			tmp[1] = (x)*(1.0 - y)*(1.0 - z);
			tmp[2] = (1.0 - x)*(y)*(1.0 - z);
			tmp[3] = (x)*(y)*(1.0 - z);
			tmp[4] = (1.0 - x)*(1.0 - y)*(z);
			tmp[5] = (x)*(1.0 - y)*(z);
			tmp[6] = (1.0 - x)*(y)*(z);
			tmp[7] = (x)*(y)*(z);
			break;
		case 2: // second order function in 3D
	        result.resize(20);
	        tmp.resize(20);

			// general form of the Lagrange polynom

			// 1.0 + 1.0*x + 1.0*y + 1.0*z + 1.0*x*y + 1.0*x*z + 1.0*y*z + 1.0*x*y*z
		   //+ 1.0*(x*x - 1)*(y-1)*(z-1) + 1.0*(x*x - 1)*(y-1)*(z+1) + 1.0*(x*x - 1)*(y+1)*(z-1) + 1.0*(x*x - 1)*(y+1)*(z+1)
		   //+ 1.0*(y*y - 1)*(x-1)*(z-1) + 1.0*(y*y - 1)*(x-1)*(z+1) + 1.0*(y*y - 1)*(x+1)*(z-1) + 1.0*(y*y - 1)*(x+1)*(z+1)
		   //+ 1.0*(z*z - 1)*(x-1)*(y-1) + 1.0*(z*z - 1)*(x-1)*(y+1) + 1.0*(z*z - 1)*(x+1)*(y-1) + 1.0*(z*z - 1)*(x+1)*(y+1);

             // 7  -7  -7  -7   7   7   7  -7
			tmp[0] = 7.0 - 7.0*x - 7.0*y - 7.0*z + 7.0*x*y + 7.0*x*z + 7.0*y*z - 7.0*x*y*z
			// 2   0   0   0
				   + 2.0*(x*x - 1)*(y-1)*(z-1)
		    // 2   0   0   0
				   + 2.0*(y*y - 1)*(x-1)*(z-1)
		    // 2   0   0   0
				   + 2.0*(z*z - 1)*(x-1)*(y-1);

			// 2   3  -2  -2  -3  -3   2   3
			tmp[1] = 2.0 + 3.0*x - 2.0*y - 2.0*z - 3.0*x*y - 3.0*x*z + 2.0*y*z + 3.0*x*y*z
            //2   0   0   0
					+ 2.0*(x*x - 1)*(y-1)*(z-1)
		    //-1   0  -1   0
					- 1.0*(y*y - 1)*(x-1)*(z-1)                             - 1.0*(y*y - 1)*(x+1)*(z-1)
			//-1   0  -1    0
                    - 1.0*(z*z - 1)*(x-1)*(y-1)                             - 1.0*(z*z - 1)*(x+1)*(y-1);

            // 2  -2   3  -2  -3   2  -3   3
			tmp[2] = 2.0 - 2.0*x + 3.0*y - 2.0*z - 3.0*x*y + 2.0*x*z - 3.0*y*z + 3.0*x*y*z
				   // -1   0  -1   0
				   - 1.0*(x*x - 1)*(y-1)*(z-1)                             - 1.0*(x*x - 1)*(y+1)*(z-1)
				   //  2   0   0   0
				   + 2.0*(y*y - 1)*(x-1)*(z-1)
				   // -1  -1   0   0
				   - 1.0*(z*z - 1)*(x-1)*(y-1) - 1.0*(z*z - 1)*(x-1)*(y+1);

			 //0   2   2   0  -1  -2  -2  1
			tmp[3] =       2.0*x + 2.0*y       - 1.0*x*y - 2.0*x*z - 2.0*y*z + 1.0*x*y*z
					//-1   0  -1   0
					- 1.0*(x*x - 1)*(y-1)*(z-1)                             - 1.0*(x*x - 1)*(y+1)*(z-1)
                    //-1   0  -1   0
					- 1.0*(y*y - 1)*(x-1)*(z-1)                             - 1.0*(y*y - 1)*(x+1)*(z-1)
					//0.5   0.5   0.5   0.5
					+ 0.5*(z*z - 1)*(x-1)*(y-1) + 0.5*(z*z - 1)*(x-1)*(y+1) + 0.5*(z*z - 1)*(x+1)*(y-1) + 0.5*(z*z - 1)*(x+1)*(y+1);

			// 2  -2  -2   3   2  -3  -3   3
			tmp[4] = 2.0 - 2.0*x - 2.0*y + 3.0*z + 2.0*x*y - 3.0*x*z - 3.0*y*z + 3.0*x*y*z
					 //-1  -1   0   0
					 - 1.0*(x*x - 1)*(y-1)*(z-1) - 1.0*(x*x - 1)*(y-1)*(z+1)
					 //-1  -1   0   0
					 - 1.0*(y*y - 1)*(x-1)*(z-1) - 1.0*(y*y - 1)*(x-1)*(z+1)
					 //2   0   0   0
					 + 2.0*(z*z - 1)*(x-1)*(y-1);

            // 0   2   0   2  -2  -1  -2  1
			tmp[5] =      + 2.0*x +      + 2.0*z - 2.0*x*y - 1.0*x*z - 2.0*y*z + 1.0*x*y*z
					//-1  -1   0  0
					- 1.0*(x*x - 1)*(y-1)*(z-1) - 1.0*(x*x - 1)*(y-1)*(z+1)
					//0.5   0.5  0.5  0.5
					+ 0.5*(y*y - 1)*(x-1)*(z-1) + 0.5*(y*y - 1)*(x-1)*(z+1) + 0.5*(y*y - 1)*(x+1)*(z-1) + 0.5*(y*y - 1)*(x+1)*(z+1)
					// -1   0  -1   0
					- 1.0*(z*z - 1)*(x-1)*(y-1)                             - 1.0*(z*z - 1)*(x+1)*(y-1);

			//0   0   2   2  -2  -2  -1   1
			tmp[6] =            2.0*y + 2.0*z - 2.0*x*y - 2.0*x*z - 1.0*y*z + 1.0*x*y*z
				   //0.5   0.5   0.5   0.5
			       + 0.5*(x*x - 1)*(y-1)*(z-1) + 0.5*(x*x - 1)*(y-1)*(z+1) + 0.5*(x*x - 1)*(y+1)*(z-1) + 0.5*(x*x - 1)*(y+1)*(z+1)
			       //-1  -1   0   0
				   - 1.0*(y*y - 1)*(x-1)*(z-1) - 1.0*(y*y - 1)*(x-1)*(z+1)
				   //-1  -1   0   0
				   - 1.0*(z*z - 1)*(x-1)*(y-1) - 1.0*(z*z - 1)*(x-1)*(y+1);

			//   0   0   0   0   2   2   2  -5
			tmp[7] =                             2.0*x*y + 2.0*x*z + 2.0*y*z - 5.0*x*y*z
				   //0.5   0.5   0.5   0.5
				   + 0.5*(x*x - 1)*(y-1)*(z-1) + 0.5*(x*x - 1)*(y-1)*(z+1) + 0.5*(x*x - 1)*(y+1)*(z-1) + 0.5*(x*x - 1)*(y+1)*(z+1)
				   //0.5   0.5  0.5   0.5
				   + 0.5*(y*y - 1)*(x-1)*(z-1) + 0.5*(y*y - 1)*(x-1)*(z+1) + 0.5*(y*y - 1)*(x+1)*(z-1) + 0.5*(y*y - 1)*(x+1)*(z+1)
				   //0.5   0.5   0.5   0.5
				   + 0.5*(z*z - 1)*(x-1)*(y-1) + 0.5*(z*z - 1)*(x-1)*(y+1) + 0.5*(z*z - 1)*(x+1)*(y-1) + 0.5*(z*z - 1)*(x+1)*(y+1);

            // -4   4   4   4  -4  -4  -4   4
			tmp[8] = -4.0 + 4.0*x + 4.0*y + 4.0*z - 4.0*x*y - 4.0*x*z - 4.0*y*z + 4.0*x*y*z
				   //-4   0   0   0
				   - 4.0*(x*x - 1)*(y-1)*(z-1);
				   //0   0   0   0
				   //0   0   0   0

			// -4   4   4   4  -4  -4  -4   4
			tmp[9] = -4.0 + 4.0*x + 4.0*y + 4.0*z - 4.0*x*y - 4.0*x*z - 4.0*y*z + 4.0*x*y*z
				   // 0   0   0   0
				   //-4   0   0   0
				   - 4.0*(y*y - 1)*(x-1)*(z-1);
				   // 0   0   0   0

			//  -4   4   4   4  -4  -4  -4   4
			tmp[10] = -4.0 + 4.0*x + 4.0*y + 4.0*z - 4.0*x*y - 4.0*x*z - 4.0*y*z + 4.0*x*y*z
					//-0   0   0  -0
				    //-0   0   0  -0
					//-4   0   0   0
					- 4.0*(z*z - 1)*(x-1)*(y-1);

            // 0  -4   0   0   4   4   0  -4
			tmp[11] =     - 4.0*x           + 4.0*x*y + 4.0*x*z          - 4.0*x*y*z
					 //0   0   0   0
					 //2   0   2  -0
					 + 2.0*(y*y - 1)*(x-1)*(z-1)                             + 2.0*(y*y - 1)*(x+1)*(z-1);
					 // 0   0   0   0

            //0  -4   0   0   4   4   0  -4
			tmp[12] =     - 4.0*x             + 4.0*x*y + 4.0*x*z         - 4.0*x*y*z
					 //0   0   0  -0
					 //-0   0   0  -0
					 //2   0   2   0
					 + 2.0*(z*z - 1)*(x-1)*(y-1)                         + 2.0*(z*z - 1)*(x+1)*(y-1);

			//0   0  -4   0   4   0   4  -4
			tmp[13] =            - 4.0*y       + 4.0*x*y        + 4.0*y*z - 4.0*x*y*z
					// 2   0   2   0
					+ 2.0*(x*x - 1)*(y-1)*(z-1)                      + 2.0*(x*x - 1)*(y+1)*(z-1);
					// 0   0   0   0
					// 0   0   0   0

			//0  -0  -4   0   4   0   4  -4
			tmp[14] =              - 4.0*y         + 4.0*x*y        + 4.0*y*z - 4.0*x*y*z
					 //-0  -0   0  -0
					 //-0   0   0  -0
					 // 2   2   0   0
					 + 2.0*(z*z - 1)*(x-1)*(y-1) + 2.0*(z*z - 1)*(x-1)*(y+1);
			//0   0   0   0  -4   0   0   4
			tmp[15] =                            - 4.0*x*y                    + 4.0*x*y*z
					//-0   0   0  -0
					//-0   0   0  -0
					//-1  -1  -1  -1
					- 1.0*(z*z - 1)*(x-1)*(y-1) - 1.0*(z*z - 1)*(x-1)*(y+1) - 1.0*(z*z - 1)*(x+1)*(y-1) - 1.0*(z*z - 1)*(x+1)*(y+1);
			//0   0   0  -4   0   4   4  -4
			tmp[16] =                 - 4.0*z       + 4.0*x*z + 4.0*y*z - 4.0*x*y*z
					// 2   2   0   0
					+ 2.0*(x*x - 1)*(y-1)*(z-1) + 2.0*(x*x - 1)*(y-1)*(z+1);
					// 0   0   0   0
					// 0   0   0   0

			//0   0   0  -4   0   4   4  -4
			tmp[17] =             - 4.0*z        + 4.0*x*z + 4.0*y*z - 4.0*x*y*z
					//  0   0   0   0
					//  2   2   0  -0
					+ 2.0*(y*y - 1)*(x-1)*(z-1) + 2.0*(y*y - 1)*(x-1)*(z+1);
					//  0   0   0   0

			//0   0   0  -0   0  -4   0   4
			tmp[18] =                                        - 4.0*x*z         + 4.0*x*y*z
					//0   0  -0   0
					//-1  -1  -1  -1
					- 1.0*(y*y - 1)*(x-1)*(z-1) - 1.0*(y*y - 1)*(x-1)*(z+1) - 1.0*(y*y - 1)*(x+1)*(z-1) - 1.0*(y*y - 1)*(x+1)*(z+1);
					//0   0   0   0

			// 0   0   0   0   0   0  -4   4
			tmp[19] =                                                - 4.0*y*z + 4.0*x*y*z
					//-1  -1  -1  -1
					- 1.0*(x*x - 1)*(y-1)*(z-1) - 1.0*(x*x - 1)*(y-1)*(z+1) - 1.0*(x*x - 1)*(y+1)*(z-1) - 1.0*(x*x - 1)*(y+1)*(z+1);
					//0   0   0   0
					//0   0   0   0

			break;
		default:
#ifndef TRILINOS_7
			SUNDANCE_ERROR("Lagrange::evalOnBrick poly order > 2");
#else
			SUNDANCE_ERROR7("Lagrange::evalOnBrick poly order > 2");
#endif
		}

	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else
				result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}



