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
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;



void P1NC::print(ostream& os) const 
{
  os << "P1NC()";
}

int P1NC::nNodes(int spatialDim, 
                 const CellType& cellType) const
{
  switch(cellType)
    {
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

void P1NC::getLocalDOFs(const CellType& cellType,
                            Array<Array<Array<int> > >& dofs) const 
{
  switch(cellType)
    {
    case LineCell:
      dofs.resize(2);
      dofs[0] = tuple(tuple(0), tuple(1));
      dofs[1] = tuple(makeRange(2, order()));
      return;
    case TriangleCell:
      int n = order()-1;
      dofs.resize(3);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple(tuple(0), tuple(1), tuple(2));
      dofs[2] = tuple(Array<int>());
      return;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in P1NC basis");
    }
}



void P1NC::refEval(int spatialDim, 
                   const CellType& cellType,
                   const Array<Point>& pts,
                   const MultiIndex& deriv,
                   Array<Array<double> >& result) const
{
  result.resize(pts.length());

  switch(cellType)
    {
    case LineCell:
      for (int i=0; i<pts.length(); i++)
        {
          evalOnLine(pts[i], deriv, result[i]);
        }
      return;
    case TriangleCell:
      for (int i=0; i<pts.length(); i++)
        {
          evalOnTriangle(pts[i], deriv, result[i]);
        }
      return;
    default:
      SUNDANCE_ERROR("P1NC::refEval() unimplemented for cell type "
                     << cellType);

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

  SUNDANCE_OUT(this->verbosity() > VerbHigh, "x=" << x.value() << " y="
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
			SUNDANCE_ERROR("P1NC::evalOnTriangle poly order > 2");
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


void P1NC::evalOnTet(const Point& pt, 
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
			SUNDANCE_ERROR("P1NC::evalOnTet poly order > 2");
		}

	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else 
				result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}
