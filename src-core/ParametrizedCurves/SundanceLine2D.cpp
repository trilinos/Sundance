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

#include "SundanceLine2D.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"
#include "SundancePolygon2D.hpp"

using namespace Sundance;

Line2D::Line2D(double slope, double b, double a1, double a2, bool flipD ) :
	CurveBase(1, a1, a2, flipD), slope_(slope), b_(b)
{
}

Line2D::~Line2D()
{
}

Expr Line2D::getParams() const
{
	// return the parameters of the 2D circle
	return Expr(List(slope_, b_ ));
}

double Line2D::curveEquation(const Point& evalPoint) const
{
	TEST_FOR_EXCEPTION(evalPoint.dim() != 2, RuntimeError,
			"Line2D::curveEquation() evaluation point dimension must be 2");

	// y = slope *x + b
	return flipDomains_*(slope_ * evalPoint[0] + b_ - evalPoint[1]);
}

void Line2D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
	Array<double> t;
	returnIntersect(start, end, nrPoints, t);

	result.resize(2);

	// Return coordinates instead of t values
	for (int i = 0; i < nrPoints; i++)
	{
		result[i] = start + t[i] * (end - start);
	}
}

void Line2D::returnIntersect(const Point& start, const Point& end, int& nrPoints, Array<
		double>& result) const
{
	 // we calculate two different lines which have the same slope , and take the free term
     double b1 = start[1] - slope_* start[0];
     double b2 = end[1] - slope_* end[0];

     nrPoints = 0;

     // test if b_ is between b1 and b2
     if ( ((b_ >= b1) && (b_ <= b2)) || ((b_ <= b1) && (b_ >= b2)) ){
    	 result.resize(1);
    	 result[0] = fabs((b1-b_)/(b1-b2));
    	 nrPoints = 1;
     }
}

const RCP<CurveBase> Line2D::getPolygon(const Mesh& mesh , double resolution) const {

	int verb = 0;

	/* y = slope*x + b */

	// dummy implementation
	// have always 100 points, and x goes from 0 to 1
	int nrPoints = 100;
	SUNDANCE_MSG3( verb , " Line2D::getPolygon nrPoints = " << nrPoints );
	Array<Point> points(nrPoints);


	for (int pI = 0 ; pI < nrPoints ; pI++){
		Point p( (-5.0 +10.0*(double)(pI)/(double)(nrPoints-1)) , slope_*(-5.0 +10.0*(double)(pI)/(double)(nrPoints-1)) + b_);
		SUNDANCE_MSG3( verb , " Line2D::getPolygon add pint " << pI << " p=" << p );
		points[pI] = p;
	}

	// return the polygon
	return rcp(new Polygon2D( mesh , points , _alpha1 , _alpha2 , (flipDomains_<0)));
}

