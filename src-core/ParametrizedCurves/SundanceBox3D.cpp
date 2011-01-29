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

#include "SundanceBox3D.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"

using namespace Sundance;

Box3D::Box3D(double px, double py, double pz,
		     double ox, double oy, double oz, double a1, double a2, bool flipD ) :
	CurveBase(2, a1, a2,flipD), px_(px), py_(py), pz_(pz), ox_(ox), oy_(oy), oz_(oz)
{
}

Box3D::~Box3D()
{
}

Expr Box3D::getParams() const
{
	// return the parameters of the box
	return Expr(List(px_, py_, py_, ox_ , oy_, oz_));
}

double Box3D::curveEquation(const Point& evalPoint) const
{
	int verb = 0;
	TEST_FOR_EXCEPTION(evalPoint.dim() != 3, RuntimeError,
			"Box3D::curveEquation() evaluation point dimension must be 3");

	// calculate the distance compared to the middle point
	double distX =  fabs(px_ + 0.5*ox_ - evalPoint[0]) - 0.5*ox_;
	double distY =  fabs(py_ + 0.5*oy_ - evalPoint[1]) - 0.5*oy_;
	double distZ =  fabs(pz_ + 0.5*oz_ - evalPoint[2]) - 0.5*oz_;
	distX = (distX > distY) ? distX : distY ;
	distX = (distX > distZ) ? distX : distZ ;

	SUNDANCE_OUT(verb > 3, " Box3D::curveEquation for:" << evalPoint << " is: " << distX);

	return flipDomains_*distX;
}

void Box3D::returnIntersect(const Point& start, const Point& end, int& nrPoints,
		Array<double>& result) const
{

	Array<Point> t;
	returnIntersectPoints(start, end, nrPoints, t);

	result.resize(nrPoints);

	// Return coordinates instead of t values
	for (int i = 0; i < nrPoints; i++)
	{
		Point tmp( end[0]-t[i][0]/(end[0]-start[0]) , end[1]-t[i][1]/(end[1]-start[1]) );
		result[i] =  sqrt( tmp*tmp );
	}
}

void Box3D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
	TEST_FOR_EXCEPTION(true, RuntimeError,
				"Box3D::returnIntersectPoints() not implemented yet");
	// todo: implement this
}

