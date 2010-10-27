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

#include "SundancePlane3D.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"

using namespace Sundance;

Plane3D::Plane3D(double a, double b, double c, double a1, double a2) :
	CurveBase(2, a1, a2), a_(a) , b_(b) , c_(c)
{
}

Plane3D::~Plane3D()
{
}

Expr Plane3D::getParams() const
{
	// return the parameters of the 3D plane z = a*x + b*y + c
	return Expr(List( a_ , b_ , c_ ));
}

double Plane3D::curveEquation(const Point& evalPoint) const
{
	TEST_FOR_EXCEPTION(evalPoint.dim() != 3, RuntimeError,
			"Plane3D::curveEquation() evaluation point dimension must be 3");

	// z = a*x + b*y + c
	return ( a_*evalPoint[0] + b_*evalPoint[1] + c_ -  evalPoint[2]);
}

void Plane3D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
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

void Plane3D::returnIntersect(const Point& start, const Point& end, int& nrPoints, Array<
		double>& result) const
{
	 // we calculate two different planes which have the same slope , and take the free term
	// z = a*x + b*y + c
    double c1 = start[2] - a_* start[0] - b_* start[1];
    double c2 = end[2] - a_* end[0] - b_* end[1];

    nrPoints = 0;

    // test if c_ is between c1 and c2
    if ( ((c_ >= c1) && (c_ <= c2)) || ((c_ <= c1) && (c_ >= c2)) ){
    	result.resize(1);
    	result[0] = fabs((c1-c_)/(c1-c2));
    	nrPoints = 1;
    }
}

