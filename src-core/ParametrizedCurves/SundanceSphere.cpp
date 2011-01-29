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

#include "SundanceSphere.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"

using namespace Sundance;

Sphere::Sphere(double centerx, double centery, double centerz, double radius, double a1,
		double a2 , bool flipD ) :
	CurveBase(2, a1, a2, flipD), _centerx(centerx), _centery(centery), _centerz(centerz), _radius(radius)
{
}

Sphere::~Sphere()
{
}

Expr Sphere::getParams() const
{
	// return the parameters of the 2D sphere
	return Expr(List(_centerx, _centery, _centerz, _radius));
}

double Sphere::curveEquation(const Point& evalPoint) const
{
	TEST_FOR_EXCEPTION(evalPoint.dim() != 3, RuntimeError,
			"Sphere::curveEquation() evaluation point dimension must be 3");

	Point center( _centerx , _centery , _centerz );

	// the sphere equation is (x-cx)^2 + (y-cy)^2 + (z-cz)^2 - r^2 = 0
	return flipDomains_*((evalPoint - center) * (evalPoint - center)) - _radius * _radius;
}

void Sphere::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
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

void Sphere::returnIntersect(const Point& start, const Point& end, int& nrPoints, Array<
		double>& result) const
{
	// direction vector; line == start + t * vec
	Point vec = end - start;

	// center point of sphere
	Point center( _centerx, _centery , _centerz );

	nrPoints = 0;
	result.resize(2);

	// intersection of sphere and line leads to
	// norm(start + t * vec - center) == r
	// i.e. following coefficients of quadratic equation
	double a = vec * vec;
	double b = 2 * ((start - center) * vec);
	double c = ((start - center) * (start - center)) - _radius * _radius;

	double det = b * b - 4 * a * c;
	if (det < 0.0)
		return; // no solution
	det = sqrt(det);

	// numerically more favorable than usual usage of solution formula
	int sign = b < 0.0 ? -1 : 1;
	double q = -0.5 * (b + sign * det);

	// first solution
	double t = c / q;

	// Only consider intersection points between 'start' and 'end'
	if (t >= 0.0 && t <= 1.0)
		result[nrPoints++] = t;

	// second solution (only considered if different from first)
	if (det < 1.e-14)
		return;
	t = q / a;
	if (t >= 0.0 && t <= 1.0)
	{
		// Sorting: first solution is nearest to 'start'
		if (t < result[0])
		{
			result[nrPoints++] = result[0];
			result[0] = t;
		}
		else
		{
			result[nrPoints++] = t;
		}
	}
}

