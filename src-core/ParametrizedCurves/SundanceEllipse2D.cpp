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

#include "SundanceEllipse2D.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"

using namespace Sundance;

Ellipse2D::Ellipse2D(double px, double py, double a, double b, double a1, double a2) :
	CurveBase(1, a1, a2), px_(px), py_(px), a_(a), b_(b)
{
}

Ellipse2D::~Ellipse2D()
{
}

Expr Ellipse2D::getParams() const
{
	// return the parameters of the box
	// the ellipse equation is ((x-px)/a)^2 + ((y-py)/b)^2 - 1 = 0
	return Expr(List(px_, py_, a_ , b_ ));
}

double Ellipse2D::curveEquation(const Point& evalPoint) const
{
	int verb = 0;
	TEST_FOR_EXCEPTION(evalPoint.dim() != 2, RuntimeError,
			"Ellipse2D::curveEquation() evaluation point dimension must be 2");

	// calculate the distance compared to the middle point
	// the ellipse equation is ((x-px)/a)^2 + ((y-py)/b)^2 - 1 = 0
	double distX =  (evalPoint[0] - px_)*(evalPoint[0] - px_)/(a_*a_);
	double distY =  (evalPoint[1] - py_)*(evalPoint[1] - py_)/(b_*b_);
	distX = distX + distY - 1.0;
	SUNDANCE_OUT(verb > 3, " Ellipse2D::curveEquation for:" << evalPoint << " is: " << distX);
	return distX;
}

void Ellipse2D::returnIntersect(const Point& start, const Point& end, int& nrPoints, Array<
		double>& result) const
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


void Ellipse2D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
    // set initial values
	int verb = 0;
	nrPoints = 0;
	result.resize(2);
	Point tmpP(0.0,0.0);
	SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints STARTS , start:" << start << " , end:" << end );
	// we have two main cases
	if ( fabs(end[0] - start[0]) < 1e-9){
	    // we have a vertical line, where x is constant
        double term = (b_/a_)* (end[0] - px_);
         // (y - py_)^2 = (b-(b/a)*(x-px))*(b+(b/a)*(x-px))
        double sqrt_term = ( b_ - term) * ( b_ + term);

		SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints CASE I , sqrt_term =" << sqrt_term );
        if (sqrt_term < 0.0) return;
        sqrt_term = sqrt(sqrt_term);

        tmpP[0] = end[0]; tmpP[1] = py_ - sqrt_term;
		SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints CASE I , P1 =" << tmpP );
        if ( (tmpP[0] <= end[0]) && ( tmpP[0] >= start[0]) && (tmpP[1] <= end[1]) && (tmpP[1] >= start[1])){
        	result[nrPoints++] = tmpP;
        	SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints  ADD1 tmpP:" << tmpP << " nrPoints:" << nrPoints);
        }

        if (sqrt_term > 1e-10){
            tmpP[0] = end[0]; tmpP[1] = py_ + sqrt_term;
    		SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints CASE I , P2 =" << tmpP );
            if ( (tmpP[0] <= end[0]) && ( tmpP[0] >= start[0]) && (tmpP[1] <= end[1]) && (tmpP[1] >= start[1])){
            	result[nrPoints++] = tmpP;
            	SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints  ADD2 tmpP:" << tmpP << " nrPoints:" << nrPoints);
            }
        }
	}
	else{
		double alpha = (end[1] - start[1])/(end[0] - start[0]);
		double beta = start[1] - alpha*start[0];
		// now set up the complicated second oreder equation
		// a X^X + b X + C = 0 , for the intersection points
		double a2 = a_*a_;
		double b2 = b_*b_;
		double a = b2 + alpha*alpha* a2;
		double b = -2*px_*b2 + 2*alpha*(beta-py_)*a2;
		double c = a2*b2*(px_*px_/a2 + ((beta-py_)*(beta-py_)/b2) - 1);

		double det = b * b - 4 * a * c;

		SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints CASE II , det:" << det );

		if (det < 0.0)
			return; // no solution
		det = sqrt(det);

		// numerically more favorable than usual usage of solution formula
		int sign = b < 0.0 ? -1 : 1;
		double q = -0.5 * (b + sign * det);

		// first solution
		double x = c / q; //(-b - det) / (2*a); //c / q;
		double y = alpha*x + beta;
		SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints CASE II , P1 x=" << x << " y=" << y );

		// Only consider intersection points between 'start' and 'end'
		if ( (x <= end[0]) && ( x >= start[0]) && (y <= end[1]) && ( y >= start[1])){
		   tmpP[0] = x; tmpP[1] = y;
		   result[nrPoints++] = tmpP;
	        SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints  ADD1 tmpP:" << tmpP << " nrPoints:" << nrPoints);
		}

		// second solution (only considered if different from first)
		if (det < 1.e-14)
			return;

		x = q / a; //(-b + det) / (2*a); //q / a;
		y = alpha*x + beta;
		SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints CASE II , P2 x=" << x << " y=" << y );

		if ((x <= end[0]) && ( x >= start[0]) && (y <= end[1]) && ( y >= start[1]))
		{
			tmpP[0] = x; tmpP[1] = y;
			// Sorting: first solution is nearest to 'start'
			if ( (nrPoints > 0) && (x < result[0][0]) )
			{
				result[nrPoints++] = result[0];
				result[0] = tmpP;
		        SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints  ADD21 tmpP:" << tmpP << " nrPoints:" << nrPoints);
			}
			else
			{
				result[nrPoints++] = tmpP;
		        SUNDANCE_OUT(verb > 3, " Ellipse2D::returnIntersectPoints  ADD22 tmpP:" << tmpP << " nrPoints:" << nrPoints);
			}
		}
	}
}
