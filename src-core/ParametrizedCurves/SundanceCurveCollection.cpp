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

#include "SundanceCurveCollection.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"

using namespace Sundance;

CurveCollection::CurveCollection(double a1, double a2,int dim) :
	CurveBase(dim, a1, a2, false) , nrCurves_(0)
{
}

CurveCollection::~CurveCollection()
{
}

Expr CurveCollection::getParams() const
{
	// no parameters
	return Expr(List(0.0));
}

double CurveCollection::curveEquation_intern(const Point& evalPoint) const
{
	int verb = 0;

	// calculate the distance compared to the middle point
	double distX = 1e+100;

	SUNDANCE_OUT(verb > 3, " CurveCollection::curveEquation before:" << evalPoint << " is: " << distX);
	if ( nrCurves_ > 0) distX = curves_[0]->curveEquation(evalPoint);
	for (int c = 1 ; c < nrCurves_; c++){
		double tmp = curves_[c]->curveEquation(evalPoint);
		distX = (tmp > distX) ? distX : tmp;
	}

	SUNDANCE_OUT(verb > 3, " CurveCollection::curveEquation for:" << evalPoint << " is: " << distX);

	return distX;
}

void CurveCollection::returnIntersect(const Point& start, const Point& end, int& nrPoints,
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

void CurveCollection::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
    int verb = 0;

    result.resize(0); nrPoints = 0;
	for (int c = 0 ; c < nrCurves_; c++){
		Array<Point> pnts;
		int nrP;
		curves_[c]->returnIntersectPoints( start, end, nrP, pnts);
		for (int p = 0 ; p < nrP ; p++){
			result.append( pnts[p] );
			nrPoints++;
		}
	}

	SUNDANCE_OUT(verb > 3, " CurveCollection::returnIntersectPoints END , nrPoints:" << nrPoints
			<< " , start:" << start << " , end:" << end);
}
