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

#ifndef SUNDANCEELLIPSE2D_H_
#define SUNDANCEELLIPSE2D_H_

#include "SundanceCurveBase.hpp"

namespace Sundance
{

/** This class models the Ellipse in 2D.<br>
 * the ellipse equation is ((x-px)/a)^2 + ((y-py)/b)^2 - 1 = 0*/
class Ellipse2D: public Sundance::CurveBase
{
public:

	/** The constructor of a 2D Ellipse
	 * @param px center X coordinate
	 * @param py center Y coordinate
	 * @param a
	 * @param c
	 * @param a1 alpha1 coefficient for FCM
	 * @param a2 alpha2 coefficient for FCM */
	Ellipse2D(double px, double py, double a, double b, double a1, double a2);

	virtual ~Ellipse2D();

	/** @return Expr The parameters of the curve which uniquely defines the curve*/
	virtual Expr getParams() const;

	/**
	 * This function should be implemented
	 * @param evalPoint point where the box equation is evaluated <br>
	 * @return double the value of the curve equation at the evaluation point  */
	virtual double curveEquation(const Point& evalPoint) const;

	/**
	 * This function is important for nonstructural mesh integration.<br>
	 * The function returns the intersection points with a given line (only 2D at the moment)<br>
	 * The line is defined with a starting point and an end point<br>
	 * The resulting points should be between the two points
	 * @param start, the start point of the line
	 * @param end, the end point of the line
	 * @param nrPoints , number of resulted (intersected) points
	 * @param result , the resulted points of intersections */
	virtual void returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
			Array<Point>& result) const;

	/**
	 * As above, but, instead of coordinates, returns intersection values t in [0,1]
	 * and the line is defined by "start + t * (end-start)"
	 */
	virtual void returnIntersect(const Point& start, const Point& end, int& nrPoints,
			Array<double>& result) const;

	/** Return a ref count pointer to self */
	virtual RCP<CurveBase> getRcp()
	{
		return rcp(this);
	}

	/** This curve is a real curve */
	virtual bool isCurveValid() const
	{
		return true;
	}

	/** In the case of simple geometries the geometry it can be transformed to a polygon, which
	 * reflects the original geometry
	 * @param mesh
	 * @param resolution , the global resolution */
	virtual const RCP<CurveBase> getPolygon(const Mesh& mesh ,double resolution) const ;

private:

	/** */
	double px_;

	/** */
	double py_;

	/** */
	double a_;

	/** */
	double b_;

};

}
#endif /* SUNDANCEELLIPSE2D_H_ */
