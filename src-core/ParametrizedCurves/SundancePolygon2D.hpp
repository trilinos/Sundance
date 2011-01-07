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

#ifndef SUNDANCEPOLYGON2D_H_
#define SUNDANCEPOLYGON2D_H_

#include "SundanceCurveBase.hpp"
#include "SundanceParametrizedCurve.hpp"

namespace Sundance
{

/**  */
class Polygon2D: public Sundance::CurveBase
{
public:

	/** Ctor with the points which form the polygon (the direction of the points matters for IN & OUT )*/
	Polygon2D(const Mesh& mesh , const Array<Point>& points , double a1 , double a2);

	/** Ctor to read in the polygon from file from file */
	Polygon2D(const Mesh& mesh , const std::string& filename , double a1 , double a2);

	/** */
	virtual ~Polygon2D() {;}

	/** @return Expr The parameters of the curve which uniquely defines the curve*/
	virtual Expr getParams() const;

	/**
	 * This function should be implemented
	 * @param evalPoint point where the polygon equation is evaluated <br>
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

	/** Writes the geometry into a VTK file for visualization purposes
	 * @param filename */
	virtual void writeToVTK(const std::string& filename) const;

	/** generate the polygon which is the unification (not intersection) of two polygons */
    static RCP<CurveBase> unite(ParametrizedCurve& c1 , ParametrizedCurve& c2);

private:

	/** looks for each point in which cell it is contained (later for evaluation purposes)*/
	void computeMaxCellLIDs();

	/** the mesh will be needed for point location, since the geometry is
	 * non-conform with respect to the mesh */
	const Mesh& mesh_;

	/** The points in the specified order which form */
	Array<Point> polyPoints_;

	/** each points should know in which cell is evaluable (DiscreteSpace)*/
	Array<int> pointsMaxCellLID_;

	/** the index of the line which is at last the intersection point */
	static int intersectionEdge_;

};

}
#endif /* SUNDANCEPOLYGON2D_H_ */
