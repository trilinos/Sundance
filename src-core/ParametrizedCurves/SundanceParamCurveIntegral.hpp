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
/*
 * SundanceParamCurveIntegral.hpp
 *
 *  Created on: Sep 1, 2010
 *      Author: benk
 */

#ifndef SUNDANCEPARAMCURVEINTEGRAL_HPP_
#define SUNDANCEPARAMCURVEINTEGRAL_HPP_

#include "SundanceCurveBase.hpp"

namespace Sundance {

   // forward declaration
   class ParametrizedCurve;

/** class to wrap one other class */
class ParamCurveIntegral : public Sundance::CurveBase
{

public:

	/** The constructor */
	ParamCurveIntegral(ParametrizedCurve& paramcurve);

	virtual ~ParamCurveIntegral() {;}

	/** See upper class */
	virtual Expr getParams() const { return paramcurve_.getParams();}

	/** return the control points of the parameterized curve */
	virtual Array<Point>& getControlPoints() { return paramcurve_.getControlPoints(); }

	/** updates the state of the curve if the control parameters are changed */
	virtual void update() { return paramcurve_.update(); }

protected:
	/** See upper class */
	virtual double curveEquation_intern(const Point& evaluationPoint) const { return paramcurve_.curveEquation(evaluationPoint); }
public:

	/** See upper class */
	virtual void returnIntersectPoints(const Point& startEdgePoint, const Point& endEdgePoint,
			                      int& nrPoints ,Array<Point>& result) const {
		return paramcurve_.returnIntersectPoints( startEdgePoint , endEdgePoint , nrPoints , result );
	}

	/** See upper class */
	virtual void returnIntersect(const Point& startEdgePoint, const Point& endEdgePoint,
				                      int& nrPoints ,Array<double>& result) const {
		return paramcurve_.returnIntersect( startEdgePoint , endEdgePoint , nrPoints , result);
	}

	/** See upper class */
	virtual bool isCurveValid() const {
		return paramcurve_.isCurveValid();
	}

	/** This is the only place where this function should return true */
	virtual bool isCurveIntegral() const { return true; }

	/** Returns the underlying object */
	virtual const CurveBase* getUnderlyingCurveObj() const { return paramcurve_.ptr().get() ; }

	/** Return a ref count pointer to self */
	virtual RCP<CurveBase> getRcp()
	{
		return rcp(this);
	}

private:


	/** Parameterized curve for which this class is a wrapper and so this class signals
	 * that this will be a curve integral */
	ParametrizedCurve paramcurve_;
};

}

#endif /* SUNDANCEPARAMCURVEINTEGRAL_HPP_ */
