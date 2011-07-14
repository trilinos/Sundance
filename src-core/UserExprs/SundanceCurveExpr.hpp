/*
 * SundanceCurveExpr.hpp
 *
 *  Created on: Jun 7, 2011
 *      Author: benk
 */

#ifndef SUNDANCECURVEEXPR_HPP_
#define SUNDANCECURVEEXPR_HPP_

#include "SundancePointwiseUserDefFunctor.hpp"
#include "SundanceParametrizedCurve.hpp"

namespace Sundance{

/** class which serves as expression for interface (curve or surface) values */
class CurveExpr : public PointwiseUserDefFunctor0
{
public:
    /** Ctor the input is a ParamCurve and an index for the scalar field */
	CurveExpr(const ParametrizedCurve& interface , int scalarFieldIndex)
    : PointwiseUserDefFunctor0("CurveExpr", interface.getCurveDim()+1 , 1)
      , interface_(interface) , scalarFieldIndex_(scalarFieldIndex){;}

    /** to evaluate the expression just pass the control to the interface object */
    void eval0(const double* vars, double* f) const {
    	interface_.eval0( vars, f , scalarFieldIndex_ ) ;
    }

private:
  const ParametrizedCurve& interface_;
  int scalarFieldIndex_;
};

//Expr CurveExpression(const ParametrizedCurve& interface , int scalarFieldIndex, const Expr& x, const Expr& y)
//{
//  return new UserDefOp( List(x,y), rcp(new CurveExpr(interface, scalarFieldIndex)));
//}

}

#endif /* SUNDANCECURVEEXPR_HPP_ */
