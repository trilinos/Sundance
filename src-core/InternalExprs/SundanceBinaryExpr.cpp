/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBinaryExpr.hpp"
#include "SundanceExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;


BinaryExpr::BinaryExpr(const RefCountPtr<ScalarExpr>& left,
                       const RefCountPtr<ScalarExpr>& right, int sign)
	: ExprWithChildren(tuple(left, right)), 
    sign_(sign)
{}





ostream& BinaryExpr:: toText(ostream& os, bool paren) const 
{
	if (Expr::showAllParens() || (paren && parenthesizeSelf())) os << "(";
	leftScalar()->toText(os, parenthesizeOperands()) ;
	os	 << opChar();
	if (leftScalar()->isHungryDiffOp())
		{
			rightScalar()->toText(os, true);
		}
	else
		{
			rightScalar()->toText(os, parenthesizeOperands());
		}
	if (Expr::showAllParens() || (paren && parenthesizeSelf())) os << ")";

	return os;
}

ostream& BinaryExpr:: toLatex(ostream& os, bool paren) const 
{
	if (Expr::showAllParens() || (paren && parenthesizeSelf())) os << "(";
	leftScalar()->toLatex(os, parenthesizeOperands()) ;
	os	 << opChar();

	if (leftScalar()->isHungryDiffOp())
		{
			rightScalar()->toText(os, true);
		}
	else
		{
			rightScalar()->toText(os, parenthesizeOperands());
		}
	if (Expr::showAllParens() || (paren && parenthesizeSelf())) os << ")";

	return os;
}

XMLObject BinaryExpr::toXML() const 
{
	XMLObject rtn(xmlTag());
	rtn.addChild(left().toXML());
	rtn.addChild(right().toXML());

	return rtn;
}


