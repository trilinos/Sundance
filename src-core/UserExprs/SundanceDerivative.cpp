/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDerivative.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

Derivative::Derivative(int direction)
	: ScalarExpr(), m_()
{
	m_[direction] = 1;
}

ostream& Derivative::toText(ostream& os, bool /* paren */) const 
{
	os << "D[" << m_.toString() << "]";
	return os;
}

ostream& Derivative::toLatex(ostream& os, bool /* paren */) const 
{
	os << "D^{" << m_.toString() << "}";
	return os;
}

XMLObject Derivative::toXML() const 
{
	XMLObject rtn("Derivative");
	rtn.addAttribute("m", m_.toString());
	return rtn;
}


