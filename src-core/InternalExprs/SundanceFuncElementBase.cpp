/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFuncElementBase.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

FuncElementBase::FuncElementBase(const string& name)
	: ScalarExpr(), name_(name), id_(nextID()++)
{}

ostream& FuncElementBase::toText(ostream& os, bool /* paren */) const 
{
	os << name_;
	return os;
}

ostream& FuncElementBase::toLatex(ostream& os, bool /* paren */) const 
{
	os << name_;
	return os;
}



  
