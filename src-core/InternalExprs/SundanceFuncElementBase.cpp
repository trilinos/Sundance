/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFuncElementBase.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

FuncElementBase::FuncElementBase(const string& rootName,
                                 const string& suffix)
	: ScalarExpr(), name_(rootName + suffix), rootName_(rootName),
    suffix_(suffix), id_(nextID()++)
{}

FuncElementBase::FuncElementBase(const string& rootName)
	: ScalarExpr(), name_(rootName), rootName_(rootName),
    suffix_(), id_(nextID()++)
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



  
