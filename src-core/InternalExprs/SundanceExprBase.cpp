/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExprBase.hpp"
#include <typeinfo>

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

ExprBase::ExprBase()
{;}

string ExprBase::typeName() const 
{
	return typeid(*this).name();
}

string ExprBase::toString() const 
{
	TeuchosOStringStream ss;
	toText(ss, false);
	return TEUCHOS_OSTRINGSTREAM_GET_C_STR(ss);
}
