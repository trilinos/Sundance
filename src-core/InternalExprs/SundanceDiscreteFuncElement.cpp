/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

DiscreteFuncElement::DiscreteFuncElement(DiscreteFunctionStub* master, 
                                         const string& name,
                                         int myIndex)
	: LeafExpr(), 
    FuncElementBase(name),
    master_(master),
    myIndex_(myIndex)
{}

bool DiscreteFuncElement::hasNonzeroDeriv(const MultipleDeriv& md) const 
{
  TimeMonitor t(nonzeroDerivCheckTimer());

  if (md.order()==0) return true;

  MultipleDeriv::const_iterator iter;
  
  for (iter=md.begin(); iter != md.end(); iter++)
    {
      const Deriv& d = *iter;
      if (!d.isCoordDeriv()) return false;
    }
  return true;
}

XMLObject DiscreteFuncElement::toXML() const 
{
	XMLObject rtn("DiscreteFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}

