/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUnaryMinus.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

UnaryMinus::UnaryMinus(const RefCountPtr<ScalarExpr>& arg)
  : UnaryExpr(arg)
{
  typedef Set<int>::const_iterator setIter;
  if (isEvaluatable(arg.get()))
    {
      for (int d=0; d<MultiIndex::maxDim(); d++) 
        {
          setOrderOfDependency(d, evaluatableArg()->orderOfSpatialDependency(d));
        }

      setFuncIDSet(evaluatableArg()->funcIDSet());
      for (setIter i=funcIDSet().begin(); i != funcIDSet().end(); i++)
        {
          int argOD = evaluatableArg()->orderOfFunctionalDependency(*i);
          setOrderOfFunctionalDependency(*i, argOD);
        }
    }
}

ostream& UnaryMinus::toText(ostream& os, bool paren) const 
{
  if (paren) os << "(";
  os << "-" << arg().toString();
  if (paren) os << ")";
  return os;
}

ostream& UnaryMinus::toLatex(ostream& os, bool paren) const 
{
  return toText(os, paren);
}

XMLObject UnaryMinus::toXML() const
{
  XMLObject rtn("UnaryMinus");
  rtn.addChild(arg().toXML());
  return rtn;
}
