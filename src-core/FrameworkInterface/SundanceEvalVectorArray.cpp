/* @HEADER@ */
/* @HEADER@ */


#include "SundanceEvalVectorArray.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

EvalVectorArray::EvalVectorArray()
  : Array<RefCountPtr<EvalVector> >()
{;}

void EvalVectorArray::copy(const RefCountPtr<EvalVectorArray>& other)
{
  resize(other->size());

  for (int i=0; i<other->size(); i++)
    {
      (*this)[i]->copy((*other)[i]);
    }
}

void EvalVectorArray::steal(const RefCountPtr<EvalVectorArray>& other)
{
  resize(other->size());

  for (int i=0; i<other->size(); i++)
    {
      (*this)[i] = (*other)[i];
    }
}

ostream& EvalVectorArray::print(ostream& os, const DerivSet& derivs) const
{
  Tabs tab;
  TEST_FOR_EXCEPTION(derivs.size() != size(),
                     InternalError,
                     "mismatch between deriv set size=" << derivs.size()
                     << "and result vector size " << size()
                     << "in EvalVectorArray::print");

  DerivSet::const_iterator iter;

  int maxlen = 25;
  for (iter=derivs.begin(); iter != derivs.end(); iter++)
    {
      int s = (*iter).toString().length();
      if (s > maxlen) maxlen = s;
    }
  
  int i = 0;
  for (iter=derivs.begin(); iter != derivs.end(); i++, iter++)
    {
      os << tab;
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << (*iter).toString() << "\t\t";
      (*this)[i]->print(os);
      os << endl;
    }
  return os;
}
