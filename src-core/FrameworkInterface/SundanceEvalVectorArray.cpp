/* @HEADER@ */
/* @HEADER@ */


#include "SundanceEvalVectorArray.hpp"
#include "SundanceSparsitySuperset.hpp"
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

ostream& EvalVectorArray::print(ostream& os, 
                                const SparsitySuperset* derivs) const
{
  Tabs tab;
  TEST_FOR_EXCEPTION(derivs->numDerivs() != size(),
                     InternalError,
                     "mismatch between deriv set size=" << derivs->numDerivs()
                     << "and result vector size " << size()
                     << "in EvalVectorArray::print");

  int maxlen = 25;
  for (int i=0; i<derivs->numDerivs(); i++)
    {
      int s = derivs->deriv(i).toString().length();
      if (s > maxlen) maxlen = s;
    }
  
  for (int i=0; i<derivs->numDerivs(); i++)
    {
      os << tab;
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << derivs->deriv(i).toString() << "\t\t";
      (*this)[i]->print(os);
      os << endl;
    }
  return os;
}
