/* @HEADER@ */
/* @HEADER@ */

#include "SundanceWeakFormBatch.hpp"
#include "SundanceOut.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceTabs.hpp"
#include "SundanceFuncWithBasis.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceTestFunction.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;

WeakFormBatch::WeakFormBatch(const MultipleDeriv& functionalDeriv, 
                             int derivIndex)
  : isOneForm_(),
    testBasis_(),
    unkBasis_(),
    miTest_(),
    miUnk_(),
    testID_(),
    unkID_(),
    derivIndex_(),
    trans_()
{
  int testID;
  int unkID;

  extractWeakForm(functionalDeriv, testBasis_, unkBasis_,
                  miTest_, miUnk_, testID, unkID, isOneForm_);

  testID_.append(testID);
  unkID_.append(unkID);
  derivIndex_.append(derivIndex);
  trans_.append(0);
}

bool WeakFormBatch::tryToAdd(const MultipleDeriv& functionalDeriv, 
                             int derivIndex)
{
  bool isOneForm;
      
  BasisFamily testBasis;

  BasisFamily unkBasis;

  MultiIndex miTest;

  MultiIndex miUnk;

  int testID;

  int unkID;

  extractWeakForm(functionalDeriv, testBasis, unkBasis,
                  miTest, miUnk, testID, unkID, isOneForm);

  if (!(isOneForm_ == isOneForm)) return false;

  if (isOneForm)
    {
      if (testBasis_ == testBasis && miTest_==miTest)
        {
          trans_.append(false);
          testID_.append(testID);
          derivIndex_.append(derivIndex);
          return true;
        }
      else
        {
          return false;
        }
    }
  else
    {
      if (testBasis_ == testBasis && miTest_==miTest
          && unkBasis_==unkBasis && miUnk_==miUnk)
        {
          trans_.append(false);
          testID_.append(testID);
          unkID_.append(unkID);
          derivIndex_.append(derivIndex);
          return true;
        }
      else if (testBasis_ == unkBasis && miTest_==miUnk
          && unkBasis_==testBasis && miUnk_==miTest)
        {
          trans_.append(true);
          testID_.append(testID);
          unkID_.append(unkID);
          derivIndex_.append(derivIndex);
          return true;
        }
      else
        {
          return false;
        }
    }
      
}


void WeakFormBatch::extractWeakForm(const MultipleDeriv& functionalDeriv,
                                    BasisFamily& testBasis, 
                                    BasisFamily& unkBasis,
                                    MultiIndex& miTest, MultiIndex& miUnk,
                                    int& testID, int& unkID, 
                                    bool& isOneForm) const
{
  MultipleDeriv::const_iterator iter;

  isOneForm = false;  

  if (functionalDeriv.size()==0) return;
  TEST_FOR_EXCEPTION(functionalDeriv.size()>2, InternalError,
                     "WeakFormBatch::extractWeakForm detected a functional "
                     "derivative of order > 2: " 
                     << functionalDeriv.toString());

  bool foundTest = false;
  bool foundUnk = false;

  for (iter = functionalDeriv.begin(); iter != functionalDeriv.end(); iter++)
    {
      const Deriv& d = *iter;
      
      TEST_FOR_EXCEPTION(!d.isFunctionalDeriv(), InternalError,
                         "WeakFormBatch::extractWeakForm "
                         "detected a non-functional derivative: "
                         << functionalDeriv.toString());
      
      const FunctionalDeriv* f = d.funcDeriv();
      
      const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(f->func());
      if (u != 0)
        {
          foundUnk = true;
          unkID = f->funcID();
          const UnknownFunction* uf 
            = dynamic_cast<const  UnknownFunction*>(u->master());
          TEST_FOR_EXCEPTION(uf==0, InternalError, 
                             "WeakFormBatch::extractWeakForm failed to cast "
                             "UnknownFunctionStub to UnknownFunction");
          const FuncWithBasis* fb 
            = dynamic_cast<const  FuncWithBasis*>(uf);
          TEST_FOR_EXCEPTION(fb==0, InternalError, 
                             "WeakFormBatch::extractWeakForm failed to cast "
                             "UnknownFunction to funcWithBasis");
          int myIndex = u->myIndex();
          unkBasis = fb->basis()[myIndex];
          miUnk = f->multiIndex();
          continue;
        }
      
      const TestFuncElement* t 
        = dynamic_cast<const TestFuncElement*>(f->func());
      if (t != 0)
        {
          foundTest = true;
          testID = f->funcID();
          
          const TestFunction* tf 
            = dynamic_cast<const  TestFunction*>(t->master());
          TEST_FOR_EXCEPTION(tf==0, InternalError, 
                             "WeakFormBatch::extractWeakForm failed to cast "
                             "TestFunctionStub to TestFunction");
          const FuncWithBasis* fb 
            = dynamic_cast<const  FuncWithBasis*>(tf);
          TEST_FOR_EXCEPTION(fb==0, InternalError, 
                             "WeakFormBatch::extractWeakForm failed to cast "
                             "TestFunction to funcWithBasis");

          int myIndex = t->myIndex();
          testBasis = fb->basis()[myIndex];
          miTest = f->multiIndex();
          continue;
        }
    }

  if (!foundUnk) isOneForm = true;
}

void WeakFormBatch::print(ostream& os) const 
{
  os << "WeakFormBatch" << endl;
  {
    Tabs tab1;
    os << tab1 << "positions in results array: " << derivIndex_ << endl;
    os << tab1 << "test basis: " << testBasis_ << endl;
    os << tab1 << "test multiindex: " << miTest_.toString() << endl;
    os << tab1 << "test funcID: " << testID_ << endl;
  }
  if (!isOneForm_)
    {
      Tabs tab1;
      os << tab1 << "unk basis: " << unkBasis_ << endl;
      os << tab1 << "unk multiindex: " << miUnk_.toString() << endl;
      os << tab1 << "unk funcID: " << unkID_ << endl;
      os << tab1 << "is transposed: " << trans_ << endl;
    }
}
