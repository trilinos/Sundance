/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_LISTFUNCTION_H
#define SUNDANCE_LISTFUNCTION_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFunc.hpp"



namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  
  namespace Internal

  using std::string;
  using std::ostream;

  /** 
   * ListFunction 
   */
  template <class FuncElemType>
  class ListFunction : public ListExpr
  {
  public:
    /** Construct a list function with a name and number of elements */
    ListFunction(const string& name, int nElems)
      : ListExpr()
    {
      for (int i=0; i<nElems; i++)
        {
          string elemName = name + "[" + Teuchos::toString(i) + "]";
          append(new FuncElemType(this, elemName, i));
        }
    }

    /** virtual dtor */
    virtual ~ListFunction() {;}

  };
}
}


#endif
