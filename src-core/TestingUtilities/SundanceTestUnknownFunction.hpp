/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TESTUNKNOWNFUNCTION_H
#define SUNDANCE_TESTUNKNOWNFUNCTION_H

#include "SundanceDefs.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceADField.hpp"

namespace SundanceTesting
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** 
   *
   */
  class TestUnknownFunction : public UnknownFunctionStub
  {
  public:
    /** */
    TestUnknownFunction(const ADField& field, const string& name="")
      : UnknownFunctionStub(name, 1), field_(field)
    {;}
    
    /** virtual destructor */
    virtual ~TestUnknownFunction() {;}

    /** */
    Expr createDiscreteFunction() const ;

    /* boilerplate */
    GET_RCP(ExprBase);

  private:
    ADField field_;
  };

}



#endif
