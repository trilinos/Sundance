/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TESTDISCRETEFUNCTION_H
#define SUNDANCE_TESTDISCRETEFUNCTION_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceADField.hpp"
#include "TSFVector.hpp"

namespace SundanceTesting
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** 
   *
   */
  class TestDiscreteFunction : public DiscreteFunctionStub
  {
  public:
    /** */
    TestDiscreteFunction(const ADField& field, 
                         const string& name="")
      :  DiscreteFunctionStub(name, 1), field_(field)
    {;}
   
    /** virtual destructor */
    virtual ~TestDiscreteFunction() {;}

    /** */
    ADReal evaluate() const ;

    /** */
    const ADField& field() const {return field_;}

    /* boilerplate */
    GET_RCP(ExprBase);


  private:
    ADField field_;
  };

}



#endif
