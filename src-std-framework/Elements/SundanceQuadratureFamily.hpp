/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_QUADRATUREFAMILY_H
#define SUNDANCE_QUADRATUREFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdFwk
{
  using namespace TSFExtended;
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** */
  class QuadratureFamily : public Handle<QuadratureFamilyStub>
  {
  public:
    /* */
    HANDLE_CTORS(QuadratureFamily, QuadratureFamilyStub);
    /** */
    XMLObject toXML() const ;

    /** */
    int order() const ;
  private:
  };
}

#endif
