/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BASISFAMILY_H
#define SUNDANCE_BASISFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceBasisFamilyBase.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdFwk::Internal;
  using namespace TSFExtended;

  /** */
  class BasisFamily : public Handle<BasisFamilyBase>
    {
    public:
      /* */
      HANDLE_CTORS(BasisFamily, BasisFamilyBase);

      /** */
      XMLObject toXML() const ;

      /** */
      int order() const ;
    };
}

#endif
