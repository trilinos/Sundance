/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ADBASIS_H
#define SUNDANCE_ADBASIS_H

#include "SundanceDefs.hpp"
#include "SundanceADReal.hpp"

namespace SundanceTesting
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** */
  class ADBasis
  {
  public:
    /** */
    ADBasis(int order = 0);

    /** */
    ADReal evaluate(const Point& x) const ;
    
  private:
    int m_;
  };

}



#endif
