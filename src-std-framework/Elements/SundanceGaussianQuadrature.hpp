/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_GAUSSIANQUADRATURE_H
#define SUNDANCE_GAUSSIANQUADRATURE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamilyStub.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Teuchos;

  /** */
  class GaussianQuadrature : public QuadratureFamilyStub 
  {
  public:
    /** */
    GaussianQuadrature(int order);

    /** */
    virtual ~GaussianQuadrature(){;}

    /** */
    virtual XMLObject toXML() const ;

    /* handleable boilerplate */
    GET_RCP(QuadratureFamilyStub);
      
  private:
    int order_;
      
  };
}


#endif
