/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CLOSEDNEWTONCOTES_H
#define SUNDANCE_CLOSEDNEWTONCOTES_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Teuchos;
  using namespace Internal;

  /** 
   * Closed Newton-Cotes quadrature rules for lines and triangles
   */
  class ClosedNewtonCotes : public QuadratureFamilyBase
  {
  public:
    /** */
    ClosedNewtonCotes();

    /** */
    virtual ~ClosedNewtonCotes(){;}

    /** */
    virtual XMLObject toXML() const ;

    /* handleable boilerplate */
    GET_RCP(QuadratureFamilyStub);

  protected:
    /** compute a rule for the reference line cell */
    virtual void getLineRule(Array<Point>& quadPoints,
                             Array<double>& quadWeights) const ;

    /** compute a rule for the reference triangle cell */
    virtual void getTriangleRule(Array<Point>& quadPoints,
                                 Array<double>& quadWeights) const ;

  };
}


#endif
