/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_GAUSSIANQUADRATURE_H
#define SUNDANCE_GAUSSIANQUADRATURE_H

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
   * Family of optimal Gaussian integration rules, e.g., Gauss-Legendre on 
   * lines, Dunavant on triangles. 
   */
  class GaussianQuadrature : public QuadratureFamilyBase
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

  protected:
    /** compute a rule for the reference line cell */
    virtual void getLineRule(Array<Point>& quadPoints,
                             Array<double>& quadWeights) const ;

    /** compute a rule for the reference triangle cell */
    virtual void getTriangleRule(Array<Point>& quadPoints,
                                 Array<double>& quadWeights) const ;

    /** compute a rule for the reference tet cell */
    virtual void getTetRule(Array<Point>& quadPoints,
                            Array<double>& quadWeights) const ;

  };
}


#endif
