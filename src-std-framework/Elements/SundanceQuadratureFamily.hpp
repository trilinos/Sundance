/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_QUADRATUREFAMILY_H
#define SUNDANCE_QUADRATUREFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdFwk
{
  using namespace TSFExtended;
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** 
   * QuadratureFamily is a geometry-independent specification of
   * a method by which quadrature is to be carried out. For example,
   * a GaussianQuadrature family will generate Gaussian
   * quadrature points on any cell type.
   */
  class QuadratureFamily : public Handle<QuadratureFamilyStub>
  {
  public:
    /* */
    HANDLE_CTORS(QuadratureFamily, QuadratureFamilyStub);
    /** */
    XMLObject toXML() const ;

    /** */
    int order() const ;

    /** Get the quadrature points and weights for the given cell type */
    void getPoints(const CellType& cellType, 
                   Array<Point>& quadPoints,
                   Array<double>& quadWeights) const ;
  private:
  };
}

#endif
