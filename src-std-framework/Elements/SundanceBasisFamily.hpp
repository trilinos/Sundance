/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BASISFAMILY_H
#define SUNDANCE_BASISFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceBasisFamilyBase.hpp"
#include "SundanceOrderedHandle.hpp"
#include "Teuchos_Array.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceStdFwk::Internal;
  using namespace TSFExtended;

  /** 
   * BasisFamily is the user-level handle class for specifying the
   * basis with which a test, unknown, or discrete function is represented.
   * Basis functions can be vector-valued, as is the case with, for example,
   * the Nedelec basis in electromagnetics; the dim() method
   * returns the spatial dimension of the basis functions. Scalar-valued
   * bases naturally have dim()=1. 
   */
  class BasisFamily : public OrderedHandle<BasisFamilyBase>
    {
    public:
      /* handle ctor boilerplate */
      ORDERED_HANDLE_CTORS(BasisFamily, BasisFamilyBase);

      /** write to XML */
      XMLObject toXML() const ;

      /** return the polynomial order to which the basis is complete */
      int order() const ;

      /** return the spatial dimension of the basis elements */
      int dim() const ;

      /** return the number of nodes for this basis on the given cell type */
      int nNodes(const CellType& cellType) const ;

      /** */
      bool operator==(const BasisFamily& other) const ;

      /** Sum up the dim() values for array of bases. */
      static int size(const Array<BasisFamily>& b) ;

      /** Extract the basis from an expression */
      static BasisFamily getBasis(const Expr& expr);
    };
}

#endif
