/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CONSTANTEXPR_H
#define SUNDANCE_CONSTANTEXPR_H

#include "SundanceSpatiallyConstantExpr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      class ConstantExpr : public SpatiallyConstantExpr
        {
        public:
          ConstantExpr(const double& value);
          virtual ~ConstantExpr() {;}

          /** */
          virtual ostream& toText(ostream& os, bool paren) const ;

          /** */
          virtual ostream& toLatex(ostream& os, bool paren) const ;

          /** */
          virtual XMLObject toXML() const ;



          /**
           * Indicate whether the given functional derivative is nonzero.
           * A constant expression has a nonzero derivative only if the
           * order of the derivative is zero. 
           */
          virtual bool hasNonzeroDeriv(const MultipleDeriv& d) const
          {return d.order()==0;}

          /**
           * Find all functions and their derivatives beneath my level
           * in the tree. A constant expr has no functions beneath it,
           * so this method does nothing.
           */
          virtual void getRoughDependencies(Set<Deriv>& /* funcs */) const {;}

          /** */
          virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}
        protected:
        private:
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
