/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FUNCTIONALDERIV_H
#define SUNDANCE_FUNCTIONALDERIV_H

#include "SundanceDefs.hpp"
#include "SundanceDerivBase.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMultiIndex.hpp"
#include "SundanceFuncElementBase.hpp"
#include "Teuchos_RefCountPtr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY



namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      /**
       * FunctionalDeriv represents a single first-order derivative
       * with respect to the derivative of a function. It is used internally
       * to represent
       * functional differentiation.
       *
       * @see Deriv
       */
      class FunctionalDeriv : public DerivBase
        {
        public:
          /** */
          FunctionalDeriv(const FuncElementBase* func,
                          const MultiIndex& mi);

          /** */
          virtual ~FunctionalDeriv(){;}

          /** */
          virtual bool lessThan(const Deriv& other) const ;

          /** */
          virtual string toString() const ;

          /** */
          virtual int funcID() const {return func_->funcID();}

          /** */
          const MultiIndex& multiIndex() const {return mi_;}

          /** */
          Deriv derivWrtMultiIndex(const MultiIndex& mi) const ;

          /** */
          const FuncElementBase* func() const {return func_;}

          /** */
          virtual RefCountPtr<DerivBase> getRcp()
            {return rcp(this);}
        private:
          const FuncElementBase* func_;

          MultiIndex mi_;
        };
    }
}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
