/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_NONLINEARUNARYOP_H
#define SUNDANCE_NONLINEARUNARYOP_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"



#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
  {
    /**
     *
     */
    class NonlinearUnaryOp : public UnaryExpr
    {
    public:
      /** construct with an argument and the functor defining the operation */
      NonlinearUnaryOp(const RefCountPtr<ScalarExpr>& arg, 
                       const RefCountPtr<UnaryFunctor>& op);

      /** virtual destructor */
      virtual ~NonlinearUnaryOp() {;}

      /** Write a simple text description suitable
       * for output to a terminal */
      virtual ostream& toText(ostream& os, bool paren) const ;

      /** Write in a form suitable for LaTeX formatting */
      virtual ostream& toLatex(ostream& os, bool paren) const ;

      /** Write in XML */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

      /** The function is nonlinear, so 
       * all mixed partials are assumed to be nonzero as long as
       * their constituent single derivs are nonzero. */
      virtual bool hasNonzeroDeriv(const MultipleDeriv& d) const ;


      /** Access to the operator */
      const UnaryFunctor* op() const {return op_.get();}
    private:
      RefCountPtr<UnaryFunctor> op_;

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
