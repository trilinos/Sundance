/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_USERDEFOP_H
#define SUNDANCE_USERDEFOP_H

#include "SundanceDefs.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"



namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

    /**
     *
     */
    class UserDefOp : public ExprWithChildren
    {
    public:
      /** construct with an argument and the functor defining the operation */
      UserDefOp(const Expr& arg,
                const RefCountPtr<UserDefFunctor>& op);

      /** virtual destructor */
      virtual ~UserDefOp() {;}

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
      const UserDefFunctor* op() const {return op_.get();}
    private:
      RefCountPtr<UserDefFunctor> op_;

      /** */
      static Array<RefCountPtr<ScalarExpr> > getScalarArgs(const Expr& args);

    };
}

#endif
