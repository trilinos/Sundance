/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_LISTEXPR_H
#define SUNDANCE_LISTEXPR_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "Teuchos_Array.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using std::string;

  namespace Internal
    {
      /** */
      class ListExpr : public ExprBase
        {
        public:
          /** */
          ListExpr();

          /** */
          ListExpr(const Array<Expr>& elements);

          /** */
          virtual ~ListExpr() {;}

          /** */
          const Expr& element(int i) const {return elements_[i];}

          /** */
          void append(const Expr& expr);

          /** */
          Expr flatten() const ;

          /** */
          Expr join(const Expr& other) const ;

          /** */
          int size() const ;

          /** */
          int totalSize() const ;

          /** Write a simple text description suitable 
           * for output to a terminal */
          virtual ostream& toText(ostream& os, bool paren) const ;

          /** Write in a form suitable for LaTeX formatting */
          virtual ostream& toLatex(ostream& os, bool paren) const ;


          /** Write in XML */
          virtual XMLObject toXML() const ;

          /** */
          virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

        private:
          Array<Expr> elements_;
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
