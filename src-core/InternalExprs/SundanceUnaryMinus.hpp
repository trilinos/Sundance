/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_UNARYMINUS_H
#define SUNDANCE_UNARYMINUS_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceUnaryMinusEvaluator.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using std::string;

  namespace Internal
    {
      /** */
      class UnaryMinus : public UnaryExpr,
                         GenericEvaluatorFactory<UnaryMinus, UnaryMinusEvaluator>
        {
        public:
          /** construct with the argument */
          UnaryMinus(const RefCountPtr<ScalarExpr>& arg);

          /** virtual dtor */
          virtual ~UnaryMinus() {;}


          /** */
          virtual ostream& toText(ostream& os, bool paren) const ;

          /** */
          virtual ostream& toLatex(ostream& os, bool paren) const ;

          /** */
          virtual XMLObject toXML() const ;

          /** */
          virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

        };
    }
}
                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  



#endif
