/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SCALAREXPR_H
#define SUNDANCE_SCALAREXPR_H


#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include "SundanceExprBase.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
    {
      /** */
      class ScalarExpr : public ExprBase
        {
        public:
          /** empty ctor */
          ScalarExpr();

          /** virtual destructor */
          virtual ~ScalarExpr() {;}


          /** Indicate whether this expression is constant in space */
          virtual bool isConstant() const {return false;}

          /** Indicate whether this expression is a "hungry"
           * differential operator that is awaiting an argument. */
          virtual bool isHungryDiffOp() const {return false;}

        protected:
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
