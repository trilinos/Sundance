/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SYMBOLICTRANSFORMATION_H
#define SUNDANCE_SYMBOLICTRANSFORMATION_H

#include "SundanceDefs.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "Teuchos_RefCountPtr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY 

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace TSFExtended;
  using std::string;
  


  namespace Internal
    {
      class ScalarExpr;

      /** 
       * SumTransformation is a base class for any transformation
       * which takes the two operands of a sum (left, right) and produces
       * a new expression mathematically equivalent to the original
       * sum. This will be used to effect simplification
       * transformations on sum expressions.
       */
      class SymbolicTransformation : public ObjectWithVerbosity<SymbolicTransformation>
        {
        public:
          /** */
          SymbolicTransformation();

          /** */
          virtual ~SymbolicTransformation(){;}

          /** Returns -expr if sign == -1, otherwise returns expr */
          static RefCountPtr<ScalarExpr> chooseSign(int sign, 
                                                    const RefCountPtr<ScalarExpr>& expr);

          /** Returns -expr if sign == -1, otherwise returns expr */
          static Expr chooseSign(int sign, 
                                 const Expr& expr);

          /** extract the underlying ScalarExpr from an Expr. */
          static RefCountPtr<ScalarExpr> getScalar(const Expr& expr); 
        };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
