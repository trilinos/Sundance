/* @HEADER@ */
/* @HEADER@ */



#ifndef SUNDANCE_PRODUCTEXPR_H
#define SUNDANCE_PRODUCTEXPR_H

#include "SundanceBinaryExpr.hpp"


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
    class ProductExpr : public BinaryExpr
    {
    public:
      /** */
      ProductExpr(const RefCountPtr<ScalarExpr>& left,
                  const RefCountPtr<ScalarExpr>& right);

      /** */
      virtual ~ProductExpr() {;}

      /** Indicate whether this expression is a "hungry"
       * differential operator that is awaiting an argument. */
      virtual bool isHungryDiffOp() const ;

      /**
       * Indicate whether the given functional derivative is nonzero
       */
      virtual bool hasNonzeroDeriv(const MultipleDeriv& d) const ;

      /** Return the set of derivatives required by the operands
       * of this expression given that this expression
       * requires the set d. */
      virtual Array<DerivSet>
      derivsRequiredFromOperands(const DerivSet& d) const ;


      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}
    protected:
      /** */
      virtual bool parenthesizeSelf() const {return true;}
      /** */
      virtual bool parenthesizeOperands() const {return true;}
      /** */
      virtual const string& xmlTag() const ;
      /** */
      virtual const string& opChar() const ;

    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
