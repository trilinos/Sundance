/* @HEADER@ */
/* @HEADER@ */



#ifndef SUNDANCE_PRODUCTEXPR_H
#define SUNDANCE_PRODUCTEXPR_H

#include "SundanceBinaryExpr.hpp"
#include "SundanceNonlinearExpr.hpp"
#include "SundanceProductEvaluator.hpp"

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
     * ProductExpr represents a product of two scalar-valued expressions
     */
    class ProductExpr : public BinaryExpr,
                        public NonlinearExpr,
                        public GenericEvaluatorFactory<ProductExpr, ProductEvaluator>
    {
    public:
      /** */
      ProductExpr(const RefCountPtr<ScalarExpr>& left,
                  const RefCountPtr<ScalarExpr>& right);

      /** virtual dtor */
      virtual ~ProductExpr() {;}

      /** Indicate whether this expression is a "hungry"
       * differential operator that is awaiting an argument. */
      virtual bool isHungryDiffOp() const ;

      /** Preprocessing step to determine which functional 
       * derivatives are nonzero */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                bool regardFuncsAsConstant) const ;

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

    private:

      void findChildActiveFuncs(const Set<MultiSet<int> >& funcIDs,
                                Set<MultiSet<int> >& leftFuncs,
                                Set<MultiSet<int> >& rightFuncs) const ;

      /** */
      void findChildMultiIndexSets(const Set<MultiIndex>& miSet,
                                   Set<MultiIndex>& miLeft,
                                   Set<MultiIndex>& miRight) const ;

    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
