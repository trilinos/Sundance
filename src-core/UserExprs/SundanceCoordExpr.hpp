/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_COORDEXPR_H
#define SUNDANCE_COORDEXPR_H

#include "SundanceFuncElementBase.hpp"
#include "SundanceLeafExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceCoordExprEvaluator.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace SundanceCore::Internal;

  /** 
   * CoordExpr is an 
   * expression that returns one of the cartesian coordinates for
   * each point at which it evaluated. Which coordinate (i.e., \f$x\f$,
   * \f$y\f$, or \f$z\f$) to be represented is specified by the first
   * argement to the constructor. 
   */
  class CoordExpr : public FuncElementBase,
                    public GenericEvaluatorFactory<CoordExpr, CoordExprEvaluator>,
                    virtual public LeafExpr
    {
    public:
      /** */
      CoordExpr(int dir, const string& name="");

      /** */
      virtual ~CoordExpr() {;}

      /** */
      virtual XMLObject toXML() const ;

#ifndef DOXYGEN_DEVELOPER_ONLY
      /** */
      int dir() const {return dir_;}


    

      /** 
       * Determine which functional and spatial derivatives are nonzero in the
       * given context. We also keep track of which functional derivatives
       * are known to be constant, which can simplify evaluation. 
       */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                const Set<int>& allFuncIDs,
                                bool regardFuncsAsConstant) const ;

      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

      static string coordName(int dir, const string& name);

    private:
      int dir_;
#endif  /* DOXYGEN_DEVELOPER_ONLY */

    };
}

#endif
