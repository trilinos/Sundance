/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SPATIALLYCONSTANTEXPR_H
#define SUNDANCE_SPATIALLYCONSTANTEXPR_H

#include "SundanceLeafExpr.hpp"
#include "SundanceConstantEvaluator.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      /** */
      class SpatiallyConstantExpr : public virtual LeafExpr,
                                    public GenericEvaluatorFactory<SpatiallyConstantExpr, ConstantEvaluator>
        {
        public:
          /** */
          SpatiallyConstantExpr(const double& value) ;

          /** */
          virtual ~SpatiallyConstantExpr(){;}

          /** */
          const double& value() const {return value_;}

          /** */
          void setValue(const double& value) {value_ = value;}

          /** */
          virtual bool isConstant() const {return true;}

      /** 
       * Determine which functional and spatial derivatives are nonzero in the
       * given context. We also keep track of which derivatives
       * are known to be constant, which can simplify evaluation. 
       */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                bool regardFuncsAsConstant) const ;

        private:
          double value_;
        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
