/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SPATIALLYCONSTANTEXPR_H
#define SUNDANCE_SPATIALLYCONSTANTEXPR_H

#include "SundanceLeafExpr.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      /** */
      class SpatiallyConstantExpr : public virtual LeafExpr
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



        private:
          double value_;
        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
