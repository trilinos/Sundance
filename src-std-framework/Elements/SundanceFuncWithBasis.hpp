/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FUNCWITHBASIS_H
#define SUNDANCE_FUNCWITHBASIS_H

#include "SundanceDefs.hpp"
#include "SundanceBasisFamily.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  namespace Internal
  {
    /** 
     * 
     */
    class FuncWithBasis
    {
    public:
      /** */
      FuncWithBasis(const BasisFamily& basis) ;

      /** */
      const Array<BasisFamily>& basis() const {return basis_;}
    private:
      /** */
      Array<BasisFamily> basis_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
