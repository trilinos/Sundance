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
      FuncWithBasis(const BasisFamily& basis) : basis_(tuple(basis)) {;}

      /** */
      FuncWithBasis(const Array<BasisFamily>& basis) :basis_(basis) {;}

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
