/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BASISFAMILY_H
#define SUNDANCE_BASISFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceBasisFamilyBase.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  using namespace SundanceStdFwk::Internal;

  /** */
  class BasisFamily 
    {
    public:
      /** */
      BasisFamily();
      /** */
      BasisFamily(BasisFamilyBase* ptr);

      /** */
      XMLObject toXML() const ;

      /** */
      int order() const ;

#ifndef DOXYGEN_DEVELOPER_ONLY
      /** */
      RefCountPtr<BasisFamilyBase> ptr() const {return ptr_;}

    private:

      RefCountPtr<BasisFamilyBase> ptr_;

#endif  /* DOXYGEN_DEVELOPER_ONLY */

    };
}

#endif
