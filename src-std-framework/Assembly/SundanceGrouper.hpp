/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_GROUPER_H
#define SUNDANCE_GROUPER_H

#include "SundanceDefs.hpp"
#include "SundanceElementIntegral.hpp"
#include "TSFObjectWithVerbosity.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * Grouper
     */
    class GrouperBase
      : public TSFExtended::ObjectWithVerbosity<GrouperBase>
    {
    public:
      /** */
      GrouperBase() {;}

      /** */
      virtual ~GrouperBase(){;}

      /** */
      virtual void findGroups(const RefCountPtr<SparsityPattern>& sparsity,
                              RefCountPtr<LocalMatrixContainer>& A) const = 0 ;
                              
    };

  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
