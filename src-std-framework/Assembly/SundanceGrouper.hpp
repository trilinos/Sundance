/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_GROUPER_H
#define SUNDANCE_GROUPER_H

#include "SundanceDefs.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceIntegralGroup.hpp"


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
      virtual void findGroups(const RefCountPtr<SparsitySuperset>& sparsity,
                              Array<IntegralGroup>& groups) const = 0 ;
                              
    };

  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
