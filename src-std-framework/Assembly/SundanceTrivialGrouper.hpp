/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TRIVIALGROUPER_H
#define SUNDANCE_TRIVIALGROUPER_H

#include "SundanceDefs.hpp"
#include "SundanceGrouperBase.hpp"

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
    class TrivialGrouper : public GrouperBase
    {
    public:
      /** */
      TrivialGrouper() : GrouperBase() {;}

      /** */
      virtual ~TrivialGrouper(){;}

      /** */
      virtual void findGroups(const EquationSet& eqn,
                              const CellType& cellType,
                              int cellDim,
                              const QuadratureFamily& quad,
                              const RefCountPtr<SparsityPattern>& sparsity,
                              Array<IntegralGroup>& groups) const ;
                              
    };

  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
