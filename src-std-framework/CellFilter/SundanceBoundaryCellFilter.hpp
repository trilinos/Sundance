/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BOUNDARYCELLFILTER_H
#define SUNDANCE_BOUNDARYCELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceSideCellFilter.hpp"
#include "SundanceBoundaryCellPredicate.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace TSFExtended;
  namespace Internal 
  {
    using namespace Teuchos;

    /** */
    class BoundaryCellFilter : public SubsetCellFilter
    {
    public:
      /** */
      BoundaryCellFilter() : SubsetCellFilter(new SideCellFilter(),
                                              new BoundaryCellPredicate()){;}

      /** */
      virtual ~BoundaryCellFilter(){;}

      /* */
      GET_RCP(CellFilterStub);

    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
