/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BOUNDARYCELLFILTER_H
#define SUNDANCE_BOUNDARYCELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceSideCellFilter.hpp"
#include "SundanceBoundaryCellPredicate.hpp"


namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace TSFExtended;
  using namespace Teuchos;

  /** 
   * BoundaryCellFilter identifies all cells of dimension \f$D-1\f$
   on the boundary. The boundary cells can be identified topologically
   as those cells of \f$D-1\f$ having only one cofacet.
  */
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


#endif
