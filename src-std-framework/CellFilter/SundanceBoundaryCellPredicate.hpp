/* @HEADER@ */
/* @HEADER@ */


#ifndef SUNDANCE_BOUNDARYCELLPREDICATE_H
#define SUNDANCE_BOUNDARYCELLPREDICATE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellPredicateBase.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;

  
  namespace Internal
  {
    /** 
     * BoundaryCellPredicate tests whether a cell is on the 
     * boundary of the domain
     */
    class BoundaryCellPredicate : public CellPredicateBase 
    {
    public:
      /** Empty ctor */
      BoundaryCellPredicate() {;}

      /** virtual dtor */
      virtual ~BoundaryCellPredicate(){;}
      
      /** Test whether the cell with the given LID satisfies the condition */
      virtual bool test(int cellLID) const ;

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** comparison */
      virtual bool lessThan(const CellPredicateBase* other) const ;

      /* */
      GET_RCP(CellPredicateBase);

    private:
    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
