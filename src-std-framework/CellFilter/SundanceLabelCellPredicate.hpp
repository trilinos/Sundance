/* @HEADER@ */
/* @HEADER@ */


#ifndef SUNDANCE_LABELCELLPREDICATE_H
#define SUNDANCE_LABELCELLPREDICATE_H

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
     * LabelCellPredicate tests whether a cell's label is equal to
     * a reference label. 
     */
    class LabelCellPredicate : public CellPredicateBase 
    {
    public:
      /** Construct with a label string */
      LabelCellPredicate(int label) 
        : CellPredicateBase(), labelIndex_(label){;}

      /** virtual dtor */
      virtual ~LabelCellPredicate(){;}
      
      /** Test whether the cell with the given LID satisfies the condition */
      virtual bool test(int cellLID) const ;

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** comparison */
      virtual bool lessThan(const CellPredicateBase* other) const ;

      /* */
      GET_RCP(CellPredicateBase);

    private:
      string label_;

      mutable int labelIndex_;

    };
  }
}

#endif
