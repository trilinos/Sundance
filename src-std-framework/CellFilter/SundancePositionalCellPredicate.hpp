/* @HEADER@ */
/* @HEADER@ */


#ifndef SUNDANCE_POSITIONALCELLPREDICATE_H
#define SUNDANCE_POSITIONALCELLPREDICATE_H


#include "SundanceDefs.hpp"
#include "SundanceCellPredicateBase.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;

  /** Prototype for predicates used to filter points. */
  typedef bool positionalPredicate(const Point& x);
  
  namespace Internal
  {
    /** 
     * PositionalCellPredicate tests whether the cell's nodes satisfy
     * a condition on their positions.
     */
    class PositionalCellPredicate : public CellPredicateBase 
    {
    public:
      /** Construct with a function of positions */
      PositionalCellPredicate(positionalPredicate* func) 
        : CellPredicateBase(), func_(func) {;}

      /** virtual dtor */
      virtual ~PositionalCellPredicate(){;}
      
      /** Test whether the cell with the given LID satisfies the condition */
      virtual bool test(int cellLID) const ;

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** comparison */
      virtual bool lessThan(const CellPredicateBase* other) const ;

      /* */
      GET_RCP(CellPredicateBase);

    private:
      positionalPredicate* func_;
    };
  }
}


#endif
