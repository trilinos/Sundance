/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSubsetCellFilter.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

SubsetCellFilter::SubsetCellFilter(const CellFilter& superset,
                                   const CellPredicate& predicate)
  : CellFilterBase(), superset_(superset), predicate_(predicate)
{;}


XMLObject SubsetCellFilter::toXML() const 
{
  XMLObject rtn("SubsetCellFilter");
  rtn.addAttribute("id", Teuchos::toString(id()));
  rtn.addChild(predicate_.toXML());
  return rtn;
}

bool SubsetCellFilter::lessThan(const CellFilterStub* other) const
{
  const SubsetCellFilter* S 
    = dynamic_cast<const SubsetCellFilter*>(other);

  TEST_FOR_EXCEPTION(S==0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to SubsetCellFilter::lessThan() should be "
                     "a SubsetCellFilter pointer.");

  cerr << "comparing subset cell filter\n me=" 
       << toXML() << endl
       << "you=" << other->toXML() << endl;

  return OrderedPair<CellFilter, CellPredicate>(superset_, predicate_)
    < OrderedPair<CellFilter, CellPredicate>(S->superset_, S->predicate_);
}

CellSet SubsetCellFilter::internalGetCells(const Mesh& mesh) const
{
  SUNDANCE_OUT(verbosity() > VerbLow,
                   "SubsetCellFilter::internalGetCells()");
  CellSet super = superset_.getCells(mesh);

  int dim = superset_.dimension(mesh);

  CellType cellType = mesh.cellType(dim);

  predicate_.setMesh(mesh, dim);

  ExplicitCellSet* rtn = new ExplicitCellSet(mesh, dim, cellType);

  Set<int>& cells = rtn->cells();

  const CellPredicateBase* pred = predicate_.ptr().get();

  
  for (CellIterator i=super.begin(); i != super.end(); i++)
    {
      int LID = *i;
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "SubsetCellFilter is testing " << LID);
      if (pred->test(LID)) 
        {
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       "accepted " << LID);
          cells.insert(LID);
        }
      else
        {
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       "rejected " << LID);
        }
    }

  return rtn;
}
