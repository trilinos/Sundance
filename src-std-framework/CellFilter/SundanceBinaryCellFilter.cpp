/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBinaryCellFilter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

BinaryCellFilter::BinaryCellFilter(const CellFilter& left,
                                   const CellFilter& right,
                                   const CellFilterOpType& op)
  : CellFilterBase(), left_(left), right_(right), op_(op)
{;}

int BinaryCellFilter::dimension(const Mesh& mesh) const
{
  int d1 = left_.dimension(mesh);
  int d2 = right_.dimension(mesh);

  TEST_FOR_EXCEPTION(d1 != d2, RuntimeError,
                     "BinaryCellFilter::dimension() mismatched dimensions. "
                     "Left filter has dimension d1=" << d1 << " but "
                     "right filter has dimension d2=" << d2);

  return d1;
}

CellSet BinaryCellFilter::internalGetCells(const Mesh& mesh) const
{
  Tabs tab;

  SUNDANCE_OUT(verbosity() > VerbMedium,
               "cell filter " << toXML().toString() << " is getting its cells for mesh " 
               << mesh.id());

  CellSet L = left_.getCells(mesh);
  CellSet R = right_.getCells(mesh);


  SUNDANCE_OUT(verbosity() > VerbMedium,
               "cell filter " << toXML().toString() << " is performing its operation");
  
  switch(op_)
    {
    case Union:
      return L.setUnion(R);
    case Intersection:
      return L.setIntersection(R);
    case Difference:
      return L.setDifference(R);
    }
}

string BinaryCellFilter::opName() const 
{
  switch(op_)
    {
    case Union:
      return "UnionCellFilter";
    case Intersection:
      return "IntersectionCellFilter";
    case Difference:
      return "DifferenceCellFilter";
    }
}

XMLObject BinaryCellFilter::toXML() const 
{
  XMLObject rtn(opName());
  rtn.addChild(left_.toXML());
  rtn.addChild(right_.toXML());
  return rtn;
}

bool BinaryCellFilter::lessThan(const CellFilterStub* other) const
{
  const BinaryCellFilter* B 
    = dynamic_cast<const BinaryCellFilter*>(other);

  TEST_FOR_EXCEPTION(B==0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to BinaryCellFilter::lessThan() should be "
                     "a BinaryCellFilter pointer.");

  return OrderedTriple<CellFilterOpType, CellFilter, CellFilter>(op_, left_, right_) 
    < OrderedTriple<CellFilterOpType, CellFilter, CellFilter>(B->op_, B->left_, B->right_) ;
}
