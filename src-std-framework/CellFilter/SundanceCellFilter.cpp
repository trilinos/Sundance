/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellFilter.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceBinaryCellFilter.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceNullCellFilterStub.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

bool CellFilter::isNullCellFilter() const 
{
  return dynamic_cast<const NullCellFilterStub*>(ptr().get()) != 0;
}

CellSet CellFilter::getCells(const Mesh& mesh) const
{
  if (isNullCellFilter())
    {
      return new ExplicitCellSet(mesh, -1, 
                                 NullCell);
    }
  return cfbPtr()->getCells(mesh);
}



int CellFilter::dimension(const Mesh& mesh) const
{
  if (isNullCellFilter())
    {
      return -1;
    }
  return cfbPtr()->dimension(mesh);
}



CellFilter CellFilter::operator+(const CellFilter& other) const 
{
  return new BinaryCellFilter(*this, other, BinaryCellFilter::Union);
}



CellFilter CellFilter::operator-(const CellFilter& other) const 
{
  return new BinaryCellFilter(*this, other, BinaryCellFilter::Difference);
}



CellFilter CellFilter::intersection(const CellFilter& other) const 
{
  return new BinaryCellFilter(*this, other, BinaryCellFilter::Intersection);
}



CellFilter CellFilter::labeledSubset(int label) const
{
  CellPredicate pred = new LabelCellPredicate(label);
  return new SubsetCellFilter(*this, pred);
}


CellFilter CellFilter::subset(const CellPredicate& pred) const
{
  return new SubsetCellFilter(*this, pred);
}



XMLObject CellFilter::toXML() const 
{
  return ptr()->toXML();
}

const CellFilterBase* CellFilter::cfbPtr() const
{
  const CellFilterBase* rtn = dynamic_cast<CellFilterBase*>(ptr().get());
  TEST_FOR_EXCEPTION(rtn==0, InternalError, "CellFilter::cfbPtr() cast failed");
  return rtn;
}

