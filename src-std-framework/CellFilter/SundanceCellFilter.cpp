/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellFilter.hpp"
#include "SundanceBinaryCellFilter.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceLabelCellPredicate.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;


CellSet CellFilter::getCells(const Mesh& mesh) const
{
  return ptr()->getCells(mesh);
}



int CellFilter::dimension(const Mesh& mesh) const
{
  return ptr()->dimension(mesh);
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



bool CellFilter::operator<(const CellFilter& other) const 
{
  string me = ptr()->typeName();
  string you = other.ptr()->typeName();
  if (me < you) return true;
  if (you < me) return false;
  return ptr()->lessThan(other.ptr().get());
}
