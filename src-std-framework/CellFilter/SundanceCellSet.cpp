/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellSet.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceImplicitCellSet.hpp"
#include "SundanceOut.hpp"
#include "SundanceExceptions.hpp"
#include <algorithm>
#include <iterator>

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CellSet CellSet::setUnion(const CellSet& other) const
{
  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("union", other);
  
  Set<int>& cells = rtn->cells();

  std::set_union(begin(), end(), other.begin(), other.end(), 
                 std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}

CellSet CellSet::setIntersection(const CellSet& other) const
{
  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("intersection", other);
  
  Set<int>& cells = rtn->cells();

  std::set_intersection(begin(), end(), other.begin(), other.end(), 
                        std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}

CellSet CellSet::setDifference(const CellSet& other) const
{
  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("difference", other);
  
  Set<int>& cells = rtn->cells();

  std::set_difference(begin(), end(), other.begin(), other.end(), 
                      std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}


void CellSet::checkCompatibility(const string& op, const CellSet& other) const 
{
  TEST_FOR_EXCEPTION(meshID() != other.meshID(), RuntimeError,
                     "CellSet::checkCompatibility(): "
                     "incompatible mesh ID numbers in " << op
                     << ". LHS=" << meshID() << " RHS=" << other.meshID());

  TEST_FOR_EXCEPTION(dimension() != other.dimension(), RuntimeError,
                     "CellSet::checkCompatibility() incompatible dimensions in " << op
                     << "LHS has "
                     "dimension=" << dimension() << " but RHS has dimension="
                     << other.dimension());
  
  TEST_FOR_EXCEPTION(cellType() != other.cellType(), RuntimeError,
                     "CellSet::checkCompatibility() incompatible cell types. "
                     " in " << op << " LHS has "
                     "cellType=" << cellType() << " but RHS has cellType="
                     << other.cellType());

  SUNDANCE_OUT(verbosity() > VerbMedium,
               "Set operation: " << op << endl
               << "LHS cells: " << *this << endl
               << "RHS cells: " << other);
               
}



