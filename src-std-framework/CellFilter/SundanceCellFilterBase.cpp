/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellFilterBase.hpp"
#include "SundanceOut.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CellFilterBase::CellFilterBase()
  : cellSetCache_()
{;}

CellSet CellFilterBase::getCells(const Mesh& mesh) const
{
  if (!cellSetCache_.containsKey(mesh.id()))
    {
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "cell set " << toXML() << " is computing cell cache");
      cellSetCache_.put(mesh.id(), internalGetCells(mesh));
    }
  else
    {
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "cell set " << toXML() << " is reusing cell cache");
    }
  return cellSetCache_.get(mesh.id());
}



