/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellFilterBase.hpp"
#include "SundanceOut.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CellFilterBase::CellFilterBase()
  : CellFilterStub(), cellSetCache_()
{;}

CellSet CellFilterBase::getCells(const Mesh& mesh) const
{
  if (cellSetCache_.ptr().get()==0)
    {
      cellSetCache_ = internalGetCells(mesh);
    }
  return cellSetCache_;
//   print(cerr);

  // if (!cellSetCache_.containsKey(mesh.id()))
//     {
//       SUNDANCE_OUT(verbosity() > VerbMedium,
//                    "cell set " << toXML() << " is computing cell cache");

//       CellSet cells = internalGetCells(mesh);
//       cellSetCache_.put(mesh.id(), cells);
//     }
//   else
//     {
//       SUNDANCE_OUT(verbosity() > VerbMedium,
//                    "cell set " << toXML() << " is reusing cell cache");
//     }
//  return cellSetCache_.get(mesh.id());

}



