/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDOFMapBase.hpp"
#include "Teuchos_MPIContainerComm.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;



DOFMapBase::DOFMapBase(const Mesh& mesh)
  : localProcID_(mesh.comm().getRank()),
    mesh_(mesh),
    cellSets_(),
    funcIDOnCellSets_(),
    cellDimOnCellSets_(),
    lowestLocalDOF_(),
    numDOFs_(),
    dofsHaveBeenAssigned_()
{;}

