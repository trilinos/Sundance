/* @HEADER@ */
/* @HEADER@ */

#include "SundanceIdentityReorderer.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance::StdMesh::Internal;
using namespace Sundance;
using namespace Teuchos;

IdentityReordererImplem::IdentityReordererImplem(const MeshBase* mesh) 
        : CellReordererImplemBase(mesh) {;}
