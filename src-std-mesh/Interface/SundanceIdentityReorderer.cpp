/* @HEADER@ */
/* @HEADER@ */

#include "SundanceIdentityReorderer.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdMesh::Internal;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace SundanceUtils;

IdentityReordererImplem::IdentityReordererImplem(const MeshBase* mesh) 
        : CellReordererImplemBase(mesh) {;}
