/* @HEADER@ */
/* @HEADER@ */

#include "SundanceMesh.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance::StdMesh::Internal;
using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;


MeshCreationInterface* Mesh::creatableMesh()
{
  MeshCreationInterface* mci 
    = dynamic_cast<MeshCreationInterface*>(ptr().get());
  TEST_FOR_EXCEPTION(mci==0, RuntimeError, 
                     "Mesh::creatableMesh() could not convert mesh to "
                     "a type deriving from MeshCreationInterface");

  return mci;
}



