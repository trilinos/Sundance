/* @HEADER@ */
/* @HEADER@ */

#include "SundanceMesh.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdMesh::Internal;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace SundanceUtils;
using namespace TSFExtended;


CreatableMesh* Mesh::creatableMesh()
{
  CreatableMesh* mci 
    = dynamic_cast<CreatableMesh*>(ptr().get());
  TEST_FOR_EXCEPTION(mci==0, RuntimeError, 
                     "Mesh::creatableMesh() could not convert mesh to "
                     "a type deriving from CreatableMesh");

  return mci;
}



