/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellReordererBase.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdMesh::Internal;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace SundanceUtils;


CellReordererImplemBase::CellReordererImplemBase(const MeshBase* mesh)
  : mesh_(mesh) 
{;}



int CellReordererImplemBase::end() const
{
  return mesh_->numCells(mesh_->spatialDim());
}
