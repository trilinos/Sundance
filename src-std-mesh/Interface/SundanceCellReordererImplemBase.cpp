/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellReordererBase.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance::StdMesh::Internal;
using namespace Sundance;
using namespace Teuchos;


CellReordererImplemBase::CellReordererImplemBase(const MeshBase* mesh)
  : mesh_(mesh) 
{;}



int CellReordererImplemBase::end() const
{
  return mesh_->numCells(mesh_->spatialDim());
}
