/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteSpace.hpp"
#include "SundanceHomogeneousDOFMap.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
                             const VectorType<double>& vecType)
  : map_(), mesh_(mesh), basis_(tuple(basis)), vecSpace_(), vecType_(vecType)
{
  map_ = rcp(new HomogeneousDOFMap(mesh, basis_[0], basis_.size()));

  int nDof = map_->numLocalDOFs();
  int lowDof = map_->lowestLocalDOF();

  Array<int> dofs(nDof);
  for (int i=0; i<nDof; i++) dofs[i] = lowDof + i;

  vecSpace_ = vecType_.createSpace(map_->numDOFs(),
                                   map_->numLocalDOFs(),
                                   &(dofs[0]));
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const Array<BasisFamily>& basis,
                             const VectorType<double>& vecType)
  : map_(), mesh_(mesh), basis_(basis), vecSpace_(), vecType_(vecType)
{
  map_ = rcp(new HomogeneousDOFMap(mesh, basis_[0], basis_.size()));

  int nDof = map_->numLocalDOFs();
  int lowDof = map_->lowestLocalDOF();

  Array<int> dofs(nDof);
  for (int i=0; i<nDof; i++) dofs[i] = lowDof + i;

  vecSpace_ = vecType_.createSpace(map_->numDOFs(),
                                   map_->numLocalDOFs(),
                                   &(dofs[0]));
}
