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
  : map_(),
    mesh_(mesh), 
    basis_(tuple(basis)), 
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  map_ = rcp(new HomogeneousDOFMap(mesh, basis_[0], basis_.size()));

  int nDof = map_->numLocalDOFs();
  int lowDof = map_->lowestLocalDOF();

  Array<int> dofs(nDof);
  for (int i=0; i<nDof; i++) dofs[i] = lowDof + i;

  vecSpace_ = vecType_.createSpace(map_->numDOFs(),
                                   map_->numLocalDOFs(),
                                   &(dofs[0]));

  RefCountPtr<Array<int> > ghostIndices = map_->ghostIndices();
  int nGhost = ghostIndices->size();
  int* ghosts = 0;
  if (nGhost!=0) ghosts = &((*ghostIndices)[0]);

  ghostImporter_ = vecType_.createGhostImporter(vecSpace_, nGhost, ghosts);
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const Array<BasisFamily>& basis,
                             const VectorType<double>& vecType)
  : map_(), mesh_(mesh), basis_(basis), vecSpace_(), vecType_(vecType),
    ghostImporter_()
{
  map_ = rcp(new HomogeneousDOFMap(mesh, basis_[0], basis_.size()));

  int nDof = map_->numLocalDOFs();
  int lowDof = map_->lowestLocalDOF();

  Array<int> dofs(nDof);
  for (int i=0; i<nDof; i++) dofs[i] = lowDof + i;

  vecSpace_ = vecType_.createSpace(map_->numDOFs(),
                                   map_->numLocalDOFs(),
                                   &(dofs[0]));

  RefCountPtr<Array<int> > ghostIndices = map_->ghostIndices();
  int nGhost = ghostIndices->size();
  int* ghosts = 0;
  if (nGhost!=0) ghosts = &((*ghostIndices)[0]);

  ghostImporter_ = vecType_.createGhostImporter(vecSpace_, nGhost, ghosts);
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const Array<BasisFamily>& basis,
                             const RefCountPtr<DOFMapBase>& map,
                             const VectorType<double>& vecType)
  : map_(map), mesh_(mesh), basis_(basis), vecSpace_(), vecType_(vecType),
    ghostImporter_()
{
  int nDof = map_->numLocalDOFs();
  int lowDof = map_->lowestLocalDOF();

  Array<int> dofs(nDof);
  for (int i=0; i<nDof; i++) dofs[i] = lowDof + i;

  vecSpace_ = vecType_.createSpace(map_->numDOFs(),
                                   map_->numLocalDOFs(),
                                   &(dofs[0]));

  RefCountPtr<Array<int> > ghostIndices = map_->ghostIndices();
  int nGhost = ghostIndices->size();
  int* ghosts = 0;
  if (nGhost!=0) ghosts = &((*ghostIndices)[0]);

  ghostImporter_ = vecType_.createGhostImporter(vecSpace_, nGhost, ghosts);
}

void DiscreteSpace::importGhosts(const Vector<double>& x,
                                 RefCountPtr<GhostView<double> >& ghostView) const
{
  ghostImporter_->importView(x, ghostView);
}
