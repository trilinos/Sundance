/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceDiscreteSpace.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceMaximalCellFilter.hpp"

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
    regions_(maximalRegions(basis_.size())),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init();
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                             const VectorType<double>& vecType)
  : map_(), 
    mesh_(mesh), 
    basis_(basis), 
    regions_(maximalRegions(basis_.size())),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init();
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                             const Array<Set<CellFilter> >& regions,
                             const VectorType<double>& vecType)
  : map_(), 
    mesh_(mesh), 
    basis_(basis), 
    regions_(regions),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init();
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                             const RefCountPtr<DOFMapBase>& map,
                             const VectorType<double>& vecType)
  : map_(map), 
    mesh_(mesh), 
    basis_(basis), 
    regions_(maximalRegions(basis_.size())),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init();
}



void DiscreteSpace::init()
{
  if (map_.get()==0) 
    {
      map_ = DOFMapBuilder::makeMap(mesh_, basis_, regions_);
    }
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

Array<Set<CellFilter> > DiscreteSpace::maximalRegions(int n) const
{
  CellFilter cf = new MaximalCellFilter();
  Array<Set<CellFilter> > rtn(n);
  for (int i=0; i<n; i++) 
    {
      Set<CellFilter> s;
      s.put(cf);
      rtn[i] = s;
    }
  return rtn;
}



void DiscreteSpace::importGhosts(const Vector<double>& x,
                                 RefCountPtr<GhostView<double> >& ghostView) const
{
  ghostImporter_->importView(x, ghostView);
}
