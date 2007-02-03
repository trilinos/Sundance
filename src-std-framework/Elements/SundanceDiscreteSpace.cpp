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
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init(maximalRegions(1), tuple(basis));
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                             const VectorType<double>& vecType)
  : map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init(maximalRegions(basis.size()), basis);
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                             const Array<CellFilter>& funcDomains,
                             const VectorType<double>& vecType)
  : map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init(funcDomains, basis);
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
                             const CellFilter& funcDomains,
                             const VectorType<double>& vecType)
  : map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init(tuple(funcDomains), tuple(basis));
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                             const CellFilter& funcDomains,
                             const VectorType<double>& vecType)
  : map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init(Array<CellFilter>(basis.size(), funcDomains), basis);
}



DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                             const RefCountPtr<DOFMapBase>& map,
                             const VectorType<double>& vecType)
  : map_(map), 
    mesh_(mesh),
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init(map->funcDomains(), basis);
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
                             const SpectralBasis& spBasis,
                             const VectorType<double>& vecType)
  : map_(),
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init(maximalRegions(spBasis.nterms()), 
       replicate(basis, spBasis.nterms()));
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                             const SpectralBasis& spBasis,
                             const VectorType<double>& vecType)
  : map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init(maximalRegions(basis.size() * spBasis.nterms()), 
       replicate(basis, spBasis.nterms()));
}



void DiscreteSpace::init(const Array<CellFilter>& regions,
                         const BasisArray& basis)
{
  basis_ = basis;
  subdomains_ = regions;
  if (map_.get()==0) 
    {
      Array<Set<CellFilter> > cf(regions.size());
      for (unsigned int i=0; i<regions.size(); i++) cf[i] = makeSet(regions[i]);
      map_ = DOFMapBuilder::makeMap(mesh_, basis, cf);
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

Array<CellFilter> DiscreteSpace::maximalRegions(int n) const
{
  CellFilter cf = new MaximalCellFilter();
  Array<CellFilter> rtn(n, cf);
  return rtn;
}



void DiscreteSpace::importGhosts(const Vector<double>& x,
                                 RefCountPtr<GhostView<double> >& ghostView) const
{
  ghostImporter_->importView(x, ghostView);
}
