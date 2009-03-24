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
#include "SundanceOut.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "TSFProductVectorSpaceImpl.hpp"
#include "Teuchos_MPIContainerComm.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

const int* vecPtr(const Array<int>& x)
{
  static const int* dum = 0;
  if (x.size()==0) return dum;
  else return &(x[0]);
}

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
  init(maximalRegions(1), BasisArray(tuple(basis)));
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
  init(tuple(funcDomains), BasisArray(tuple(basis)));
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
  const RefCountPtr<Array<int> >& bcIndices,
  const VectorType<double>& vecType)
  : map_(map), 
    mesh_(mesh),
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
{
  init(map->funcDomains(), basis, bcIndices, true);
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



void DiscreteSpace::init(
  const Array<CellFilter>& regions,
  const BasisArray& basis)
{
  RefCountPtr<Array<int> > dummyBCIndices;
  init(regions, basis, dummyBCIndices, false);
}

void DiscreteSpace::init(
  const Array<CellFilter>& regions,
  const BasisArray& basis,
  const RefCountPtr<Array<int> >& isBCIndex, 
  bool partitionBCs)
{
  basis_ = basis;
  subdomains_ = regions;
  if (map_.get()==0) 
    {
      Array<Set<CellFilter> > cf(regions.size());
      for (unsigned int i=0; i<regions.size(); i++) cf[i] = makeSet(regions[i]);
      DOFMapBuilder b;
      map_ = b.makeMap(mesh_, basis, cf);
    }

  int nDof = map_->numLocalDOFs();
  int lowDof = map_->lowestLocalDOF();

  if (partitionBCs)
  {
    TEST_FOR_EXCEPT(isBCIndex.get() == 0);

    int nBCDofs = 0;
    for (int i=0; i<nDof; i++)
    {
      if ((*isBCIndex)[i]) nBCDofs++;
    }
    
    int nTotalBCDofs = nBCDofs;
    mesh().comm().allReduce(&nBCDofs, &nTotalBCDofs, 1, MPIComm::INT, MPIComm::SUM);
    int nTotalInteriorDofs = map_->numDOFs() - nTotalBCDofs;

    Array<int> interiorDofs(nDof - nBCDofs);
    Array<int> bcDofs(nBCDofs);
    int iBC = 0;
    int iIn = 0;
    for (int i=0; i<nDof; i++)
    {
      if ((*isBCIndex)[i]) bcDofs[iBC++] = lowDof+i;
      else interiorDofs[iIn++] = lowDof+i;
    }
    const int* bcDofPtr = vecPtr(bcDofs);
    const int* intDofPtr = vecPtr(interiorDofs);
    Out::os() << "creating BC space" << endl;
    VectorSpace<double> bcSpace = vecType_.createSpace(nTotalBCDofs, nBCDofs,
      bcDofPtr);
    Out::os() << "bc space = " << bcSpace << endl;
    VectorSpace<double> interiorSpace = vecType_.createSpace(nTotalInteriorDofs, nDof-nBCDofs,
      intDofPtr);
    Out::os() << "int space = " << interiorSpace << endl;

    vecSpace_ = productSpace<double>(interiorSpace, bcSpace);
  }
  else
  {
    Array<int> dofs(nDof);
    for (int i=0; i<nDof; i++) dofs[i] = lowDof + i;
    
    vecSpace_ = vecType_.createSpace(map_->numDOFs(),
      map_->numLocalDOFs(),
      &(dofs[0]));
  }

  if (!partitionBCs)
  {
    RefCountPtr<Array<int> > ghostIndices = map_->ghostIndices();
    int nGhost = ghostIndices->size();
    int* ghosts = 0;
    if (nGhost!=0) ghosts = &((*ghostIndices)[0]);
    ghostImporter_ = vecType_.createGhostImporter(vecSpace_, nGhost, ghosts);
  }
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
