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

#ifndef SUNDANCE_UNIFORM_REFINEMENT_PAIR_H
#define SUNDANCE_UNIFORM_REFINEMENT_PAIR_H

#include "SundanceDefs.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceArrayOfTuples.hpp"

namespace Sundance
{
/**
 *
 */
class UniformRefinementPair
{
public:
  /** */
  UniformRefinementPair(const MeshType& meshType,
    const Mesh& coarse);

  /** */
  const Mesh& fine() const {return fine_;}

  /** */
  const Mesh& coarse() const {return coarse_;}

  /** */
  const Array<int>& oldToNewVertMap() const {return oldToNewVertMap_;}

  /** */
  const Array<int>& newVertToOldLIDMap() const {return newVertToOldLIDMap_;}

  /** */
  const Array<int>& newVertIsOnEdge() const {return newVertIsOnEdge_;}

  /** */
  const Array<int>& oldEdgeToNewVertMap() const {return oldEdgeToNewVertMap_;}

  /** */
  const ArrayOfTuples<int>& oldToNewElemMap() const 
    {return oldToNewElemMap_;}

  /** */
  const Array<int>& newToOldElemMap() const 
    {return newToOldElemMap_;}

  /** */
  const Array<Array<int> >& oldEdgeChildren() const 
    {return oldEdgeChildren_;}

  /** */
  const Array<Array<int> >& oldEdgeParallels() const 
    {return oldEdgeParallels_;}
  /** */
  const Array<int>& newEdgeParents() const 
    {return newEdgeParents_;}

  /** */
  const Array<int>& newEdgeParallels() const 
    {return newEdgeParallels_;}


  /** */
  const ArrayOfTuples<int>& interiorEdgesOfCoarseElems() const 
    {return interiorEdges_;}


  /** Run a consistency check on the pair of meshes. Returns the number
   * of errors detected. */
  int check() const ;

protected:
  void refineTriMesh();

  int lookupEdge(const Mesh& mesh, int v1, int v2) const ;
  
private:
  MeshType meshType_;
  Mesh coarse_;
  Mesh fine_;

  Array<int> oldToNewVertMap_;
  Array<int> newVertIsOnEdge_;
  Array<int> newVertToOldLIDMap_;
  Array<int> oldEdgeToNewVertMap_;

  ArrayOfTuples<int> oldToNewElemMap_;
  Array<int> newToOldElemMap_;

  Array<Array<int> > oldEdgeChildren_;
  Array<Array<int> > oldEdgeParallels_;
  Array<int> newEdgeParents_;
  Array<int> newEdgeParallels_;

  ArrayOfTuples<int> interiorEdges_;
};

}


#endif
