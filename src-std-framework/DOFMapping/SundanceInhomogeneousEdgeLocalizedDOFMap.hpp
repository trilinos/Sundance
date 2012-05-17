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

#ifndef SUNDANCE_INHOMOGENEOUSEDGELOCALIZEDDOFMAP_H
#define SUNDANCE_INHOMOGENEOUSEDGELOCALIZEDDOFMAP_H

#include "SundanceDOFMapBase.hpp"

#include "SundanceDefs.hpp"
#include "SundanceCellFilter.hpp"

#include "SundanceSet.hpp"
#include "SundanceMap.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{

using Teuchos::RCP;
using Teuchos::Array;

/** 
 * 
 */
class InhomogeneousEdgeLocalizedDOFMap : public DOFMapBase
{
public:
  /** */
  InhomogeneousEdgeLocalizedDOFMap(const Mesh& mesh, 
    const Array<Map<Set<int>, CellFilter> >& funcSetToDomainMap, 
    int setupVerb);
  
  /** */
  RCP<const MapStructure> 
  getDOFsForCellBatch(int cellDim,
    const Array<int>& cellLID,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verb) const;
  
  /** */
  RCP<const Set<int> >
  allowedFuncsOnCellBatch(int cellDim,
    const Array<int>& cellLID) const;
  
  /** */
  const Array<CellFilter>& funcDomains() const { return funcDomains_; }
  
  /** */
  virtual void print(std::ostream& os) const ;

private:
  Array<CellFilter> funcDomains_;
  Array<Array<int> > edgeDofs_;

  int meshDimension() const;

  Array<int> getEdgeLIDs(const CellFilter &filter) const;
  
  void getDOFsForEdgeBatch(const Array<int> &cellLID,
    const Set<int> &requestedFuncSet,
    Array<Array<int> > &dofs,
    int verb) const;
  
  RCP<Set<int> > allowedFuncsOnEdgeBatch(const Array<int> &edgeLIDs) const;
  RCP<Set<int> > allFuncIDs() const;
};

} // namespace Sundance

#endif /* SUNDANCE_INHOMOGENEOUSEDGELOCALIZEDDOFMAP_H */
