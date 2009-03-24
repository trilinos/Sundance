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

#ifndef SUNDANCE_SPATIALLYHOMOGENEOUSDOFMAP_H
#define SUNDANCE_SPATIALLYHOMOGENEOUSDOFMAP_H


#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * 
     */
    class SpatiallyHomogeneousDOFMapBase : public DOFMapBase
    {
    public:
      /** */
      SpatiallyHomogeneousDOFMapBase(const Mesh& mesh, int nTotalFuncs, 
        const ParameterList& verbParams = *DOFMapBase::defaultVerbParams()) ;

      /** */
      virtual ~SpatiallyHomogeneousDOFMapBase() {;}

      /** */
      virtual RefCountPtr<const MapStructure> mapStruct() const = 0 ;

      /** */
      virtual RefCountPtr<const Set<int> >
      allowedFuncsOnCellBatch(int cellDim,
                              const Array<int>& cellLID) const 
      {return allowedFuncs_;}

      /** */
      bool isHomogeneous() const {return true;}

      /** */
      void print(ostream& os) const ;

      /** */
      const Array<CellFilter>& funcDomains() const {return funcDomains_;}

    private:
      RefCountPtr<const Set<int> > allowedFuncs_;
      Array<CellFilter> funcDomains_;
    };
  }
}

#endif
