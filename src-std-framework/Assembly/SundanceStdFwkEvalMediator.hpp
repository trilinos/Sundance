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

#ifndef SUNDANCE_STDFWKEVALMEDIATOR_H
#define SUNDANCE_STDFWKEVALMEDIATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceDiscreteFunction.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceStdFwk

{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using SundanceUtils::Map;

  namespace Internal
  {
    using namespace Teuchos;

    /**

     * StdFwkEvalMediator evaluates mesh-dependent functions in the
     * standard framework. A number of subtypes are supported: 
     * QuadratureEvalMediator, which does evaluation on quadrature points,
     * and NodalEvalMediator, which does evaluation at nodal points. 
     */
    class StdFwkEvalMediator : public AbstractEvalMediator,
                               public TSFExtended::Printable
    {
    public:
      /** */
      StdFwkEvalMediator(const Mesh& mesh, int cellDim)
        : mesh_(mesh),
          cellDim_(cellDim),
          cellType_(NullCell),
          cellLID_(),
          cacheIsValid_(false)
      {;}

      /** */
      virtual ~StdFwkEvalMediator(){;}

      /** */
      void setCellBatch(const RefCountPtr<Array<int> >& cellLID,
                        RefCountPtr<CellJacobianBatch>& J) ;

      /** */

      virtual void setCellType(const CellType& cellType) 
      {cellType_=cellType; cacheIsValid() = false; jCacheIsValid_=false;}


    protected:
      const Mesh& mesh() const {return mesh_;}

      int cellDim() const {return cellDim_;}

      const CellType& cellType() const {return cellType_;}

      const RefCountPtr<Array<int> >& cellLID() const {return cellLID_;}

      bool& cacheIsValid() const {return cacheIsValid_;}

      const RefCountPtr<CellJacobianBatch>& J() const {return J_;}


      /** */
      Map<const DiscreteFunctionData*, RefCountPtr<Array<double> > >& fCache() const {return fCache_;}
      /** */
      Map<const DiscreteFunctionData*, RefCountPtr<Array<double> > >& dfCache() const {return dfCache_;}
      /** */
      Map<const DiscreteFunctionData*, RefCountPtr<Array<double> > >& localValueCache() const {return localValueCache_;}
      /** */
      Map<const DiscreteFunctionData*, bool>& fCacheIsValid() const {return fCacheIsValid_;}
      /** */
      Map<const DiscreteFunctionData*, bool>& dfCacheIsValid() const {return dfCacheIsValid_;}
      /** */
      Map<const DiscreteFunctionData*, bool>& localValueCacheIsValid() const {return localValueCacheIsValid_;}
      
    private:
      Mesh mesh_;

      int cellDim_;

      CellType cellType_;

      RefCountPtr<Array<int> > cellLID_;

      mutable RefCountPtr<CellJacobianBatch> J_;

      mutable bool cacheIsValid_;

      mutable bool jCacheIsValid_;

      /** */
      mutable Map<const DiscreteFunctionData*, RefCountPtr<Array<double> > > fCache_;
      /** */
      mutable Map<const DiscreteFunctionData*, RefCountPtr<Array<double> > > dfCache_; 
      /** */
      mutable Map<const DiscreteFunctionData*, RefCountPtr<Array<double> > > localValueCache_;

      /** */
      mutable Map<const DiscreteFunctionData*, bool> fCacheIsValid_;
      /** */
      mutable Map<const DiscreteFunctionData*, bool> dfCacheIsValid_;
      /** */
      mutable Map<const DiscreteFunctionData*, bool> localValueCacheIsValid_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
