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

#ifndef SUNDANCE_QUADRATUREEVALMEDIATOR_H
#define SUNDANCE_QUADRATUREEVALMEDIATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceOrderedTuple.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  namespace Internal
  {
    using namespace Teuchos;

    /**
     * 
     */
    class QuadratureEvalMediator : public StdFwkEvalMediator
    {
    public:
      /** 
       * 
       */
      QuadratureEvalMediator(const Mesh& mesh, 
        int cellDim,
        const QuadratureFamily& quad,
        int verb);

      /** */
      virtual ~QuadratureEvalMediator(){;}

      
      /** Evaluate the given coordinate expression, putting
       * its numerical values in the given EvalVector. */
      virtual void evalCoordExpr(const CoordExpr* expr,
                                 RefCountPtr<EvalVector>& vec) const ;
      
      /** Evaluate the given discrete function, putting
       * its numerical values in the given EvalVector. */
      virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                           const Array<MultiIndex>& mi,
                                           Array<RefCountPtr<EvalVector> >& vec) const ;

      /** Evaluate the given cell diameter expression, putting
       * its numerical values in the given EvalVector. */
      virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
                                        RefCountPtr<EvalVector>& vec) const ;


      /** Evaluate the given cell vector expression, putting
       * its numerical values in the given EvalVector. */
      virtual void evalCellVectorExpr(const CellVectorExpr* expr,
				      RefCountPtr<EvalVector>& vec) const ;
            
      /** */
      virtual void setCellType(const CellType& cellType,
                               const CellType& maxCellType) ;

      /** */
      virtual void print(ostream& os) const ;

      /** */
      RefCountPtr<Array<Array<double> > > getRefBasisVals(const BasisFamily& basis, 
                                                   int diffOrder) const ;

      /** */
      RefCountPtr<Array<Array<Array<double> > > > getFacetRefBasisVals(const BasisFamily& basis) const ;

      /** */
      const Array<double>& quadWgts() const 
      {return *(refQuadWeights_.get(cellType()));}

      static double& totalFlops() {static double rtn = 0; return rtn;}

      

      static void addFlops(const double& flops) {totalFlops() += flops;}

      /** */
      int numFacetCases() const {return numFacetCases_;}

    private:


      /** */
      void fillFunctionCache(const DiscreteFunctionData* f,
                             const MultiIndex& mi) const ;

     
      /** */
      void computePhysQuadPts() const ;

      /** */
      int numFacetCases_;

      /** */
      QuadratureFamily quad_;

      /** */
      Map<CellType, RefCountPtr<Array<Point> > > refQuadPts_;

      /** */
      Map<CellType, RefCountPtr<Array<Array<Point> > > > refFacetQuadPts_;

      /** */
      Map<CellType, RefCountPtr<Array<double> > > refQuadWeights_;

      /** */
      mutable Array<Point> physQuadPts_;

      /** */
      mutable Array<Map<OrderedPair<BasisFamily, CellType>, RefCountPtr<Array<Array<double> > > > > refBasisVals_;

      mutable Map<OrderedPair<BasisFamily, CellType>, RefCountPtr<Array<Array<Array<double> > > > > refFacetBasisVals_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
