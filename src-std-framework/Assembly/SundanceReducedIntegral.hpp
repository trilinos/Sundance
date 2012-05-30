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

#ifndef SUNDANCE_REDUCED_INTEGRAL_H
#define SUNDANCE_REDUCED_INTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceElementIntegral.hpp"

namespace Sundance
{

using namespace Teuchos;

/** 
 * 
 */
class ReducedIntegral : public ElementIntegral
{
public:
  /** Construct a zero-form to be computed by reference integration
   * with coefficients that are piecewise constant */
  ReducedIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a one form to be computed by reference integration
   * with coefficients that are piecewise constant  */
  ReducedIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a two-form to be computed by reference integration
   * with coefficients that are piecewise constant */
  ReducedIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim,
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const BasisFamily& unkBasis,
    int beta,
    int unkDerivOrder,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** virtual dtor */
  virtual ~ReducedIntegral(){;}

  /** */
  void transform(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetNum,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeffs,
    RCP<Array<double> >& A) const
    {
      if (order()==2) transformTwoForm(JTrans, JVol, facetNum, cellLIDs, coeffs, A);
      else if (order()==1) transformOneForm(JTrans, JVol, facetNum, cellLIDs, coeffs, A);
      else transformZeroForm(JTrans, JVol, isLocalFlag, facetNum,
        cellLIDs, coeffs, A);
    }

  /** */
  virtual void transformZeroForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeffs,
    RCP<Array<double> >& A) const ;
      
  /** */
  virtual void transformTwoForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeffs,
    RCP<Array<double> >& A) const ;
      
  /** */
  void transformOneForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeffs,
    RCP<Array<double> >& A) const ;

private:

  

  /** */
  inline double& value(int facetCase, int testDerivDir, int testNode,
    int unkDerivDir, int unkNode)
    {return W_[facetCase][unkNode + nNodesUnk()*testNode 
        + nNodes()*(unkDerivDir 
          + nRefDerivUnk()*testDerivDir)];}

  /** */
  inline const double& value(int facetCase, 
    int testDerivDir, int testNode,
    int unkDerivDir, int unkNode) const 
    {
      return W_[facetCase][unkNode + nNodesUnk()*testNode 
        + nNodes()*(unkDerivDir 
          + nRefDerivUnk()*testDerivDir)];
    }
      
  /** */
  inline double& value(int facetCase, int testDerivDir, int testNode)
    {return W_[facetCase][nNodesTest()*testDerivDir + testNode];}

  /** */
  inline const double& value(int facetCase, 
    int testDerivDir, int testNode) const 
    {return W_[facetCase][nNodesTest()*testDerivDir + testNode];}

  static double& totalFlops() {static double rtn = 0; return rtn;}

protected:

  static void addFlops(const double& flops) {totalFlops() += flops;}
      
private:

  Array<Array<double> > W_;

};
}


#endif
