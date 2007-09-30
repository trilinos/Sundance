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


#ifndef SUNDANCE_FIATLAGRANGE_H
#define SUNDANCE_FIATLAGRANGE_H

#include "SundanceDefs.hpp"
#include "SundanceBasisFamilyBase.hpp"
#include "SundanceMap.hpp"

#ifdef HAVE_FIAT
#include "FIAT.hpp"
#endif

namespace SundanceStdFwk
{
  /** 
   * Lagrange basis, implemented through the FIAT element library
   */
  class FIATLagrange : public ScalarBasis
  {
  public:
    /** */
    FIATLagrange(int order);

    /**     
     * \brief   Inform caller as to whether a given cell type is supported 
     */
    bool supportsCellTypePair(
      const CellType& maximalCellType,
      const CellType& cellType
      ) const ;

    /** */
    virtual void print(std::ostream& os) const ;

    /** */
    virtual int order() const {return order_;}

    /** return the number of nodes for this basis on the given cell type */
    virtual int nReferenceDOFs(
      const CellType& maximalCellType,
      const CellType& cellType
      ) const ;

    /** */
    virtual void getReferenceDOFs(
      const CellType& maximalCellType,
      const CellType& cellType,
      Array<Array<Array<int> > >& dofs) const ;

    /** */
    virtual void refEval(
      const CellType& maximalCellType,
      const CellType& cellType,
      const Array<Point>& pts,
      const MultiIndex& deriv,
      Array<Array<Array<double> > >& result) const ;




    /* Handleable boilerplate */
    GET_RCP(BasisFamilyBase);

  private:
    int order_;
#ifdef HAVE_FIAT
    // RCK
    // this has to be mutable because lookups on maps
    // are non-const :/ and otherwise seemingly const methods
    // like tabulation and nNodes die at compile-time
    mutable SundanceUtils::Map<CellType,FIAT::RCP<FIAT::LagrangeElement> > fiatElems_;
#endif
  };


}

#endif

