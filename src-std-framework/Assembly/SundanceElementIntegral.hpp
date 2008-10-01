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

#ifndef SUNDANCE_ELEMENTINTEGRAL_H
#define SUNDANCE_ELEMENTINTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "Teuchos_Array.hpp"

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
 * ElementIntegral encapsulates the common data needed for the
 * integration of groups of related zero-forms, one-forms, and two-forms.
 */
class ElementIntegral 
  : public TSFExtended::ParameterControlledObjectWithVerbosity<ElementIntegral>,
    public TSFExtended::Printable
{
public:
  /** Construct a zero-form */
  ElementIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const ParameterList& verbParams);

  /** Construct a one-form */
  ElementIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const ParameterList& verbParams);

  /** Construct a two-form */
  ElementIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim,
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const BasisFamily& unkBasis,
    int beta,
    int unkDerivOrder,
    const ParameterList& verbParams);

  /** virtual dtor */
  virtual ~ElementIntegral(){;}

  /** Indicate whether this element integral is a 
   *  zero, one, or two form */
  int order() const {return order_;}
      
  /** Return the number of nodes associated with the test function */
  int nNodesTest() const {return nNodesTest_;}
      
  /** Return the number of nodes associated with the test function */
  int nNodesUnk() const {return nNodesUnk_;}

  /** Return the total number of elements in this local stiffness
   * matrix */
  int nNodes() const {return nNodes_;}

  /** Return the number of different facets for which integrals must be tabulated
   * in the cases where an integral must be done by referring back to a maximal
   * cell */
  int nFacetCases() const {return nFacetCases_;}


  /** */
  static int& transformationMatrixIsValid(int alpha);


  /** */
  static int& transformationMatrixIsValid(int alpha, int beta);

  /** */
  static void invalidateTransformationMatrices();

  /** */
  static double& totalFlops() {static double rtn = 0; return rtn;}


  /** */
  static RefCountPtr<ParameterList> defaultVerbParams()
    {
      static RefCountPtr<ParameterList> rtn = rcp(new ParameterList("Integration"));
      static int first = true;
      if (first)
      {
        rtn->set<int>("setup", 0);
        rtn->set<int>("transformation", 0);
        rtn->set<int>("integration", 0);
        rtn->set<int>("extract weak form", 0);
        rtn->set<int>("find groups", 0);
        first = false;
      }
      return rtn;
    }

protected:
  /** */
  static void addFlops(const double& flops) {totalFlops() += flops;}

  /** The dimension of the cell being integrated */
  int dim() const {return dim_;}

  /** The dimension of the space in which the cell is embedded */
  int spatialDim() const {return spatialDim_;}
      
  /** Number of test function derivatives wrt reference coordinates that
   * are needed to evaluate this integral. Will always be equal to 
   * ipow(element dimension, differentiation order). */
  int nRefDerivTest() const {return nRefDerivTest_;}
      
  /** Number of unknown function derivatives wrt reference coordinates that
   * are needed to evaluate this integral. Will always be equal to 
   * ipow(element dimension, differentiation order). */
  int nRefDerivUnk() const {return nRefDerivUnk_;}

  /** The order to which the test function is differentiated in this
   * integral. */
  int testDerivOrder() const {return testDerivOrder_;}

  /** The order to which the unknown function is differentiated in this
   * integral. */
  int unkDerivOrder() const {return unkDerivOrder_;}

  /** */
  int alpha() const {return alpha_;}

  /** */
  int beta() const {return beta_;}

  /** Workspace for element transformations involving one derivative*/
  static Array<double>& G(int gamma) ;

  /** Workspace for element transformations involving two derivatives */
  static Array<double>& G(int gamma, int delta) ;

  /** return base to the given power */
  static int ipow(int base, int power);

  /** The value below which chop() sets numbers to zero */
  static double chopVal() {static double rtn=1.0e-14; return rtn;}

  /** Chop a number to zero if it is smaller in magnitude than
   * the value chopVal() */
  static double chop(const double& x) 
    {
      if (::fabs(x) > chopVal()) return x;
      else return 0.0;
    }


  /** */
  void createTwoFormTransformationMatrix(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol) const;
  /** */
  void createOneFormTransformationMatrix(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol) const;

private:

  int spatialDim_;

  int dim_;

  int nFacetCases_;

  int testDerivOrder_;

  int nRefDerivTest_;

  int nNodesTest_;

  int unkDerivOrder_;

  int nRefDerivUnk_;

  int nNodesUnk_;

  int nNodes_;

  int order_;

  int alpha_;

  int beta_;
      
};
}
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
